import bz2
import gzip
import logging

import numba as nb
import numpy as np
import os
import simplejson as json
import subprocess

from hashlib import sha256

from json2html import json2html
from typing import Optional, TextIO, Dict
from typing.re import Pattern as tPattern
from mimetypes import guess_type
from .exceptions import ApplicationException

logger = logging.getLogger(__name__)


def simple_observed_fraction(obs_extent, mean_frag_size, n_fragments):
    """
    A simple estimate of the observed fraction seen from the reads over a set of analyzed
    fragments. Here, the total fragment extent is estimated by N_frag * mean(frag). The user
    is expected to computed the total extent of all observed reads.

    The mean fragment size should be empirically deteremined from data, rather than rely on
    a value quoted by a sequencing facility. These are often wrong. Software that assumes
    shotgun data will also do poorly at estimating the mean fragment size when inspecting
    Hi-C reads.

    :param obs_extent: total extent of reads analyzed - Sum(len(read_i))
    :param mean_frag_size: the mean fragment size
    :param n_fragments: the number of unique fragments analyzed
    :return: the observed fraction
    """
    return obs_extent / (mean_frag_size * n_fragments)


@nb.jit(nopython=True, cache=True)
def make_observable_mask(read_len: int, insert_len: int, paired: bool,
                         left_margin: int, right_margin: int) -> np.ndarray:
    """
    Use a mask to calculate the proportion of a fragment which can be interrogated
    due to the need for a sliding window around any observed junction sequence.

    :param read_len: the average read-length
    :param insert_len: the average insert (or fragment) length
    :param paired: True=fragment observed from both ends, otherwise one end obs only
    :param left_margin: left-side margin which cannot be observed due to algorithm constraints
    :param right_margin: right-side margin which cannot observed due to algorithm constraints
    :return: mask
    """
    frag_mask = np.zeros(insert_len, np.uint8)
    read_mask = np.zeros(read_len, np.uint8)
    # create a read mask that represents the region of the read which can be interrogated
    x_min, x_max = left_margin, read_len - right_margin
    read_mask[x_min:x_max] = 1
    # create a fragment mask by transferring this silhouette to either end of the fragment
    # handling the edge-case where the insert length is less than the read length
    a = read_len if read_len < insert_len else insert_len
    frag_mask[:a] += read_mask[:a]
    if paired:
        frag_mask[-a:] += read_mask[::-1][-a:]
    # return the fragment mask
    return frag_mask


def observed_fraction(read_len: int, insert_len: int, method: str, paired: bool,
                      left_margin: int = 0, right_margin: int = 0) -> float:
    """
    Calculate an estimate of the observed fraction. Here, read-pairs provide a means of inspecting
    the sequenced fragments for Hi-C junctions. Additionally, the k-mer and junction size affect
    how much of each read can be inspected.

    :param read_len: the average read-length
    :param insert_len: the average insert (or fragment) length
    :param method: "additive" or "binary" determines how the mean of the mask is calculated.
    :param paired: True=fragment observed from both ends, otherwise one end obs only
    :param left_margin: left-side margin which cannot be observed due to algorithm constraints
    :param right_margin: right-side margin which cannot observed due to algorithm constraints
    :return: estimated fraction of extent observed depending on method
    """
    if method == 'additive':
        return make_observable_mask(read_len, insert_len, paired, left_margin, right_margin).mean()
    elif method == 'binary':
        return (make_observable_mask(read_len, insert_len, paired, left_margin, right_margin) > 0).mean()
    else:
        raise ApplicationException('unknown method {}'.format(method))


def clean_output(msg: bytes):
    """
    Convert bytes to string and strip ends.
    :param msg: the bytes message to decode and strip
    :return: str version of message.
    """
    if msg is None:
        return ''
    return msg.decode().strip()


def count_sequences(file_name: str, fmt: str, max_cpu: int = 1, wait_timeout: float = None) -> int:
    """
    Estimate the number of fasta or fastq sequences in a file by counting headers. Decompression
    is automatically attempted for files ending in .gz. Counting and decompression is by
    way of subprocess calls to grep and gzip. Uncompressed files are also handled. This
    is about 8 times faster than parsing a file with BioPython and 6 times faster than
    reading all lines in Python.
    :param file_name: the file to inspect
    :param fmt: fasta or fastq
    :param max_cpu: the number of cpus if pigz exists
    :return: the estimated number of records
    """

    if fmt == 'fasta':
        pattern = r'^>'
    elif fmt == 'fastq':
        pattern = r'^@'
    else:
        raise ApplicationException('sequence format {} is unknown'.format(fmt))

    if not os.path.exists(file_name):
        raise FileNotFoundError(file_name)

    if max_cpu < 1:
        raise ApplicationException('max_cpus must 1 or greater')

    if wait_timeout is not None and wait_timeout <= 0:
        raise ApplicationException('wait_timeout must be None or greater than 0')

    if file_name.endswith('.gz'):
        pigz_path = check_executable_exists('pigz')
        if max_cpu > 1 and pigz_path is not None:
            proc_uncomp = subprocess.Popen([pigz_path, '-p{}'.format(max_cpu), '-cd', file_name],
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name],
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_grep = subprocess.Popen(['grep', '-c', pattern], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_grep = subprocess.Popen(['grep', '-c', pattern, file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    # wait for the process to complete and then pull the results in with communicate
    retval = proc_grep.wait(wait_timeout)
    outs, errs = proc_grep.communicate()

    if retval > 1:
        # grep returns 1 when there are no matching lines, therefore only raise
        # errors on values larger.
        raise OSError('retval=[{}] stdout=[{}] stderr=[{}]'.format(
            retval, clean_output(outs), clean_output(errs)))

    return int(outs)


def modification_hash(file_path: str) -> str:
    """
    Create a string hash for a given file using base file name,
    modification time and size as a proxy for file changes.

    :param file_path: the path to the target file
    :return: string hash
    """
    stat = os.stat(file_path)
    hasher = sha256()
    hasher.update(os.path.basename(file_path).encode())
    hasher.update(stat.st_size.to_bytes(16, 'little', signed=False))
    hasher.update(hex(stat.st_mtime_ns).encode())
    return hasher.hexdigest()


def open_input(file_name: str) -> TextIO:
    """
    Open a file, transparently handling compressed files, regardless
    of suffix.

    :param file_name: the file path to open
    :return: a TextIO object
    """
    if not os.path.exists(file_name):
        raise FileNotFoundError(file_name)
    _encoding = guess_type(file_name)[1]
    if _encoding == 'gzip':
        return gzip.open(file_name, 'rt')
    if _encoding == 'bzip2':
        return bz2.open(file_name, 'rt')
    else:
        return open(file_name, 'rt')


def warn_if(_test: bool) -> int:
    """
    Conditional help for logging warnings.
    :param _test: if true, level is that of WARNING
    :return: return WARNING if true, otherwise INFO
    """
    if _test:
        return logging.WARNING
    else:
        return logging.INFO


def check_executable_exists(prog_name: str) -> Optional[str]:
    """
    Test whether a program exists on the system. This can be either a full path
    or just the executable name. This is case sensitive.
    :param prog_name: the program name with/without path
    :return: the full path or None
    """

    def is_exe(fp):
        return os.path.isfile(fp) and os.access(fp, os.X_OK)

    fpath, fname = os.path.split(prog_name)
    if fpath:
        if is_exe(prog_name):
            return prog_name
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, prog_name)
            if is_exe(exe_file):
                return exe_file

    return None


def init_random_state(seed: int = None) -> np.random.RandomState:
    """
    Helper function to initialise a random state using a given seed, or
    from a random seed. The used seed value is sent to log.
    :param seed: an integer or None. If none, a randomly selected seed is used.
    :return: an initialised numpy.random.RandomState
    """
    # set up random number generation
    if seed is None:
        seed = np.random.randint(1000000, 99999999)
        logger.info('Random seed was not set, using {}'.format(seed))
    else:
        logger.info('Random seed was {}'.format(seed))
    return np.random.RandomState(seed)


def to_json_string(obj) -> str:
    """
    Convert an object to a JSON string
    :param obj: the object to convert
    :return: a JSON format string
    """
    def qc3c_default(obj):
        if isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.ndarray):
            return [qc3c_default(oi) for oi in obj]
        elif isinstance(obj, tPattern):
            return obj.pattern
        raise TypeError(repr(obj) + " is not JSON serializable")

    return json.dumps(obj, default=qc3c_default)


def write_jsonline(fpath: str, obj, suffix='jsonl', append=True):
    """
    Append an object to a JSON Lines format file. (newline delimited json). The suffix "jsonl" will be
    appended if not included
    :param fpath: the file path to open for writing
    :param obj: the object to write
    :param suffix: file name suffix
    :param append: True append to existing files, False overwrite
    """
    if not fpath.endswith(suffix):
        fpath = '{}.{}'.format(fpath, suffix)
    _mode = 'at+' if append else 'at'
    with open(fpath, _mode) as fp:
        print(to_json_string(obj), file=fp)


def write_html_report(fpath: str, report: Dict):
    """
    Write with truncation an HTML format version of the qc3C report. The suffix "html" will be appended
    if not included.
    :param fpath: the file path to open for writing
    :param report: a qc3C report dictionary
    """
    if not fpath.endswith(".html"):
        fpath = '{}.html'.format(fpath)
    with open(fpath, 'wt') as fp:
        fp.write("""
            <!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN">
            <html lang="en">
            <head>
                    <title>qc3C report</title>
                    <style>
                        table {border: 1px solid black; border-collapse: collapse;}
                        td {padding: 5px; border: 1px solid black;}
                        th {padding: 5px; border: 1px solid black; background-color: lightblue;}
                        ul {list-style-type: none; padding: 0; margin: 0;}
                        li {padding: 0 0 5px 0;}
                    </style>
                </head>
                <body>
        """)
        _html_table = json2html.convert(to_json_string(report), table_attributes='', escape=True)
        fp.write('        {}\n'.format(_html_table))
        fp.write("    </body>\n</html>")


def guess_quality_encoding(fastq_path: str, n_reads: int = 5000, solexa: bool = False) -> int:
    """
    Look at a sample of reads within a FastQ format file, to guess
    whether its encoding begins at 33 or 64. The more reads inspected
    the greater the odds that one will contain low enough qualities
    to demonstrate that it is base-33.
    :param fastq_path: the path to the FastQ format file
    :param n_reads: the number of reads to inspect
    :param solexa: if true, report when its possibly solexa encoding
    :return: 33 if base-33, 64 if base-64 (or 59 if solexa = True)
    """
    with open_input(fastq_path) as in_h:
        lowest = 999
        for n, (_id, _desc, _seq, _qual) in enumerate(read_seq(in_h), 1):
            min_qual = np.asarray(bytearray(_qual.encode())).min()
            if min_qual < lowest:
                lowest = min_qual
            if n == n_reads:
                break

    # Phred33
    if lowest < 59:
        enc = 33
    # possibly Solexa
    elif solexa and lowest < 64:
        enc = 59
    # Phred64
    else:
        enc = 64
    return enc


def read_seq(fp: TextIO) -> (str, Optional[str], str, Optional[str]):
    """
    Method to quickly read FastA or FastQ files using a generator function.
    Originally sourced from https://github.com/lh3/readfq
    :param fp: input file object
    :return: tuple (name, desc, seq, qual)
    """
    last = None  # this is a buffer keeping the last unprocessed line

    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break

        parts = last[1:].partition(" ")
        name, desc, seqs, last = parts[0], parts[2] if len(parts) > 1 else None, [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield name, desc, ''.join(seqs), None  # yield a fasta record
            if not last:
                break

        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, desc, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, desc, seq, None  # yield a fasta record instead
                break
