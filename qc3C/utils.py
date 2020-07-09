import bz2
import gzip
import itertools
import logging
from collections import Counter as collectionsCounter

import numpy as np
import os
import simplejson as json
import subprocess

from hashlib import sha256

import tqdm
from json2html import json2html
from scipy.stats import beta
from typing import Optional, TextIO, Dict, List
from typing.re import Pattern as tPattern
from mimetypes import guess_type
from .exceptions import ApplicationException
from .kmer_based import choice

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


def observed_fraction(read_len: int, insert_len: int, kmer_size: int = 0,
                      junc_size: int = 0, is_phase: bool = False, is_single: bool = False) -> float:
    """
    Calculate an estimate of the observed fraction. Here, read-pairs provide a means of inspecting
    the sequenced fragments for Hi-C junctions. Additionally, the k-mer and junction size affect
    how much of each read can be inspected.

    :param read_len: the average read-length
    :param insert_len: the average insert (or fragment) length
    :param kmer_size: the requested k-mer size
    :param junc_size: the junction size
    :param is_phase: whether or not the library was constructed using Phase's kit
    :param is_single: the observational set has been reduced to a single read per fragment
    :return: the observed fraction [0, 1]
    """

    def make_observable_mask(read_len: int, insert_len: int, kmer_size: int,
                             junc_size: int, is_single: bool) -> np.ndarray:
        """
        Use a mask to calculate the proportion of a fragment which can be interrogated
        due to the need for a sliding window around any observed junction sequence.

        :param read_len: the average read-length
        :param insert_len: the average insert (or fragment) length
        :param kmer_size: the requested k-mer size
        :param junc_size: the junction size
        :param is_single: the observational set has been reduced to a single read per fragment
        :return:
        """

        frag_mask = np.zeros(insert_len, dtype=np.uint8)
        read_mask = np.zeros(read_len, dtype=np.uint8)
        x_min, x_max = kmer_size + 1, read_len - (kmer_size + junc_size + 1)
        # create a read mask that represents the region of the read which can be interrogated
        read_mask[x_min:x_max] = 1
        # create a fragment mask by transferring this silhouette to either end of the fragment
        # handling the edge-case where the insert length is less than the read length
        a = read_len if read_len < insert_len else insert_len
        frag_mask[:a] += read_mask[:a]
        if not is_single:
            # treat both ends if paired
            frag_mask[-a:] += read_mask[::-1][-a:]
        # return the fragment mask
        return frag_mask

    read_len = round(read_len)
    insert_len = round(insert_len)

    obs_mask = make_observable_mask(read_len, insert_len, kmer_size, junc_size, is_single)

    if is_phase:

        logger.debug('Calculating the observed fraction as a Phase library')

        # phase's library construction appears to produce a non-uniform distribution
        # of junction locations, preferring the interior of the fragment to the ends.
        rv = beta(2, 2)

        # find the windows across the fragment which are observed
        obs_wins = []
        start = None
        for i in range(len(obs_mask)):
            if obs_mask[i]:
                if start is None:
                    start = i
            elif start is not None:
                obs_wins.append([start, i])
                start = None
        if start is not None:
            obs_wins.append([start, len(obs_mask)])

        # sum the distribution slices from each window
        obs_frac = 0
        for wi in obs_wins:
            obs_frac += rv.cdf(wi[1] / insert_len) - rv.cdf(wi[0] / insert_len)

    else:
        logger.debug('Calculating the observed fraction as a generic library')

        # uniform treatment
        obs_frac = obs_mask.mean()

    # sanity checks that should not get invoked
    if obs_frac > 1:
        logger.warning('Observed fraction > 1, resetting to 1.')
        obs_frac = 1
    elif obs_frac < 0:
        logger.warning('Observed fraction < 0, resetting to 0.')
        obs_frac = 0

    return obs_frac


def count_sequences(file_name: str, fmt: str, max_cpu: int = 1) -> int:
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

    if file_name.endswith('.gz'):
        pigz_path = test_for_exe('pigz')
        if max_cpu > 1 and pigz_path is not None:
            proc_uncomp = subprocess.Popen([pigz_path, '-p{}'.format(max_cpu), '-cd', file_name],
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name],
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', pattern], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', pattern, file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n


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


def test_for_exe(prog_name: str) -> Optional[str]:
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
        fpath = f'{fpath}.{suffix}'
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
        fpath = f'{fpath}.html'
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


class ReadMarkovModeller:
    """
    Use a k-th order Markov model to compute the expected number of junction
    sequences that would occur in a set of sequences with the given composition
    """

    def __init__(self, kmer_size: int, seed: int, reject_prob: float, targets: List[str]):
        """
        :param kmer_size:
        :param seed:
        :param reject_prob:
        :param targets:
        """
        self.kmer_size = kmer_size
        self.counts = collectionsCounter()
        self.counts_plus_one = collectionsCounter()
        self.trans_prob = None
        self.kmers_prob = None
        self.kmers = None
        self.reject_prob = reject_prob
        self.targets = [ti.encode() for ti in targets]
        self.random_state = np.random.RandomState(seed)
        self.total_obs = 0
        self.total_extent = 0

    def average_read_length(self) -> int:
        """
        :return: Average read length accumulated by the model instance
        """
        return self.total_extent // self.total_obs

    @staticmethod
    def filter_n(kmer_dict: collectionsCounter) -> Dict[bytes, int]:
        """
        Remove records whose k-mer contains one or more Ns
        """
        return {k: v for k, v in kmer_dict.items() if b'N' not in k}

    def contains_target(self, s):
        return any(ti in s for ti in self.targets)

    def prepare_model(self):
        """
        Once some training data has been added, this method will prepare the transition
        probabilty data structure for use in simulating reads.
        """
        # filter k-mers containing ambiguous bases
        clean_counts = ReadMarkovModeller.filter_n(self.counts)

        # rather than rely on kmers from counts, we will assume there
        # are missing observations and instead explicitly iterate
        # over all possibilities.
        self.kmers = [bytes(ki) for ki in itertools.product(b'ACGT', repeat=self.kmer_size)]
        # assign probabilities for each kmer
        self.kmers_prob = [clean_counts[kmer] if kmer in clean_counts else 1 for kmer in self.kmers]
        # normalise these counts to probs
        self.kmers_prob = np.array(self.kmers_prob, dtype=np.float)
        self.kmers_prob /= self.kmers_prob.sum()

        uniform_prob = np.array([0.25] * 4)
        self.trans_prob = {}
        for kmer in self.kmers:
            if kmer not in clean_counts:
                # if we didn't see the kmer, assume uniform probability
                self.trans_prob[kmer] = uniform_prob
            else:
                # generate the four posssible k+1-mers
                plus_ones = [kmer + nt for nt in [b'A', b'C', b'G', b'T']]
                # get the counts for these possibilities
                next_counts = np.array([self.counts_plus_one[po] if po else 0 for po in plus_ones])
                # transition probabilities for this k-mer to the next 4 possibilities
                self.trans_prob[kmer] = next_counts / next_counts.sum()

    def find_target_mid(self, seq: bytes) -> Optional[int]:
        """
        Find target sequences within a read and return the modpoint of the
        matched region. Only the first match is found.

        :param seq: the sequence to search
        """
        self.random_state.shuffle(self.targets)
        for ti in self.targets:
            ix = seq.find(ti)
            if ix != -1:
                return ix + len(ti) // 2
        return None

    def __iadd__(self, other: str):
        """
        Introduce another sequence for training model. In the case of reads which contain
        a target sequence, split the reads at the target sequence midpoint prior to extracting
        k-mers from each end. This aims to suppress k-mer bias in high quality libraries which
        contain a high proportion of target sequences.

        :param other: (sequence is assumed to already be upper case)
        """
        other = other.upper().encode()

        # split reads containing junction at midpoint at given rate
        # this surpresses Hi-C signal bias which is a problem for high quality libraries
        if self.contains_target(other) and self.random_state.uniform() < self.reject_prob:
            break_point = self.find_target_mid(other)
            for oi in [other[:break_point], other[break_point:]]:
                # extract k and k+1 mers
                self.counts.update(oi[i: i+self.kmer_size] for i in range(len(oi) - self.kmer_size + 1))
                self.counts_plus_one.update(oi[i: i+self.kmer_size+1] for i in range(len(oi) - self.kmer_size))
        else:
            # extract k and k+1 mers
            self.counts.update(other[i: i+self.kmer_size] for i in range(len(other) - self.kmer_size + 1))
            self.counts_plus_one.update(other[i: i+self.kmer_size+1] for i in range(len(other) - self.kmer_size))

        self.total_obs += 1
        self.total_extent += len(other)
        return self

    def estimate_freq(self, sample_size: int) -> float:
        """
        Compute the expected frequency at which any of the target sequences appear in reads
        simulated from the model.

        :param sample_size:
        :return:
        """

        read_len = int(self.average_read_length())

        # simulate sequences at the target read length and count the fraction
        # that contain the junction sequence(s)
        acgt = bytearray(b'ACGT')
        # convert the list of junctions to a set of bytes
        containing_reads = 0
        # generate an initial starting point for each read
        initial_kmers = self.random_state.choice(self.kmers, size=sample_size, p=self.kmers_prob)
        for kmer in tqdm.tqdm(initial_kmers):
            # choose a kmer based on observed prob
            # start the simulated sequence from this
            sim_seq = bytearray(kmer + (b'N' * (read_len - self.kmer_size)))
            for i in range(self.kmer_size, len(sim_seq)):
                sim_seq[i] = choice(acgt, self.trans_prob[bytes(sim_seq[i-self.kmer_size: i])], self.random_state)
            # does our simulated sequence contain the junction?
            if self.contains_target(sim_seq):
                containing_reads += 1

        return containing_reads / sample_size