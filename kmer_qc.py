import numpy as np
import dna_jellyfish
import tqdm
import pandas
import subprocess
import bz2
import gzip
from collections import namedtuple

CovInfo = namedtuple('cov_info', ['mean_inner', 'mean_outer', 'read_type', 'read_id', 'ix'])


def collect_coverage(seq, ix, site_size, k):
    assert k <= ix <= len(seq) - (k + site_size), \
        'The site index {} is either too close to start (min {}) or ' \
        'end (max {}) of read to scan for coverage'.format(ix, k, len(seq) - (k + site_size))
    dat = []
    for i in range(-k, site_size+1):
        smer = seq[ix+i:ix+i+k]
        mer = dna_jellyfish.MerDNA(smer)
        mer.canonicalize()
        dat.append(query_jf[mer])
    return np.array(dat)


def open_input(file_name):
    """
    Open a text file for input. The filename is used to indicate if it has been
    compressed. Recognising gzip and bz2.

    :param file_name: the name of the input file
    :return: open file handle, possibly wrapped in a decompressor
    """
    suffix = file_name.split('.')[-1].lower()
    if suffix == 'bz2':
        return bz2.open(file_name, 'rt')
    elif suffix == 'gz':
        return gzip.open(file_name, 'rt')
    else:
        return open(file_name, 'rt')


def readfq(fp):  # this is a generator function

    last = None  # this is a buffer keeping the last unprocessed line

    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break

        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])

        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break

        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def next_read(reader, site, k_size):

    site_size = len(site)

    for _id, _seq, _qual in reader:
        seq_len = len(_seq)
        if seq_len < 2*k_size + site_size + 1:
            continue
        _seq = _seq.upper()
        _ix = _seq.find(site)
        # no site found
        if _ix == -1:
            yield(_seq, None, _id, seq_len)
        # only report sites which meet flank constraints
        elif k_size <= _ix <= (seq_len - (k_size + site_size)):
            yield(_seq, _ix, _id, seq_len)


def count_fastq_sequences(file_name):
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression is automatically attempted
    for files ending in .gz. Counting and decompression is by why of subprocess calls to grep and gzip. Uncompressed
    files are also handled. This is about 8 times faster than parsing a file with BioPython and 6 times faster
    than reading all lines in Python.

    :param file_name: the fasta file to inspect
    :return: the estimated number of records
    """
    if file_name.endswith('.gz'):
        proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', r'^@'], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', r'^@', file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser('Prototype kmer analyzer')
    parser.add_argument('--max-n', default=500000, type=int, help='Stop after collecting N samples')
    parser.add_argument('KMER_SIZE', type=int, help='Kmer size used in database')
    parser.add_argument('SITE', help='Junction site')
    parser.add_argument('FASTQ', help='FastQ file used in making the kmer database')
    parser.add_argument('KMER_DB', help='Jellyfish kmer database')
    parser.add_argument('OUTPUT', help='Output table name')
    args = parser.parse_args()

    max_reads = int(args.max_n)
    k_size = args.KMER_SIZE
    site = args.SITE.upper()
    flank_size = (k_size - 8) // 2
    site_size = len(site)

    # initialize jellyfish API
    dna_jellyfish.MerDNA_k(k_size)
    query_jf = dna_jellyfish.QueryMerFile(args.KMER_DB)

    # When finding junction site in reads, the flanks must be big enough to
    # accommodate the sliding window.
    # kmer_patt = re.compile('[^N]{{{0},}}({1})[^N]{{{0},}}'.format(k_size, site))

    OUTER_IX = np.array([True] * (flank_size+2) + [False] * (site_size*2 + 1 - 4) + [True] * (flank_size+2))
    INNER_IX = ~ OUTER_IX

    print('Counting FastQ reads...')
    n_reads = count_fastq_sequences(args.FASTQ)
    print('Analyzing {} reads'.format(n_reads))

    # probability of acceptance for subsampling
    p_accept = 1 if max_reads >= n_reads else max_reads / n_reads * 10
    print('Acceptance threshold: {:.2e}'.format(p_accept))

    unif = np.random.uniform

    with tqdm.tqdm(total=max_reads) as progress:

        n = 0
        cov_dat = []
        fq_reader = next_read(readfq(open_input(args.FASTQ)), site, k_size)
        while True:
            try:
                seq, ix, _id, seq_len = next(fq_reader)
            except StopIteration:
                progress.close()
                print('Reached end of FastQ file')
                break

            # subsample over entire read set to fulfill quota
            if unif() > p_accept:
                continue

            # TODO temp fudge as it seems we have a mixture of MluCI and Sau3AI
            if 'AATTAATT' in seq:
                continue

            # if the read contains no junction, we might still use it
            # as an example of a shotgun read (wgs)
            if ix is None:

                # as there are so many non-junction reads, we need to subsample
                if unif() > 0.1:
                    continue

                # try to find a randomly selected region which does not contain
                # ambiguous bases. Skip the read if we fail in a few attempts
                _attempts = 0
                while _attempts < 3:
                    ix = np.random.randint(k_size, seq_len - (k_size + site_size)+1)
                    subseq = seq[ix - k_size: ix + k_size + site_size]
                    # if no N within the subsequence, then accept it

                    if 'N' not in subseq:
                        break
                    # otherwise keep trying
                    _attempts += 1
                # too many tries occurred, abandon this sequence
                if _attempts >= 3:
                    continue

                rtype = 'wgs'

            # junction containing reads, categorised as hic for simplicity
            else:

                # abandon this sequence if it contains an N
                if 'N' in seq[ix - k_size: ix + k_size + site_size]:
                    continue

                rtype = 'hic'

            cov = collect_coverage(seq, ix, site_size, k_size)

            cov_dat.append(CovInfo(np.mean(cov[INNER_IX]), np.mean(cov[OUTER_IX]), rtype, _id, ix))
            progress.update()

            if args.max_n is not None and len(cov_dat) == args.max_n:
                break

    df = pandas.DataFrame(cov_dat)

    # remove any row which had zero coverage in inner or outer region
    z_inner = df.mean_inner == 0
    z_outer = df.mean_outer == 0
    nz_either = ~(z_inner | z_outer)
    df = df[nz_either]
    print('Rows removed with coverage: inner {}, outer {}, shared {}'.format(
        sum(z_inner), sum(z_outer), sum(z_inner & z_outer)))

    df['ratio'] = df.mean_inner / df.mean_outer
    df.sort_values('ratio', inplace=True)

    print('Initially collected:')
    print(df.groupby('read_type').size())

    wgs = df.loc[df.read_type == 'wgs'].copy()
    hic = df.loc[df.read_type == 'hic'].copy()
    hic['pvalue'] = None

    wgs.reset_index(inplace=True, drop=True)
    wgs['pvalue'] = 1.0 / len(wgs) * (wgs.index + 1)

    df = wgs.append(hic)
    df.sort_values('ratio', inplace=True)
    df.reset_index(inplace=True, drop=True)
    df.to_csv(args.OUTPUT, sep='\t')
