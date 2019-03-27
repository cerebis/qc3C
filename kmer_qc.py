#!/usr/bin/env python3
import numpy as np
import tqdm
import pandas
import subprocess
import bz2
import gzip
from collections import namedtuple

try:
    import dna_jellyfish
except ImportError:
    import jellyfish as dna_jellyfish

CovInfo = namedtuple('cov_info', ['mean_inner', 'mean_outer', 'read_type']) #, 'read_id', 'ix'])


def collect_coverage(seq, ix, site_size, k):
    """
    Collect the k-mer coverage centered around the position ix. From the left, the sliding
    window begins just before the site region and slides right until just after. Means
    are then calculated for the inner (within the junction) and outer (left and right flanks)
    :param seq: the sequence to analyze
    :param ix: the position marking the beginning of a junction or any other arbitrary location if so desired
    :param site_size: the size of the junction site
    :param k: the kmer size
    :return: mean(inner), mean(outer)
    """
    assert k <= ix <= len(seq) - (k + site_size), \
        'The site index {} is either too close to start (min {}) or ' \
        'end (max {}) of read to scan for coverage'.format(ix, k, len(seq) - (k + site_size))

    sliding_cov = np.zeros(shape=k+site_size+1, dtype=np.int)
    for i in range(-k, site_size+1):
        smer = seq[ix+i:ix+i+k]
        mer = dna_jellyfish.MerDNA(smer)
        mer.canonicalize()
        sliding_cov[i] = query_jf[mer]
    return np.mean(sliding_cov[INNER_IX]), np.mean(sliding_cov[OUTER_IX])


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
    parser.add_argument('-s', '--seed', type=int, help='Random seed used in sampling the read-set')
    parser.add_argument('--max-n', default=500000, type=int, help='Stop after collecting N samples')
    parser.add_argument('--max-coverage', default=500, type=int, help='Ignore regions with more than this coverage')
    parser.add_argument('--p-value', default=0.05, type=float, help='p-value threshold for Hi-C junctions')
    parser.add_argument('-o', '--output', help='Output the table of observations to a file')
    parser.add_argument('KMER_SIZE', type=int, help='Kmer size used in database')
    parser.add_argument('SITE', help='Junction site')
    parser.add_argument('FASTQ', help='FastQ file used in making the kmer database')
    parser.add_argument('KMER_DB', help='Jellyfish kmer database')
    args = parser.parse_args()

    max_reads = int(args.max_n)

    # We will treat all sequence in upper case
    site = args.SITE.upper()
    site_size = len(site)

    k_size = args.KMER_SIZE
    # TODO for flexible flanks could be handled if L/R flanks treated independently
    flank_size = (k_size - site_size) // 2

    # some sanity checks.
    assert k_size - site_size > 0, 'Kmer size must be larger the the junction size'
    assert (k_size - site_size) % 2 == 0, 'Kmer size and junction size should match (even/even) or (odd/odd)'

    # initialize jellyfish API
    dna_jellyfish.MerDNA_k(k_size)
    query_jf = dna_jellyfish.QueryMerFile(args.KMER_DB)

    # When finding junction site in reads, the flanks must be big enough to
    # accommodate the sliding window.
    # kmer_patt = re.compile('[^N]{{{0},}}({1})[^N]{{{0},}}'.format(k_size, site))

    OUTER_IX = np.array([True] * (flank_size+2) +
                        [False] * (site_size*2 + 1 - 4) +
                        [True] * (flank_size+2), dtype=np.bool)
    INNER_IX = ~ OUTER_IX

    print('Counting FastQ reads...')
    n_reads = count_fastq_sequences(args.FASTQ)
    print('Found {} reads'.format(n_reads))

    # probability of acceptance for subsampling
    p_accept = 1.0 if max_reads >= n_reads*10 else max_reads / n_reads * 10
    print('Acceptance threshold: {:.2f}'.format(p_accept))

    # set up random number generation
    if args.seed is None:
        args.seed = np.random.randint(1000000, 99999999)
        print('Random seed was not set, using {}'.format(args.seed))
    else:
        print('Random seed was {}'.format(args.seed))
    random_state = np.random.RandomState(args.seed)
    unif = random_state.uniform
    randint = random_state.randint

    with tqdm.tqdm(total=max_reads) as progress:

        cov_obs = []

        # set up the generator over FastQ reads
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
                    ix = randint(k_size, seq_len - (k_size + site_size)+1)
                    # if no N within the subsequence, then accept it
                    if 'N' not in seq[ix - k_size: ix + k_size + site_size]:
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

            mean_inner, mean_outer = collect_coverage(seq, ix, site_size, k_size)

            # avoid regions with pathologically high coverage
            if mean_outer > args.max_coverage:
                continue

            # record this observation
            cov_obs.append(CovInfo(mean_inner, mean_outer, rtype)) #, _id, ix))
            progress.update()

            # if we have reached a requested number of observations
            # then stop processing
            if args.max_n is not None and len(cov_obs) == args.max_n:
                break

    # lets do some tabular munging
    df = pandas.DataFrame(cov_obs)

    # remove any row which had zero coverage in inner or outer region
    z_inner = df.mean_inner == 0
    z_outer = df.mean_outer == 0
    nz_either = ~(z_inner | z_outer)
    df = df[nz_either]
    print('Rows removed with coverage: inner {}, outer {}, shared {}'.format(
        sum(z_inner), sum(z_outer), sum(z_inner & z_outer)))

    n_sampled = len(df)

    df['ratio'] = df.mean_inner / df.mean_outer
    df.sort_values('ratio', inplace=True)

    agg_type = df.groupby('read_type').size()
    print('Collected observation breakdown. WGS: {} Hi-C: {}'.format(agg_type.wgs, agg_type.hic))

    # the suspected non-hic observations
    wgs = df.loc[df.read_type == 'wgs'].copy()
    wgs.reset_index(inplace=True, drop=True)
    wgs['pvalue'] = 1.0 / len(wgs) * (wgs.index + 1)

    # the suspected hi-c observations
    hic = df.loc[df.read_type == 'hic'].copy()
    hic['pvalue'] = None

    # combine them together
    if args.output is not None:
        print('Writing observations to gzipped csv file: {}'.format(args.output))
        df = wgs.append(hic)
        df.sort_values('ratio', inplace=True)
        df.reset_index(inplace=True, drop=True)
        df.to_csv(gzip.open(args.output, 'wt'), sep='\t')

    # compute the number of HiC reads exceeding a p-value threshold
    wgs_pval_ratio = wgs.iloc[(wgs['pvalue'] < args.p_value).sum()].ratio
    n_hic = (hic.ratio < wgs_pval_ratio).sum()

    # report as a fraction of all n_reads
    print("For p-value {:.3g}, {} reads from {} are Hi-C. Estimated fraction: {:4g}".format(
        args.p_value, n_hic, n_sampled, n_hic / n_sampled))
