import numpy as np
import tqdm
import pandas
import subprocess
import bz2
import gzip
import logging

from collections import namedtuple
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance

try:
    import dna_jellyfish
except ImportError:
    import jellyfish as dna_jellyfish

logger = logging.getLogger(__name__)


CovInfo = namedtuple('cov_info', ['mean_inner', 'mean_outer', 'read_type'])


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


def readfq(fp):
    """
    Method to quickly read FastA or FastQ files using a generator function.
    Originally sourced from https://github.com/lh3/readfq
    :param fp:
    :return:
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


def count_fastq_sequences(file_name):
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression
    is automatically attempted for files ending in .gz. Counting and decompression is by
    way of subprocess calls to grep and gzip. Uncompressed files are also handled. This
    is about 8 times faster than parsing a file with BioPython and 6 times faster than
    reading all lines in Python.
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


def print_report(hic, all, site_size, mean_insert, read_length, reads_evaluated, starts_with_cutsite):

    logger.info("Fraction of reads starting with a cut site: {:.3g}".format(starts_with_cutsite / reads_evaluated))
    logger.info("Expected fraction at 50% GC: {:.3g}".format(1 / np.power(4, site_size / 2)))
    logger.info("Fraction of reads containing the junction sequence: {:.3g}".format(len(hic) / reads_evaluated))

    cur_pval = 0
    sum_pvals = 0
    var_pvals = 0
    for row in all[['pvalue', 'read_type']].itertuples():
        if pandas.notnull(row.pvalue):
            cur_pval = row.pvalue
        if row.read_type == 'hic':
            sum_pvals += 1 - cur_pval
            var_pvals += cur_pval * (1 - cur_pval)

    fraction_hic = sum_pvals / reads_evaluated
    hic_stddev = np.sqrt(var_pvals) / reads_evaluated
    if mean_insert is not None:
        unobserved_fraction = (mean_insert - (read_length / reads_evaluated) * 2) / mean_insert
        fraction_hic += fraction_hic * unobserved_fraction

    logger.info("Estimated Hi-C read fraction via p-value sum method: "
                "{:.4g} +/- {:.4g}".format(fraction_hic, hic_stddev))

    if mean_insert is None:
        logger.info("To adjust the estimate for unobserved sequence between paired-end reads,\n"
                    "specify the average library fragment size with the --mean-insert= option")
    else:
        logger.info("The estimate above has been adjusted for unobserved junction sequences\n"
                    "based on an average fragment length of {:.4g}nt".format(mean_insert))


def analyze(k_size, enzymes, kmer_db, fastq, mean_insert, output=None, seed=None, pool_size=None,
            max_reads=None, accept_all=False, max_coverage=500):

    def collect_coverage(seq, ix, site_size, k, min_cov=0):
        """
        Collect the k-mer coverage centered around the position ix. From the left, the sliding
        window begins just before the site region and slides right until just after. Means
        are then calculated for the inner (within the junction) and outer (left and right flanks)
        :param seq: the sequence to analyze
        :param ix: the position marking the beginning of a junction or any other arbitrary location if so desired
        :param site_size: the size of the junction site
        :param k: the kmer size
        :param min_cov: apply a minimum value to coverage reported by jellyfish.
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
            k_cov = query_jf[mer]
            if i == 0 and k_cov < 1:
                k_cov = 1
            sliding_cov[i] = k_cov
        return np.mean(sliding_cov[INNER_IX]), np.mean(sliding_cov[OUTER_IX])

    def next_read(filename, site, k_size, prob_accept, progress):
        """
        Create a generator function which returns sequence data site positions
        from a FastQ file.
        :param filename: the fastq file, compressed or not
        :param site: the target site to find in reads, eg GATCGATC
        :param k_size: the k-mer size used in the counting database
        :param prob_accept: probability threshold used to subsample the full read-set
        :param progress: progress bar for user feedback
        :return: yield tuple of (sequence, site_position, read_id, sequence length)
        """
        reader = readfq(open_input(filename))
        site_size = len(site)

        for _id, _seq, _qual in reader:
            progress.update()
            seq_len = len(_seq)
            # skip sequences which are too short to analyze
            if seq_len < 2 * k_size + site_size + 1:
                continue
            # subsampling read-set
            if unif() > prob_accept:
                continue
            # search for the site
            _seq = _seq.upper()
            _ix = _seq.find(site)
            # no site found
            if _ix == -1:
                yield _seq, None, _id, seq_len
            # only report sites which meet flank constraints
            elif k_size <= _ix <= (seq_len - (k_size + site_size)):
                yield _seq, _ix, _id, seq_len

    # Determine the junction, treat as uppercase
    assert len(enzymes) == 1, 'Kmer-based approach currently supports only a single enzyme'
    junction = ligation_junction_seq(get_enzyme_instance(enzymes[0]))
    site = junction.junction.upper()
    site_size = junction.junc_len

    # TODO for flexible flanks could be handled if L/R flanks treated independently
    flank_size = (k_size - site_size) // 2

    # some sanity checks.
    assert k_size - site_size > 0, 'Kmer size must be larger the the junction size'
    assert (k_size - site_size) % 2 == 0, 'Kmer size and junction size should match (even/even) or (odd/odd)'

    # initialize jellyfish API
    dna_jellyfish.MerDNA_k(k_size)
    query_jf = dna_jellyfish.QueryMerFile(kmer_db)

    OUTER_IX = np.array([True] * (flank_size+2) +
                        [False] * (site_size*2 + 1 - 4) +
                        [True] * (flank_size+2), dtype=np.bool)
    INNER_IX = ~ OUTER_IX

    # either count the reads or use the information provided by the user.
    if pool_size is None:
        logger.info('No pool size provided, counting FastQ reads...')
        n_reads = count_fastq_sequences(fastq)
        logger.info('Found {} reads in {}'.format(n_reads, fastq))
    else:
        n_reads = pool_size

    # probability of acceptance for subsampling
    if accept_all:
        prob_accept = 1.0
        logger.info('Accepting all useable reads')
    else:
        prob_accept = 1.0 if max_reads*10 >= n_reads else max_reads / n_reads * 10
        logger.info('Acceptance threshold: {:.2g}'.format(prob_accept))

    # set up random number generation
    if seed is None:
        seed = np.random.randint(1000000, 99999999)
        logger.info('Random seed was not set, using {}'.format(seed))
    else:
        logger.info('Random seed was {}'.format(seed))
    random_state = np.random.RandomState(seed)
    unif = random_state.uniform
    randint = random_state.randint

    reads_evaluated = 0
    starts_with_cutsite = 0
    read_length = 0
    half_site = site[0:site_size//2]

    with tqdm.tqdm(desc="Pool  ", total=n_reads, ncols=80) as pb_read_pool, \
            tqdm.tqdm(desc="Sample", total=max_reads, ncols=80) as pb_sample:

        cov_obs = []

        # set up the generator over FastQ reads, with sub-sampling
        fq_reader = next_read(fastq, site, k_size, prob_accept, pb_read_pool)

        while True:

            # TODO we could probably take this try/catch outside the look
            # TODO as we've more than one place where breaks can occur.
            try:
                seq, ix, _id, seq_len = next(fq_reader)
            except StopIteration:
                # we need to manually close these progress bars if we don't want
                # nearby print statements to interfere
                # TODO doing this makes the with() statement sorta useless
                pb_read_pool.close()
                pb_sample.close()
                logger.info('End of FastQ file reached before filling requested quota')
                break

            reads_evaluated += 1
            if seq[0:site_size//2] == half_site:
                starts_with_cutsite += 1
            read_length += len(seq)

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
                    read_length -= len(seq)
                    reads_evaluated -= 1
                    continue

                rtype = 'hic'

            mean_inner, mean_outer = collect_coverage(seq, ix, site_size, k_size, min_cov=1)

            # avoid regions with pathologically high coverage
            if mean_outer > max_coverage:
                continue

            # record this observation
            cov_obs.append(CovInfo(mean_inner, mean_outer, rtype))
            pb_sample.update()

            # if we have reached a requested number of observations
            # then stop processing
            if max_reads is not None and len(cov_obs) == max_reads:
                break

    # lets do some tabular munging, making sure that our categories are explicit
    all_df = pandas.DataFrame(cov_obs)
    # make read_type categorical so we will always see values for both when counting
    all_df.read_type = pandas.Categorical(all_df.read_type, ['hic', 'wgs'], ordered=False)

    # remove any row which had zero coverage in inner or outer region
    z_inner = all_df.mean_inner == 0
    z_outer = all_df.mean_outer == 0
    nz_either = ~(z_inner | z_outer)
    all_df = all_df[nz_either]
    logger.info('Rows removed with no coverage: inner {}, outer {}, shared {}'.format(
        sum(z_inner), sum(z_outer), sum(z_inner & z_outer)))

    n_sampled = len(all_df)
    if n_sampled == 0:
        raise RuntimeError('The sample set was empty. '
                           'Please check that the kmer database and fastq file are a correct match and '
                           'that --max-n is not set too small')

    all_df['ratio'] = all_df.mean_inner / all_df.mean_outer

    agg_rtype = all_df.groupby('read_type').size()
    logger.info('Collected observation breakdown. WGS: {} junction: {}'.format(agg_rtype.wgs, agg_rtype.hic))
    if agg_rtype.wgs == 0:
        raise RuntimeError('No wgs examples were contained in the collected sample. '
                           'Consider increasing --max-n')
    if agg_rtype.hic == 0:
        raise RuntimeError('No junction examples were contained in the collected sample. '
                           'Consider increasing --max-n')

    # the suspected non-hic observations
    wgs_df = all_df.loc[all_df.read_type == 'wgs'].copy()
    # we require the table to be in ascending ratio order to assign p-values
    wgs_df.sort_values('ratio', inplace=True)
    wgs_df.reset_index(inplace=True, drop=True)
    wgs_df['pvalue'] = 1.0 / len(wgs_df) * (wgs_df.index + 1)

    # the suspected hi-c observations
    hic_df = all_df.loc[all_df.read_type == 'hic'].copy()
    hic_df['pvalue'] = None

    all_df = wgs_df.append(hic_df)
    all_df.sort_values('ratio', inplace=True)
    all_df.reset_index(inplace=True, drop=True)

    print_report(hic_df, all_df, site_size, mean_insert, read_length, reads_evaluated, starts_with_cutsite)

    # combine them together
    if output is not None:
        if not output.endswith('.gz'):
            output = '{}.gz'.format(output)
        logger.info('Writing observations to gzipped csv file: {}'.format(output))
        all_df = wgs_df.append(hic_df)
        all_df.sort_values('ratio', inplace=True)
        all_df.reset_index(inplace=True, drop=True)
        all_df.to_csv(gzip.open(output, 'wt'), sep='\t')
