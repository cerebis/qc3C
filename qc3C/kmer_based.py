import numpy as np
import tqdm
import pandas
import subprocess
import bz2
import gzip
import logging

from collections import namedtuple
from typing import TextIO, Optional
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance, LigationInfo
from qc3C.utils import init_random_state, test_for_exe

try:
    import dna_jellyfish
except ImportError:
    import jellyfish as dna_jellyfish

logger = logging.getLogger(__name__)


CovInfo = namedtuple('cov_info', ['mean_inner', 'mean_outer', 'read_type'])


def open_input(file_name: str) -> TextIO:
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


def read_seq(fp: TextIO) -> (str, str, Optional[str]):
    """
    Method to quickly read FastA or FastQ files using a generator function.
    Originally sourced from https://github.com/lh3/readfq
    :param fp: input file object
    :return: tuple
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


def count_fastq_sequences(file_name: str, max_cpu: int = 1) -> int:
    """
    Estimate the number of fasta sequences in a file by counting headers. Decompression
    is automatically attempted for files ending in .gz. Counting and decompression is by
    way of subprocess calls to grep and gzip. Uncompressed files are also handled. This
    is about 8 times faster than parsing a file with BioPython and 6 times faster than
    reading all lines in Python.
    :param file_name: the fasta file to inspect
    :param max_cpu: the number of cpus if pigz exists
    :return: the estimated number of records
    """
    if file_name.endswith('.gz'):
        pigz_path = test_for_exe('pigz')
        if max_cpu > 1 and pigz_path is not None:
            proc_uncomp = subprocess.Popen([pigz_path, '-p{}'.format(max_cpu), '-cd', file_name],
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            proc_uncomp = subprocess.Popen(['gzip', '-cd', file_name],
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc_read = subprocess.Popen(['grep', r'^@'], stdin=proc_uncomp.stdout, stdout=subprocess.PIPE)
    else:
        proc_read = subprocess.Popen(['grep', r'^@', file_name], stdout=subprocess.PIPE)

    n = 0
    for _ in proc_read.stdout:
        n += 1
    return n


def print_report(hic: pandas.DataFrame, all: pandas.DataFrame, lig_info: LigationInfo,
                 mean_insert: int, cumulative_length: int, reads_evaluated: int, starts_with_cutsite: int,
                 failed_wgs: int, failed_jnc: int, failed_cov: int, max_coverage: int) -> None:
    """
    Print a report of analysis results.

    :param hic: the table of hi-c observations
    :param all:  the table of all (wgs + hi-c) observations
    :param lig_info: details of possible ligations in experiment
    :param mean_insert: the mean insert size used in library generation
    :param cumulative_length: the cumulative length of analyzed sequencing reads in bp
    :param reads_evaluated: the number of reads evaluated in the analysis
    :param starts_with_cutsite: the number of reads which began with a cut-site
    :param failed_wgs: number of wgs reads abandoned from containing ambiguous sequence
    :param failed_jnc: number of junction reads abandoned from containing ambiguous sequence
    :param failed_cov: number of redas rejected due to excessive coverage
    :param max_coverage: maximum acceptable coverage with parsing reads
    """

    logger.info('Number of reads abandoned due to ambiguous sequence. wgs: {:,} frac: {:#.4g}, junction: {:,} frac: {:#.4g}'
                .format(failed_wgs, failed_jnc, failed_wgs / reads_evaluated, failed_jnc / reads_evaluated))
    logger.info('Number of reads abandoned due to coverage > {}: {:,} frac: {:#.4g}'
                .format(max_coverage, failed_cov, failed_cov / reads_evaluated))

    logger.info('Fraction of reads starting with a cut site: {:#.4g}'.format(starts_with_cutsite / reads_evaluated))
    logger.info('Expected fraction by random chance 50% GC: {:#.4g}'.format(1 / 4 ** lig_info.site_len))
    logger.info('Fraction of reads containing the junction sequence: {:#.4g}'.format(len(hic) / reads_evaluated))

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

    logger.info('Estimated Hi-C read fraction via p-value sum method: {:#.4g} +/- {:#.4g}'.format(fraction_hic, hic_stddev))

    if mean_insert is not None:
        unobserved_fraction = (mean_insert - (cumulative_length / reads_evaluated) * 2) / mean_insert
        if unobserved_fraction < 0:
            logger.warning('Estimated unobserved fraction < 0, therefore adjustment will not be applied')
            logger.warning('Check that the supplied mean insert length is not too short')
        else:
            logger.info('Estimated unobserved fraction: {:#.4g}'.format(unobserved_fraction))
            fraction_hic += fraction_hic * unobserved_fraction
            logger.info('Adjusting for unobserved junction sequences using average fragment size: {}nt'.format(mean_insert))
            logger.info('Adjusted estimation of Hi-C read fraction: {:#.4g} +/- {:#.4g}'.format(fraction_hic, hic_stddev))


def analyze(k_size: int, enzyme: str, kmer_db: str, read_list: list, mean_insert: int, seed: int = None,
            sample_rate: float = None, max_coverage: int = 500, threads: int = 1, save_cov: bool = False) -> None:
    """
    Using a read-set and its associated Jellyfish kmer database, analyze the reads for evidence
    of proximity junctions.
        
    :param k_size: kmer size used when building the kmer database 
    :param enzyme: the enzyme used during digestion 
    :param kmer_db: the jellyfish kmer database
    :param read_list: the list of read files in FastQ format.
    :param mean_insert: mean length of inserts used in creating the library
    :param seed: random seed used in subsampling read-set
    :param sample_rate: probability of accepting an observation. If None accept all.
    :param max_coverage: ignore kmers with coverage greater than this value
    :param threads: use additional threads for supported steps
    :param save_cov: if True, write collected observations to file
    """

    def collect_coverage(seq: str, ix: int, site_size: int, k: int, min_cov: int = 0) -> (float, float):
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
            if k_cov < min_cov:
                k_cov = min_cov
            sliding_cov[i] = k_cov
        return np.mean(sliding_cov[INNER_IX]), np.mean(sliding_cov[OUTER_IX])

    def next_read(filename: str, site: str, k_size: int,
                  prob_accept: float, progress) -> (str, Optional[int], str, int):
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
        reader = read_seq(open_input(filename))
        site_size = len(site)

        for _id, _seq, _qual in reader:
            progress.update()
            seq_len = len(_seq)
            # skip sequences which are too short to analyze
            if seq_len < 2 * k_size + site_size + 1:
                continue
            # subsampling read-set
            if prob_accept is not None and prob_accept < unif():
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
    lig_info = ligation_junction_seq(get_enzyme_instance(enzyme))
    junc_site = lig_info.junction
    junc_size = lig_info.junc_len
    cut_site = lig_info.cut_site
    # TODO for flexible flanks could be handled if L/R flanks treated independently
    flank_size = (k_size - junc_size) // 2

    # some sanity checks.
    assert k_size - junc_size > 0, 'Kmer size must be larger the the junction size'
    assert (k_size - junc_size) % 2 == 0, 'Kmer size and junction size should match (even/even) or (odd/odd)'

    # for efficiency, disable sampling if rate is 1
    if sample_rate == 1:
        logger.info('Disabling sampling as requested rate was equal to 1')
        sample_rate = None

    # initialize jellyfish API
    dna_jellyfish.MerDNA_k(k_size)
    query_jf = dna_jellyfish.QueryMerFile(kmer_db)

    OUTER_IX = np.array([True] * (flank_size+2) +
                        [False] * (junc_size*2 + 1 - 4) +
                        [True] * (flank_size+2), dtype=np.bool)
    INNER_IX = ~ OUTER_IX

    # either count the reads or use the information provided by the user.
    n_reads = 0
    for reads in read_list:
        logger.info('Counting reads in {}'.format(reads))
        n_reads += count_fastq_sequences(reads, max_cpu=threads)
    logger.info('Found {:,} reads to analyse'.format(n_reads))

    # probability of acceptance for subsampling
    if sample_rate is None or sample_rate == 1:
        logger.info('Accepting all usable reads')
        sample_rate = None
    else:
        logger.info('Acceptance threshold: {:#.2g}'.format(sample_rate))

    # set up random number generation
    random_state = init_random_state(seed)
    unif = random_state.uniform
    randint = random_state.randint

    reads_evaluated = 0
    starts_with_cutsite = 0
    cumulative_length = 0
    failed_wgs = 0
    failed_jnc = 0
    failed_cov = 0

    logger.info('Beginning analysis...')

    if len(read_list) > 1:
        initial_desc = 'Progress (file 1)'
    else:
        initial_desc = 'Progress'

    with tqdm.tqdm(desc=initial_desc, total=n_reads) as progress:

        cov_obs = []
        for n, reads in enumerate(read_list, 1):

            if n > 1:
                progress.set_description('Progress (file {})'.format(n))

            # set up the generator over FastQ reads, with sub-sampling
            fq_reader = next_read(reads, junc_site, k_size, sample_rate, progress)

            while True:

                try:
                    seq, ix, _id, seq_len = next(fq_reader)
                except StopIteration:
                    break

                if seq.startswith(cut_site):
                    starts_with_cutsite += 1

                # if the read contains no junction, we might still use it
                # as an example of a shotgun read (wgs)
                if ix is None:

                    rtype = 'wgs'

                    # as there are so many non-junction reads, we need to subsample
                    if unif() > 0.1:
                        reads_evaluated += 1
                        cumulative_length += seq_len
                        continue

                    # try to find a randomly selected region which does not contain
                    # ambiguous bases. Skip the read if we fail in a few attempts
                    _attempts = 0
                    while _attempts < 3:
                        ix = randint(k_size, seq_len - (k_size + junc_size)+1)
                        # if no N within the subsequence, then accept it
                        if 'N' not in seq[ix - k_size: ix + k_size + junc_size]:
                            break
                        # otherwise keep trying
                        _attempts += 1

                    # too many tries occurred, abandon this sequence
                    if _attempts >= 3:
                        failed_wgs += 1
                        continue

                # junction containing reads, categorised as hic for simplicity
                else:

                    rtype = 'hic'

                    # abandon this sequence if it contains an N
                    if 'N' in seq[ix - k_size: ix + k_size + junc_size]:
                        failed_jnc += 1
                        continue

                mean_inner, mean_outer = collect_coverage(seq, ix, junc_size, k_size, min_cov=1)

                # avoid regions with pathologically high coverage
                if mean_outer > max_coverage:
                    failed_cov += 1
                    continue

                reads_evaluated += 1
                cumulative_length += seq_len

                # record this observation
                cov_obs.append(CovInfo(mean_inner, mean_outer, rtype))

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
                           'Please check that the kmer database and read-set are correctly matched and '
                           'that --max-n is not set too small')

    all_df['ratio'] = all_df.mean_inner / all_df.mean_outer

    agg_rtype = all_df.groupby('read_type').size()
    logger.info('Collected observation breakdown. wgs: {:,} junction: {:,}'.format(agg_rtype.wgs, agg_rtype.hic))
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

    print_report(hic_df, all_df, lig_info, mean_insert, cumulative_length,
                 reads_evaluated, starts_with_cutsite,
                 failed_wgs, failed_jnc, failed_cov, max_coverage)

    # combine them together
    if save_cov:
        logger.info('Writing observations to gzipped tsv file: {}'.format('cov_dat.tsv.gz'))
        all_df.to_csv(gzip.open('cov_dat.tsv.gz', 'wt'), sep='\t')
