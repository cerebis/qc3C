import numpy as np
import tqdm
import pandas
import subprocess
import bz2
import gzip
import logging

from collections import namedtuple
from typing import TextIO, Optional, Dict
from qc3C.exceptions import InsufficientDataException
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance
from qc3C.utils import init_random_state, test_for_exe
import jelly
import io

try:
    import dna_jellyfish
except ImportError:
    import jellyfish as dna_jellyfish

logger = logging.getLogger(__name__)


CovInfo = namedtuple('cov_info', ['mean_inner', 'mean_outer', 'read_type', 'read_pos'])


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


def assign_empirical_pvalues_uniq(df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Taking in a mixed DataFrame containing both wgs and hic observations, assign empirical
    p-values to the wgs observations. Duplicate ratio values are dropped for both WGS and Hi-C
    observations (separately). Consequently, the returned DataFrame contains only these unique rows.
    :param df: a combined DataFrame of both Hi-C and WGS observations
    :return: a DataFrame with only read_type, ratio and pvalue columns
    """
    # columns to use and return
    _cols = ('read_type', 'ratio', 'pvalue')
    # extract the WGS obs, remove duplicates and sort
    wgs_uniq = df.loc[df.read_type == 'wgs', _cols].drop_duplicates('ratio').sort_values('ratio')
    # get the underlying numpy array
    arr = wgs_uniq.to_numpy()
    # impose linear spacing of p-value bins
    pval_bin = 1 / arr.shape[0]
    # distribute linearly increasing p-values across the sorted WGS observations
    arr[:, 2] = np.arange(1, arr.shape[0]+1) * pval_bin
    # extract the Hi-C obs and remove dupes
    hic_uniq = df.loc[df.read_type == 'hic', _cols].drop_duplicates('ratio')
    hic_uniq.pvalue = None
    # recombine the WGS and Hi-C observations as a DataFrame, sort and reindex.
    return pandas.DataFrame(arr, columns=_cols).append(hic_uniq).sort_values('ratio').reset_index(drop=True)


def distribute_pvalues(df: pandas.DataFrame) -> pandas.DataFrame:
    """
    Traverse the table of observations by ascending ratio. When encountering a Hi-C
    observation (which should be unassigned) assign the last p-value seen from the
    WGS observations.
    :param df: the DataFrame of WGS and Hi-C observations, sorted by ascending ratio
    :return: a DataFrame with p-values assigned to Hi-C observations
    """
    # get underlying numpy array and column indices
    arr = df.to_numpy()
    ix_pval = df.columns.get_loc('pvalue')
    ix_rt = df.columns.get_loc('read_type')
    # the most recent p-value
    current_p = arr[0, ix_pval]
    for i in range(arr.shape[0]):
        # update p-value
        if arr[i, ix_rt] == 'wgs':
            current_p = arr[i, ix_pval]
        # distribute current p-value to hi-c obs
        elif arr[i, ix_rt] == 'hic':
            arr[i, ix_pval] = current_p
        else:
            raise RuntimeError('unknown read type {}'.format(arr[i, ix_rt]))
    return pandas.DataFrame(arr, columns=df.columns)


def pvalue_expectation(df: pandas.DataFrame) -> Dict[str, float]:
    """
    Calculate the mean p-value and variation for Hi-C observations. WGS observations
    are assigned a p-value of 0.
    :param df: the dataframe to analyze
    :return: a tuple of p-value mean and error
    """
    # number of total observations
    N = len(df)
    # Hi-C q = 1 - p-value
    q = 1 - df.loc[df.read_type == 'hic', 'pvalue'].to_numpy(np.float)

    # mean value and error
    _mean = q.sum() / N
    _err = np.sqrt((q * (1-q)).sum()) / N

    return {'mean': _mean, 'error': _err}


def analyze(enzyme: str, kmer_db: str, read_list: list, mean_insert: int, seed: int = None,
            sample_rate: float = None, max_coverage: int = 500, threads: int = 1,
            output_table: str = None, num_obs: int = None) -> None:
    """
    Using a read-set and its associated Jellyfish kmer database, analyze the reads for evidence
    of proximity junctions.
        
    :param enzyme: the enzyme used during digestion
    :param kmer_db: the jellyfish kmer database
    :param read_list: the list of read files in FastQ format.
    :param mean_insert: mean length of inserts used in creating the library
    :param seed: random seed used in subsampling read-set
    :param sample_rate: probability of accepting an observation. If None accept all.
    :param max_coverage: ignore kmers with coverage greater than this value
    :param threads: use additional threads for supported steps
    :param output_table: if not None, write the full pandas table to the specified path
    :param num_obs: the number of observations to collect before rendering report
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

    def collect_coverage_pybind(seq: str, ix: int, site_size: int, k: int, min_cov: int = 0) -> (float, float):
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
            mer = jelly.MerDNA(smer)
            mer.canonicalize()
            k_cov = query_jelly[mer]
            if k_cov < min_cov:
                k_cov = min_cov
            sliding_cov[i] = k_cov
        return np.mean(sliding_cov[INNER_IX]), np.mean(sliding_cov[OUTER_IX])

    class Counter(object):
        def __init__(self, **extended_fields):
            self.counts = {'all': 0, 'sample': 0, 'short': 0, 'flank': 0, 'ambig': 0, 'high_cov': 0}
            self.counts.update(extended_fields)

        def __getitem__(self, item: str):
            return self.counts[item]

        def __setitem__(self, key: str, value: int):
            assert key in self.counts, 'Key value {} was not found'.format(key)
            self.counts[key] = value

        def __iadd__(self, other):
            if isinstance(other, Counter):
                _counts = other.counts
            elif isinstance(other, dict):
                _counts = other
            else:
                raise RuntimeError()
            for k, v in _counts.items():
                self.counts[k] += v
            return self

        def update(self, _dict):
            for k, v in _dict.items():
                self.counts[k] = v

        def analyzed(self) -> int:
            return self.counts['all'] - self.counts['sample']

        def accepted(self) -> int:
            return self.analyzed() - self.rejected()

        def rejected(self) -> int:
            return self.counts['short'] + self.counts['flank'] + self.counts['ambig'] + self.counts['high_cov']

        def count(self, category: str) -> int:
            """
            Return the number of counts in a particular category (analyzed, rejected or ...)
            :param category: the category
            :return: the counts of category
            """
            if category == 'analyzed':
                return self.analyzed()
            elif category == 'rejected':
                return self.rejected()
            elif category == 'accepted':
                return self.accepted()
            else:
                return self.counts[category]

        def fraction(self, category: str, against: str = 'analyzed') -> float:
            """
            Return the category's fraction compared to one of (accepted, all, analyzed).
            :param category: the category (numerator)
            :param against: the denominator
            :return: a fraction [0,1]
            """
            if against == 'accepted':
                return self.count(category) / self.accepted()
            elif against == 'all':
                return self.count(category) / self.counts['all']
            elif against == 'analyzed':
                return self.count(category) / self.analyzed()
            else:
                raise RuntimeError('parameter \"against\" must be one of [accepted, analyzed or all]')

    def next_read(filename: str, site: str, k_size: int,
                  prob_accept: float, progress, counts: Counter) -> (str, Optional[int], str, int):
        """
        Create a generator function which returns sequence data site positions
        from a FastQ file.
        :param filename: the fastq file, compressed or not
        :param site: the target site to find in reads, eg GATCGATC
        :param k_size: the k-mer size used in the counting database
        :param prob_accept: probability threshold used to subsample the full read-set
        :param progress: progress bar for user feedback
        :param counts: a tracker class for various rejection conditions
        :return: yield tuple of (sequence, site_position, read_id, sequence length)
        """
        reader = read_seq(open_input(filename))
        site_size = len(site)

        _all = 0
        _flank = 0
        _sample = 0
        _short = 0
        for _id, _seq, _qual in reader:

            _all += 1
            progress.update()

            # subsampling read-set
            if prob_accept is not None and prob_accept < unif():
                _sample += 1
                continue

            # skip sequences which are too short to analyze
            seq_len = len(_seq)
            if seq_len < 2 * k_size + site_size + 1:
                _short += 1
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
            else:
                _flank += 1

        counts.update({'all': _all, 'sample': _sample, 'short': _short, 'flank': _flank})

    def get_kmer_size(filename: str) -> int:
        """
        Retrieve the k-mer size used in a library by reading the first record.
        This assumes all mers within the library have the same k-mer size.
        :param filename:
        :return: the size of the k-mer used in a Jellyfish library
        """
        iter_mer = dna_jellyfish.ReadMerFile(filename)
        current_mer = iter_mer.mer()
        k_size = current_mer.k()
        del iter_mer
        return k_size

    # Determine the junction, treat as uppercase
    lig_info = ligation_junction_seq(get_enzyme_instance(enzyme))
    junc_site = lig_info.junction
    junc_size = lig_info.junc_len
    cut_site = lig_info.cut_site

    # TODO for flexible flanks could be handled if L/R flanks treated independently
    k_size = get_kmer_size(kmer_db)
    flank_size = (k_size - junc_size) // 2

    logger.info('Found a k-mer size of {}nt for the Jellyfish library at {}'.format(k_size, kmer_db))
    logger.info('The enzyme {} with cut-site {} produces an {}nt ligation junction sequence of {}'
                .format(lig_info.enzyme_name, lig_info.elucidation, lig_info.junc_len, lig_info.junction))
    logger.info('For this k-mer size and enzyme, the flanks are {}nt'.format(flank_size))
    logger.info('Minimum usable read length for this k-mer size and enzyme is {}nt'.format(k_size + 2*flank_size))

    # some sanity checks.
    assert k_size > junc_size, 'Kmer size must be larger the the junction size'
    assert (k_size - junc_size) % 2 == 0, 'Kmer size and junction size should match (even/even) or (odd/odd)'

    # for efficiency, disable sampling if rate is 1
    if sample_rate == 1:
        logger.info('Disabling sampling as requested rate was equal to 1')
        sample_rate = None

    # initialize jellyfish API
    # dna_jellyfish.MerDNA_k(k_size)
    # query_jf = dna_jellyfish.QueryMerFile(kmer_db)
    query_jelly = jelly.QueryMerFile(kmer_db)

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

    starts_with_cutsite = 0
    cumulative_length = 0

    logger.info('Beginning analysis...')

    if len(read_list) > 1:
        initial_desc = 'Progress (file 1)'
    else:
        initial_desc = 'Progress'

    with tqdm.tqdm(desc=initial_desc, total=n_reads) as progress:

        analysis_counter = Counter()

        cov_obs = []
        for n, reads in enumerate(read_list, 1):

            if n > 1:
                progress.set_description('Progress (file {})'.format(n))

            read_counter = Counter()

            # set up the generator over FastQ reads, with sub-sampling
            fq_reader = next_read(reads, junc_site, k_size, sample_rate, progress, read_counter)

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
                        read_counter.counts['ambig'] += 1
                        continue

                # junction containing reads, categorised as hic for simplicity
                else:

                    rtype = 'hic'

                    # abandon this sequence if it contains an N
                    if 'N' in seq[ix - k_size: ix + k_size + junc_size]:
                        read_counter.counts['ambig'] += 1
                        continue

                mean_inner, mean_outer = collect_coverage_pybind(seq, ix, junc_size, k_size, min_cov=1)

                # avoid regions with pathologically high coverage
                if mean_outer > max_coverage:
                    read_counter.counts['high_cov'] += 1
                    continue

                cumulative_length += seq_len

                # record this observation
                cov_obs.append(CovInfo(mean_inner, mean_outer, rtype, ix))

            if read_counter.fraction('rejected') > 0.9:
                logger.warning('More than 90% of reads were rejected in the file {}'.format(reads))

            analysis_counter += read_counter

    #
    # Reporting
    #

    logger.info('Number of parsed reads: {:,}'
                .format(analysis_counter.counts['all']))
    logger.info('Number of analyzed reads: {:,} ({:#.4g}% of parsed)'
                .format(analysis_counter.analyzed(),
                        analysis_counter.fraction('analyzed', 'all') * 100))
    logger.info('Number of reads filtered [too short]: {:,} ({:#.4g}% of analyzed)'
                .format(analysis_counter.count('short'),
                        analysis_counter.fraction('short') * 100))
    logger.info('Number of reads filtered [insufficient flanks]: {:,} ({:#.4g}% of analyzed)'
                .format(analysis_counter.count('flank'),
                        analysis_counter.fraction('flank') * 100))
    logger.info('Number of reads filtered [ambiguous sequence]: {:,} ({:#.4g}% of analyzed)'
                .format(analysis_counter.count('ambig'),
                        analysis_counter.fraction('ambig') * 100))
    logger.info('Number of reads filtered [high coverage]: {:,} ({:#.4g}% of analyzed)'
                .format(analysis_counter.count('high_cov'),
                        analysis_counter.fraction('high_cov') * 100))

    # handle event that no reads passed through parsing.
    if analysis_counter.accepted() == 0:
        raise InsufficientDataException('No reads were accepted during parsing.')

    logger.info('Number of accepted reads: {:,} ({:#.4g}% of analyzed)'
                .format(analysis_counter.accepted(),
                        analysis_counter.fraction('accepted') * 100))

    # lets do some tabular munging, making sure that our categories are explicit
    all_df = pandas.DataFrame(cov_obs)
    # make read_type categorical so we will always see values for both when counting
    all_df.read_type = pandas.Categorical(all_df.read_type, ['hic', 'wgs'], ordered=False)

    # remove any row which had zero coverage in inner or outer region
    z_inner = all_df.mean_inner == 0
    z_outer = all_df.mean_outer == 0
    nz_either = ~(z_inner | z_outer)
    all_df = all_df[nz_either]
    logger.info('Rows removed with no coverage: inner {:,}, outer {:,}, shared {:,}'.format(
        sum(z_inner), sum(z_outer), sum(z_inner & z_outer)))

    # check that we have some of by read types
    count_rtype = all_df.groupby('read_type').size()
    logger.info('Break down of accepted reads: wgs {:,} junction {:,}'.format(count_rtype.wgs, count_rtype.hic))
    if count_rtype.wgs == 0:
        raise InsufficientDataException('No reads classified as wgs were observed.')
    if count_rtype.hic == 0:
        raise InsufficientDataException('No reads containing junctions were observed.')

    # k-mer coverage ratio. inner (junction region) vs outer (L+R flanking regions)
    all_df['ratio'] = all_df.mean_inner / all_df.mean_outer
    all_df['pvalue'] = 0

    logger.info('Number of accepted reads starting with a cut site: {:,} ({:#.4g}% of accepted)'
                .format(starts_with_cutsite,
                        starts_with_cutsite / analysis_counter.count('accepted') * 100))
    logger.info('Expected fraction by random chance 50% GC: {:#.4g}%'
                .format(1 / 4 ** lig_info.site_len * 100))
    logger.info('Number of accepted reads containing the junction sequence: {:,} ({:#.4g}% of accepted)'
                .format(count_rtype.hic,
                        count_rtype.hic / analysis_counter.count('accepted') * 100))
    logger.info('Expected fraction by random chance 50% GC: {:#.4g}%'
                .format(1 / 4 ** lig_info.junc_len * 100))

    # calculate the table of empirical p-values
    uniq_obs = assign_empirical_pvalues_uniq(all_df)

    logger.info('Number of unique observations used in empirical p-value: wgs {:,}, hic {:,}'
                .format(sum(uniq_obs.read_type == 'wgs'), sum(uniq_obs.read_type == 'hic')))

    # estimation hic fraction from empirical p-values
    hic_frac = pvalue_expectation(distribute_pvalues(uniq_obs))

    logger.info('Estimated Hi-C read fraction via p-value sum method: {:#.4g} \u00b1 {:#.4g} %'
                .format(hic_frac['mean'] * 100, hic_frac['error'] * 100))

    mean_read_len = cumulative_length / analysis_counter.count('accepted')
    logger.info('Observed mean read length for paired reads: {:.0f}nt'.format(mean_read_len))

    if mean_insert is not None:
        unobserved_fraction = (mean_insert - mean_read_len * 2) / mean_insert

        if unobserved_fraction < 0:
            unobserved_fraction = 0
            logger.warning('For supplied insert length of {:.0f}nt, estimation of the unobserved fraction '
                           'is invalid (<0). Assuming an unobserved fraction: {:#.4g}'
                           .format(mean_insert, unobserved_fraction))
        else:
            logger.info('For supplied insert length of {:.0f}nt, estimated unobserved fraction: {:#.4g}'
                        .format(mean_insert, unobserved_fraction))
            logger.info('Adjusting for unobserved junction sequences using average fragment size: {}nt'
                        .format(mean_insert))
            logger.info('Adjusted estimation of Hi-C read fraction: {:#.4g} \u00b1 {:#.4g} %'
                        .format(hic_frac['mean'] * (1 + unobserved_fraction) * 100,
                                hic_frac['error'] * (1 + unobserved_fraction) * 100))

    # combine them together
    if output_table is not None:
        import os
        # append the gzip suffix is required
        if not output_table.endswith('.gz'):
            output_table = '{}.gz'.format(output_table)
        # don't overwrite existing files
        if os.path.exists(output_table):
            logging.WARNING('The path {} already exists, output table was not written'.format(output_table))
        else:
            logger.info('Writing observations to gzipped tab-delimited file: {}'.format(output_table))
            with gzip.open(output_table, 'wt') as out_h:
                all_df.to_csv(out_h, sep='\t')
