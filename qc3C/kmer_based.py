import bz2
import gzip
import numpy as np
import pandas
import subprocess
import tqdm
import logging

from collections import namedtuple
from typing import TextIO, Optional, Dict
from qc3C.exceptions import InsufficientDataException, UnknownLibraryKitException, MaxObsLimit
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance
from qc3C.utils import init_random_state, write_jsonline, count_sequences, observed_fraction
from qc3C._version import runtime_info

try:
    import dna_jellyfish
except ImportError:
    import jellyfish as dna_jellyfish

logger = logging.getLogger(__name__)


CovInfo = namedtuple('cov_info', ['mean_inner', 'mean_outer', 'read_type', 'read_pos', 'seq'])


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

    def analysed(self) -> int:
        return self.counts['all'] - self.counts['sample']

    def accepted(self) -> int:
        return self.analysed() - self.rejected()

    def rejected(self) -> int:
        return self.counts['short'] + self.counts['flank'] + self.counts['ambig'] + self.counts['high_cov']

    def count(self, category: str) -> int:
        """
        Return the number of counts in a particular category (analysed, rejected or ...)
        :param category: the category
        :return: the counts of category
        """
        if category == 'analysed':
            return self.analysed()
        elif category == 'rejected':
            return self.rejected()
        elif category == 'accepted':
            return self.accepted()
        else:
            return self.counts[category]

    def fraction(self, category: str, against: str = 'analysed') -> float:
        """
        Return the category's fraction compared to one of (accepted, all, analysed).
        :param category: the category (numerator)
        :param against: the denominator
        :return: a fraction [0,1]
        """
        if against == 'accepted':
            return self.count(category) / self.accepted()
        elif against == 'all':
            return self.count(category) / self.counts['all']
        elif against == 'analysed':
            return self.count(category) / self.analysed()
        else:
            raise RuntimeError('parameter \"against\" must be one of [accepted, analysed or all]')


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


def assign_empirical_pvalues_all(df):
    """
    For a combined wgs/hic dataframe, assign a ranked p-value on ascending ratio using
    only the wgs observations. Finally, distribute p-values to hic observations based
    on their order in the ascending table.

    :param df: the dataframe to assign pvalues
    :return: a new 3-column dataframe ready for pvalue distribution
    """

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
        # begin with lowest assigned pvalue
        current_p = df.loc[df.read_type == 'wgs', 'pvalue'].values[0]
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

    # ascending ratio order
    df = df.sort_values('ratio')
    # extract just the relevant columns as a numpy array
    wgs_obs = df.loc[df.read_type == 'wgs'].to_numpy()

    # assign p-values overall, only incrementing when the ratio changes
    n = 1
    ix_pval = df.columns.get_loc('pvalue')
    ix_ratio = df.columns.get_loc('ratio')
    wgs_obs[0, ix_pval] = n
    last_ratio = wgs_obs[0, ix_ratio]
    for i in range(1, wgs_obs.shape[0]):
        if wgs_obs[i, ix_ratio] != last_ratio:
            n += 1
        last_ratio = wgs_obs[i, ix_ratio]
        wgs_obs[i, ix_pval] = n
    wgs_obs[:, ix_pval] /= n

    # reconstruct a dataframe for convenience by stacking together the wgs and hic obs
    df = pandas.DataFrame(np.vstack([wgs_obs, df.loc[df.read_type == 'hic'].to_numpy()]),
                          columns=df.columns)
    # reestablish order by ascending ratio
    df.sort_values('ratio', inplace=True)
    df.reset_index(drop=True, inplace=True)
    # lastly distribute p-values to hi-c obs
    df = distribute_pvalues(df)
    return df


def pvalue_expectation(df: pandas.DataFrame) -> Dict[str, float]:
    """
    Calculate the mean p-value and variation for Hi-C observations. WGS observations
    are assigned a p-value of 0.
    :param df: the dataframe to analyse
    :return: a tuple of p-value mean and error
    """
    # number of total observations
    n_obs = len(df)
    # Hi-C q = 1 - p-value
    q = 1 - df.loc[df.read_type == 'hic', 'pvalue'].to_numpy(np.float)

    # mean value and error
    _mean = q.sum() / n_obs
    _err = np.sqrt((q * (1-q)).sum()) / n_obs

    return {'mean': _mean, 'error': _err}


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


def get_kmer_frequency_cutoff(filename: str, q: float = 0.9, threads: int = 1) -> int:
    """
    Get a cutoff value for high frequency kmers by calling out to Jellyfish.
    The external call produces a histogram of kmer frequency, from which Pandas is used to
    calculate a quantile cut-off.

    :param filename: A jellyfish k-mer database
    :param q: the threshold quantile
    :param threads: number of concurrent threads to use
    :return: maximum frequency cut-off
    """
    proc = subprocess.Popen(['jellyfish', 'histo', '-t{}'.format(threads), filename],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        max_cov = pandas.read_csv(proc.stdout, sep=' ', names=['cov', 'count'])['cov'].quantile(q=q)
        return round(float(max_cov))
    except Exception as e:
        raise RuntimeError('Encountered a problem while analyzing kmer distribution. [{}]'.format(e))


def next_read(filename: str, site: str, k_size: int,
              prob_accept: float, progress, counts: Counter,
              random_state: np.random.RandomState) -> (str, Optional[int], str, int):
    """
    Create a generator function which returns sequence data site positions
    from a FastQ file.
    :param filename: the fastq file, compressed or not
    :param site: the target site to find in reads, eg GATCGATC
    :param k_size: the k-mer size used in the counting database
    :param prob_accept: probability threshold used to subsample the full read-set
    :param progress: progress bar for user feedback
    :param counts: a tracker class for various rejection conditions
    :param random_state: random state for subsampling
    :return: yield tuple of (sequence, site_position, read_id, sequence length)
    """
    reader = read_seq(open_input(filename))
    site_size = len(site)
    unif = random_state.uniform

    for _id, _seq, _qual in reader:

        counts['all'] += 1
        progress.update()

        # subsampling read-set
        if prob_accept is not None and prob_accept < unif():
            counts['sample'] += 1
            continue

        # skip sequences which are too short to analyse
        seq_len = len(_seq)
        if seq_len < 2 * k_size + site_size + 1:
            counts['short'] += 1
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
            counts['flank'] += 1


def analyse(enzyme: str, kmer_db: str, read_files: list, mean_insert: int, seed: int = None,
            sample_rate: float = None, max_freq_quantile: float = 0.9, threads: int = 1, max_obs: int = None,
            output_table: str = None, report_path: str = None, library_kit: str = 'generic') -> None:
    """
    Using a read-set and its associated Jellyfish kmer database, analyse the reads for evidence
    of proximity junctions.
        
    :param enzyme: the enzyme used during digestion
    :param kmer_db: the jellyfish kmer database
    :param read_files: the list of read files in FastQ format.
    :param mean_insert: mean length of inserts used in creating the library
    :param seed: random seed used in subsampling read-set
    :param sample_rate: probability of accepting an observation. If None accept all.
    :param max_freq_quantile: ignore kmers with kmer frequencies above this quantile
    :param threads: use additional threads for supported steps
    :param max_obs: the maximum number of reads to inspect
    :param output_table: if not None, write the full pandas table to the specified path
    :param report_path: append a report in single-line JSON format to the given path.
    :param library_kit: the type of kit used in producing the library (ie. phase, generic)
    """

    def collect_coverage(seq: str, ix: int, site_size: int, k: int, min_cov: int = 0) -> (float, float):
        """
        Collect the k-mer coverage centered around the position ix. From the left, the sliding
        window begins just before the site region and slides right until just after. Means
        are then calculated for the inner (within the junction) and outer (left and right flanks)
        :param seq: the sequence to analyse
        :param ix: the position marking the beginning of a junction or any other arbitrary location if so desired
        :param site_size: the size of the junction site
        :param k: the kmer size
        :param min_cov: apply a minimum value to coverage reported by jellyfish.
        :return: mean(inner), mean(outer)
        """
        assert k <= ix <= len(seq) - (k + site_size), \
            'The site index {} is either too close to start (min {}) or ' \
            'end (max {}) of read to scan for coverage'.format(ix, k, len(seq) - (k + site_size))

        outer = []
        for i in range(0, 5):
            smer = seq[ix-k_size+i: ix+i]
            mer = dna_jellyfish.MerDNA(smer)
            mer.canonicalize()
            k_cov = query_jf[mer]
            if k_cov < min_cov:
                k_cov = min_cov
            outer.append(k_cov)

        for i in range(4, -1, -1):
            smer = seq[ix+junc_size-i:ix+junc_size+k_size-i]
            mer = dna_jellyfish.MerDNA(smer)
            mer.canonicalize()
            k_cov = query_jf[mer]
            if k_cov < min_cov:
                k_cov = min_cov
            outer.append(k_cov)

        inner = []
        n = (k_size - junc_size) // 2
        for i in range(-(n-2), n, 2):
            smer = seq[ix+i-n: ix+i-n+24]
            mer = dna_jellyfish.MerDNA(smer)
            mer.canonicalize()
            k_cov = query_jf[mer]
            if k_cov < min_cov:
                k_cov = min_cov
            inner.append(k_cov)

        inner = np.array(inner)
        outer = np.array(outer)
        return max(1, inner.min()), max(1, outer.min())

    def set_progress_description():
        """ Convenience method for setting tqdm progress bar description """
        if len(read_files) > 1:
            progress.set_description('Progress (file {})'.format(n))
        else:
            progress.set_description('Progress')

    # Currently there is only special logic for Phase kits
    # whose proximity inserts appear to possess a non-uniform
    # distribution of junction location
    if library_kit == 'phase':
        logger.info('Phase Genomics library kit, treating junction sites as non-uniformly distributed')
        is_phase = True
    elif library_kit == 'generic':
        logger.info('Generic Hi-C library kit, treating junction sites as uniformly distributed')
        is_phase = False
    else:
        raise UnknownLibraryKitException(library_kit)

    # Determine the junction, treat as uppercase
    lig_info = ligation_junction_seq(get_enzyme_instance(enzyme))
    junc_site = lig_info.junction
    junc_size = lig_info.junc_len
    cut_site = lig_info.cut_site

    # TODO flexible flanks could be handled if L/R flanks treated independently
    k_size = get_kmer_size(kmer_db)
    flank_size = (k_size - junc_size) // 2

    logger.info('Found a k-mer size of {}nt for the Jellyfish library at {}'.format(k_size, kmer_db))
    logger.info('The enzyme {} with cut-site {} produces an {}nt ligation junction sequence of {}'
                .format(lig_info.enzyme_name, lig_info.elucidation, lig_info.junc_len, lig_info.junction))
    logger.info('For this k-mer size and enzyme, the flanks are {}nt'.format(flank_size))
    logger.info('Minimum usable read length for this k-mer size and enzyme is {}nt'.format(k_size + 2*flank_size))

    max_coverage = get_kmer_frequency_cutoff(kmer_db, max_freq_quantile, threads)
    logger.info('Maximum k-mer frequency of {} at quantile {} in Jellyfish library {}'
                .format(max_coverage, max_freq_quantile, kmer_db))

    # some sanity checks.
    assert k_size > junc_size, 'Kmer size must be larger the the junction size'
    assert (k_size - junc_size) % 2 == 0, 'Kmer size and junction size should match (even/even) or (odd/odd)'

    # for efficiency, disable sampling if rate is 1
    if sample_rate == 1:
        logger.info('Disabling sampling as requested rate was equal to 1')
        sample_rate = None

    # initialize jellyfish API
    dna_jellyfish.MerDNA_k(k_size)
    query_jf = dna_jellyfish.QueryMerFile(kmer_db)

    # probability of acceptance for subsampling
    if sample_rate is None or sample_rate == 1:
        logger.info('Accepting all usable reads')
        sample_rate = None
    else:
        logger.info('Acceptance threshold: {:#.2g}'.format(sample_rate))

    # set up random number generation
    random_state = init_random_state(seed)
    randint = random_state.randint

    starts_with_cutsite = 0
    cumulative_length = 0

    logger.info('Beginning analysis...')

    # count the available reads if not max_obs not set
    if max_obs is not None:
        logger.info('Sampling a maximum of {} observations'.format(max_obs))
        max_obs = np.linspace(0, max_obs, len(read_files) + 1, dtype=np.int)
        if len(read_files) > 1:
            logger.debug('Sampling observations equally from each of the {} read files'.format(len(read_files)))
        # prepare the progress bar
        progress = tqdm.tqdm(total=max_obs[-1])
    else:
        n_reads = 0
        for reads in read_files:
            logger.info('Counting reads in {}'.format(reads))
            n_reads += count_sequences(reads, 'fastq', max_cpu=threads)
        logger.info('Found {:,} reads to analyse'.format(n_reads))
        # prepare the progress bar
        progress = tqdm.tqdm(total=n_reads)

    analysis_counter = Counter()
    coverage_obs = []

    try:

        read_count = analysis_counter.counts

        for n, reads in enumerate(read_files, 1):

            set_progress_description()

            # set up the generator over FastQ reads, with sub-sampling
            fq_reader = next_read(reads, junc_site, k_size, sample_rate, progress,
                                  analysis_counter, random_state)

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
                        analysis_counter.counts['ambig'] += 1
                        continue

                # junction containing reads, categorised as hic for simplicity
                else:

                    rtype = 'hic'

                    # abandon this sequence if it contains an N
                    if 'N' in seq[ix - k_size: ix + k_size + junc_size]:
                        analysis_counter.counts['ambig'] += 1
                        continue

                mean_inner, mean_outer = collect_coverage(seq, ix, junc_size, k_size, min_cov=1)

                # avoid regions with pathologically high coverage
                if mean_inner > max_coverage or mean_outer > max_coverage:
                    analysis_counter.counts['high_cov'] += 1
                    continue

                cumulative_length += seq_len

                # record this observation
                coverage_obs.append(CovInfo(mean_inner, mean_outer, rtype, ix, _id))

                if max_obs is not None and read_count['all'] >= max_obs[n]:
                    break

            if max_obs is not None and read_count['all'] >= max_obs[-1]:
                raise MaxObsLimit

    except MaxObsLimit:
        if max_obs is not None and analysis_counter.count('all') >= max_obs[-1]:
            logger.info('Reached user-defined observation limit [{}]'.format(max_obs[-1]))
    finally:
        if progress is not None:
            progress.close()
            progress = None

    if analysis_counter.fraction('rejected') > 0.9:
        logger.warning('More than 90% of reads were rejected')

    #
    # Reporting
    #
    report = {
        'mode': 'kmer',
        'runtime_info': runtime_info(),
        'input_args': {'kmer_db': kmer_db,
                       'read_list': read_files,
                       'seed': seed,
                       'sample_rate': sample_rate,
                       'max_coverage': max_coverage,
                       'mean_insert': mean_insert,
                       'max_freq_quantile': max_freq_quantile,
                       'library_lit': library_kit},
        'n_parsed_reads': analysis_counter.counts['all'],
        'n_analysed_reads': analysis_counter.analysed(),
        'n_too_short': analysis_counter.count('short'),
        'n_no_flank': analysis_counter.count('flank'),
        'n_ambiguous': analysis_counter.count('ambig'),
        'n_high_cov': analysis_counter.count('high_cov'),
    }

    logger.info('Number of parsed reads: {:,}'
                .format(analysis_counter.counts['all']))
    logger.info('Number of analysed reads: {:,} ({:#.4g}% of parsed)'
                .format(analysis_counter.analysed(),
                        analysis_counter.fraction('analysed', 'all') * 100))
    logger.info('Number of reads filtered [too short]: {:,} ({:#.4g}% of analysed)'
                .format(analysis_counter.count('short'),
                        analysis_counter.fraction('short') * 100))
    logger.info('Number of reads filtered [insufficient flanks]: {:,} ({:#.4g}% of analysed)'
                .format(analysis_counter.count('flank'),
                        analysis_counter.fraction('flank') * 100))
    logger.info('Number of reads filtered [ambiguous sequence]: {:,} ({:#.4g}% of analysed)'
                .format(analysis_counter.count('ambig'),
                        analysis_counter.fraction('ambig') * 100))
    logger.info('Number of reads filtered [high coverage]: {:,} ({:#.4g}% of analysed)'
                .format(analysis_counter.count('high_cov'),
                        analysis_counter.fraction('high_cov') * 100))

    try:
        # handle event that no reads passed through parsing.
        if analysis_counter.accepted() == 0:
            raise InsufficientDataException('No reads were accepted during parsing.')

        report['n_accepted_reads'] = analysis_counter.accepted()
        logger.info('Number of accepted reads: {:,} ({:#.4g}% of analysed)'
                    .format(analysis_counter.accepted(),
                            analysis_counter.fraction('accepted') * 100))

        # lets do some tabular munging, making sure that our categories are explicit
        all_df = pandas.DataFrame(coverage_obs)
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
        report['n_without_junc'] = count_rtype.wgs
        report['n_with_junc'] = count_rtype.hic
        logger.info('Break down of accepted reads: wgs {:,} junction {:,}'.format(count_rtype.wgs, count_rtype.hic))
        if count_rtype.wgs == 0:
            raise InsufficientDataException('No reads classified as wgs were observed.')
        if count_rtype.hic == 0:
            raise InsufficientDataException('No reads containing junctions were observed.')

        mean_read_len = cumulative_length / analysis_counter.count('accepted')
        report['mean_readlen'] = mean_read_len
        logger.info('Observed mean read length: {:.0f}nt'.format(mean_read_len))

        # k-mer coverage ratio. inner (junction region) vs outer (L+R flanking regions)
        all_df['ratio'] = all_df.mean_inner / all_df.mean_outer
        all_df['pvalue'] = None

        report['cs_start'] = starts_with_cutsite
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
        all_df = assign_empirical_pvalues_all(all_df)

        # estimation hic fraction from empirical p-values
        hic_frac = pvalue_expectation(all_df)
        report['raw_fraction'] = hic_frac
        logger.info('Estimated Hi-C read fraction via p-value sum method: {:#.4g} \u00b1 {:#.4g} %'
                    .format(hic_frac['mean'] * 100, hic_frac['error'] * 100))

        if mean_insert is not None:
            unobserved_fraction = 1 - observed_fraction(is_phase, int(mean_read_len), mean_insert, k_size, junc_size)

            if unobserved_fraction < 0:
                unobserved_fraction = 0
                logger.warning('For supplied insert length of {:.0f}nt, estimation of the unobserved fraction '
                               'is invalid (<0). Assuming an unobserved fraction: {:#.4g}'
                               .format(mean_insert, unobserved_fraction))
            else:
                report['mean_insert'] = mean_insert
                logger.info('For supplied insert length of {:.0f}nt, estimated unobserved fraction: {:#.4g}'
                            .format(mean_insert, unobserved_fraction))

                report['adj_fraction'] = {'mean': hic_frac['mean'] + (hic_frac['mean'] * unobserved_fraction),
                                          'error': hic_frac['error'] + (hic_frac['error'] * unobserved_fraction)}

                logger.info('Adjusted estimation of Hi-C read fraction: {:#.4g} \u00b1 {:#.4g} %'
                            .format((hic_frac['mean'] + (hic_frac['mean'] * unobserved_fraction)) * 100,
                                    (hic_frac['error'] + (hic_frac['error'] * unobserved_fraction)) * 100))

            report['unobs_frac'] = unobserved_fraction

        # combine them together
        if output_table is not None:
            import os
            logger.info('Writing observations to gzipped tab-delimited file: {}'.format(output_table))
            with gzip.open(output_table, 'wt') as out_h:
                all_df.to_csv(out_h, sep='\t')

    finally:
        # append report to file
        if report_path is not None:
            write_jsonline(report_path, report)
