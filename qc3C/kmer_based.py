import bz2
import gzip
import numpy as np
import pandas
import subprocess
import tqdm
import logging

from collections import namedtuple, OrderedDict, defaultdict
from scipy.stats import gmean
from typing import TextIO, Optional, Dict, List, Tuple
from qc3C.exceptions import *
from qc3C.ligation import Digest
from qc3C.utils import write_jsonline, count_sequences, write_html_report, observed_fraction
from qc3C._version import runtime_info

try:
    import dna_jellyfish
except ImportError:
    import jellyfish as dna_jellyfish

logger = logging.getLogger(__name__)


CovInfo = namedtuple('cov_info', ['mean_inner', 'mean_outer', 'has_junc', 'read_pos', 'seq', 'read_len', 'junc_seq'])


class Counter(object):
    def __init__(self, **extended_fields):
        self.counts = {'all': 0, 'sample': 0, 'short': 0, 'flank': 0, 'ambig': 0,
                       'high_cov': 0, 'low_cov': 0, 'zero_cov': 0}
        self.counts.update(extended_fields)
        self.reject_keys = self.counts.keys() - {'all', 'sample'}

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
        return sum(self.counts[k] for k in self.reject_keys)

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


def assign_empirical_pvalues_all(_df):
    _cols = ['ratio', 'pvalue', 'has_junc', 'read_len', 'seq']

    # get the observations which did not contain a junction
    # these are our cases where the NULL hypothesis is true
    # NOTE as has_junc is categorical, we must use a logical comparison
    _no_junc = _df.loc[_df.has_junc == False, _cols].values
    _no_junc[:, 1] = -1

    # reduce this to the set of unique ratios. The returned
    # ratios are implicitly sorted in ascending order
    _uniq_ratios, _inv_index = np.unique(_no_junc[:, 0], return_inverse=True)

    # The number of unique ratios determines the set of p-values.
    # we use linspace sampler as means of calculating p(r,N) = (r+1)/(N+1)
    _pvals = np.linspace(0, 1, num=_uniq_ratios.shape[0] + 2)[2:]

    # reassign these pvalues using the inverse index
    _no_junc[:, 1] = _pvals[_inv_index]

    # now extract the obs which did contain the junction,
    # we will next distribute the empirical p-values
    _out = _df.loc[_df.has_junc, _cols].values
    _out[:, 1] = -1

    # join both tables and reorder by ascending ratio
    _out = np.vstack([_no_junc, _out])
    _out = _out[np.argsort(_out[:, 0]), ]

    # distribute the smallest non-zero adjacent pvalue to those which are zero
    last_pval = _pvals[0]
    for i in range(0, _out.shape[0]):
        if _out[i, 1] < 0:
            _out[i, 1] = last_pval
        else:
            last_pval = _out[i, 1]

    # remake the full table
    return pandas.DataFrame(_out, columns=_cols)


def pvalue_expectation(df: pandas.DataFrame) -> Dict[str, float]:
    """
    Calculate the mean p-value and variation for Hi-C observations. Junctionless observations
    are assigned a p-value of 0.
    :param df: the dataframe to analyse
    :return: a tuple of p-value mean and error
    """
    # number of total observations
    # n_obs = len(df)
    # Hi-C q = 1 - p-value
    q = 1 - df.loc[df.has_junc, 'pvalue'].to_numpy(np.float)

    # mean value and error
    n_frag_obs = len(set(df.seq))
    _mean = q.sum() / n_frag_obs
    _err = np.sqrt((q * (1-q)).sum()) / n_frag_obs

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


def next_read(filename: str, searcher, longest_site: int, k_size: int,
              prob_accept: float, progress, counts: Counter,
              tracker: Dict[Digest.DerivedString, int], no_user_limit: bool,
              random_state: np.random.RandomState) -> Tuple[str,
                                                            Optional[int],
                                                            str,
                                                            int,
                                                            Optional[int],
                                                            Optional[str]]:
    """
    Create a generator function which returns sequence data site positions
    from a FastQ file.
    :param filename: the fastq file, compressed or not
    :param searcher: regex search method for matching any enzyme
    :param longest_site: the length of the longest cutsite used in digestion
    :param k_size: the k-mer size used in the counting database
    :param prob_accept: probability threshold used to subsample the full read-set
    :param progress: progress bar for user feedback
    :param counts: a tracker class for various rejection conditions
    :param tracker: a tracker for junction kmer types
    :param no_user_limit: a user limit has not been set
    :param random_state: random state for subsampling
    :return: yield tuple of (sequence, site position, read id, seq length, junc length, junc seq)
    """
    reader = read_seq(open_input(filename))
    unif = random_state.uniform

    for _id, _seq, _qual in reader:

        counts['all'] += 1
        if no_user_limit:
            # no limit has been set, track all progress
            progress.update()

        # subsampling read-set
        if prob_accept is not None and prob_accept < unif():
            counts['sample'] += 1
            continue

        # skip sequences which are too short to analyse
        seq_len = len(_seq)
        if seq_len < 2 * k_size + longest_site + 2:
            counts['short'] += 1
            continue

        # search for the site
        _seq = _seq.upper()

        # TODO this (and previous approach) only matches the first instance of the site within the read.
        #  There is no reason that this would be the correct location.
        _match = searcher(_seq)
        # no site found
        if _match is None:
            yield _seq, None, _id, seq_len, None, None
        else:
            _ix, _end = _match.span()
            _junc_len = _end - _ix
            _junc_seq = _match.group()
            tracker[_junc_seq] += 1
            # only report sites which meet flank constraints
            # elif k_size <= _ix <= (seq_len - (k_size + site_size)):
            if k_size + 1 <= _ix <= seq_len - (k_size + _junc_len + 1):
                yield _seq, _ix, _id, seq_len, _junc_len, _junc_seq
            else:
                counts['flank'] += 1


def analyse(enzyme_names: List[str], kmer_db: str, read_files: List[str], mean_insert: int,
            seed: int = None, sample_rate: float = None, max_freq_quantile: float = 0.9, threads: int = 1,
            max_obs: int = None, output_table: str = None, report_path: str = None, no_json: bool = False,
            no_html: bool = False, num_sample: int = 50, frac_sample: float = 1 / 3) -> None:
    """
    Using a read-set and its associated Jellyfish kmer database, analyse the reads for evidence
    of proximity junctions.

    :param enzyme_names: the enzymes used during digestion (max 2)
    :param kmer_db: the jellyfish k-mer database
    :param read_files: list of path names to fastq reads
    :param mean_insert: mean length of inserts used in creating the library
    :param seed: random seed used in subsampling read-set
    :param sample_rate: probability of accepting an observation. If None accept all.
    :param max_freq_quantile: ignore k-mers with k-mer frequencies above this quantile
    :param threads: use additional threads for supported steps
    :param max_obs: the maximum number of reads to inspect
    :param output_table: if not None, write the full pandas table to the specified path
    :param report_path: append a report in single-line JSON format to the given path.
    :param no_json: disable json report
    :param no_html: disable html report
    :param num_sample: the number of bootstrap samples to use
    :param frac_sample: the fraction of observations to use per-bootstrap
    """

    def extract_sample(_df: pandas.DataFrame, _frac: float, _min_n: int, _max_n: int,
                       _rs: np.random.RandomState) -> pandas.DataFrame:
        """
        Extract a sample from a larger table, where random selection is upon fragment id. This
        maintains pairs across the sampling process.
        :param _df: the larger dataframe to select from
        :param _frac: the desired fraction of obs to take in a sample
        :param _min_n: minimum number of obs to take in a sample
        :param _max_n: maximum number of obs to take in a sample
        :param _rs: a numpy random state
        :return: the sampled dataframe
        """
        n_sample = int(len(_df) * _frac)
        if n_sample < _min_n:
            logger.debug('Observation pool is small, limiting sample size to {}'.format(_min_n))
            n_sample = _min_n
        if n_sample > _max_n:
            logger.debug('Observation pool is large, limiting sample size to {}'.format(_max_n))
            n_sample = _max_n

        ix = np.unique(_rs.choice(_df.seq, size=n_sample))
        return _df.set_index('seq').loc[ix].reset_index()

    def collect_coverage(seq: str, ix: int, site_size: int, k: int) -> (float, float):
        """
        Collect the k-mer coverage centered around the position ix. From the left, the sliding
        window begins just before the site region and slides right until just after. Means
        are then calculated for the inner (within the junction) and outer (left and right flanks)
        :param seq: the sequence to analyse
        :param ix: the position marking the beginning of a junction or any other arbitrary location if so desired
        :param site_size: the size of the junction site
        :param k: the kmer size
        :return: mean(inner), mean(outer)
        """
        def get_kmer_cov(mer: str) -> float:
            mer = dna_jellyfish.MerDNA(mer)
            mer.canonicalize()
            k_cov = query_jf[mer]
            return k_cov

        outer = np.empty(shape=6, dtype=np.int)
        inner = np.empty(shape=3, dtype=np.int)
        r_shift = 0
        l_flank = r_flank = (k - site_size) // 2
        if (k - site_size) % 2 != 0:
            l_flank += 1
            r_shift = 1

        assert ix - k > 0, 'smallest valid index is: {}'.format(k+1)
        assert ix + site_size + r_shift + k - 1 < len(seq), \
            'largest valid index is: {}'.format(len(seq) - site_size - k - r_shift)

        for i in range(-1, 2):
            si = seq[ix - k + i: ix + i]
            outer[i+1] = get_kmer_cov(si)
            si = seq[ix + site_size + r_shift + i: ix + site_size + k + r_shift + i]
            outer[i+4] = get_kmer_cov(si)
            si = seq[ix - l_flank + i: ix + site_size + r_flank + i]
            inner[i+1] = get_kmer_cov(si)
            if outer[i+1] == 0 or outer[i+4] == 0 or inner[i+1] == 0:
                # print("OL", si, outer[i+1], seq)
                # print("OR", si, outer[i+4], seq)
                # print("IN", si, inner[i+1], seq)
                raise ZeroCoverageException()

        cov_inner, cov_outer = gmean(inner), gmean(outer)
        if cov_outer < 1:
            raise LowCoverageException()
        return cov_inner, cov_outer

    def set_progress_description():
        """ Convenience method for setting tqdm progress bar description """
        if len(read_files) > 1:
            progress.set_description('Progress (file {})'.format(n))
        else:
            progress.set_description('Progress')

    if seed is None:
        seed = np.random.randint(1000000, 99999999)
        logger.info('Random seed was not set, using {}'.format(seed))

    assert mean_insert > 0, 'Mean insert length must be greater than zero'
    if mean_insert < 100:
        logging.warning('Mean insert length of {}bp is quite short'.format(mean_insert))

    assert 0 < len(enzyme_names) <= 2, 'only 1 or 2 enzymes can be specified'
    digest = Digest(*enzyme_names, no_ambig=False)
    longest_site = digest.longest_junction()

    k_size = get_kmer_size(kmer_db)

    logger.info('Found a k-mer size of {}nt for the Jellyfish library at {}'.format(k_size, kmer_db))

    max_coverage = get_kmer_frequency_cutoff(kmer_db, max_freq_quantile, threads)
    logger.info('Maximum k-mer frequency of {} at quantile {} in Jellyfish library {}'
                .format(max_coverage, max_freq_quantile, kmer_db))

    for li in digest.junctions.values():
        # some sanity checks.
        assert k_size > li.junc_len, 'k-mer size must be larger than any junction size'
        # assert (k_size - li.junc_len) % 2 == 0, 'k-mer size and junction sizes should match (even/even) or (odd/odd)'

        logger.info('The enzymatic combination {}/{} produces the {}nt ligation junction {}'
                    .format(li.enz5p, li.enz3p, li.junc_len, li.junction))
        logger.info('Minimum usable read length for this k-mer size and enzyme is {}nt'
                    .format(li.junc_len + 2 * k_size + 2))

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

    starts_with_cutsite = 0
    cumulative_length = 0

    logger.info('Beginning analysis...')

    # count the available reads if not max_obs not set
    if max_obs is not None:
        user_limited = True
        logger.info('Sampling a maximum of {} observations'.format(max_obs))
        sampled_obs = np.linspace(0, max_obs, len(read_files) + 1, dtype=np.int)
        if len(read_files) > 1:
            logger.debug('Sampling observations equally from each of the {} read files'.format(len(read_files)))
        # prepare the progress bar
        progress = tqdm.tqdm(total=max_obs)
    else:
        user_limited = False
        sampled_obs = None
        n_reads = 0
        for reads in read_files:
            logger.info('Counting reads in {}'.format(reads))
            n_reads += count_sequences(reads, 'fastq', max_cpu=threads)
        logger.info('Found {:,} reads to analyse'.format(n_reads))
        # prepare the progress bar
        progress = tqdm.tqdm(total=n_reads)

    analysis_counter = Counter()
    coverage_obs = []

    any_site_match = digest.cutsite_searcher('startswith')
    junction_tracker = digest.tracker('junction')

    try:

        for n, reads in enumerate(read_files, 1):

            # init same random number state for each file pass
            random_state = np.random.RandomState(seed)
            randint = random_state.randint

            set_progress_description()

            # set up the generator over FastQ reads, with sub-sampling
            fq_reader = next_read(reads, digest.junction_searcher('find'), longest_site, k_size,
                                  sample_rate, progress, analysis_counter,
                                  junction_tracker, not user_limited, random_state)

            while True:

                try:
                    seq, ix, _id, seq_len, junc_len, junc_seq = next(fq_reader)
                except StopIteration:
                    break

                if any_site_match(seq) is not None:
                    starts_with_cutsite += 1

                # if the read contains no junction, we might still use it
                # as an example of a shotgun read (no_junc)
                if ix is None:

                    has_junc = False

                    # assume longest junction for these cases
                    # TODO this should be drawn from the set of possible lengths
                    #  but will only affect Arima currently
                    junc_len = longest_site

                    # try to find a randomly selected region which does not contain
                    # ambiguous bases. Skip the read if we fail in a few attempts
                    _attempts = 0
                    while _attempts < 3:
                        ix = randint(k_size + 1, seq_len - (k_size + longest_site))
                        # if no N within the subsequence, then accept it
                        if 'N' not in seq[ix - k_size: ix + k_size + longest_site + 1]:
                            break
                        # otherwise keep trying
                        _attempts += 1

                    # too many tries occurred, abandon this sequence
                    if _attempts >= 3:
                        analysis_counter.counts['ambig'] += 1
                        continue

                # junction containing reads, categorised as junc for simplicity
                else:

                    has_junc = True

                    # abandon this sequence if it contains an N
                    if 'N' in seq[ix - k_size: ix + k_size + junc_len + 1]:
                        analysis_counter.counts['ambig'] += 1
                        continue

                try:
                    mean_inner, mean_outer = collect_coverage(seq, ix, junc_len, k_size)
                except ZeroCoverageException:
                    analysis_counter.counts['zero_cov'] += 1
                    continue
                except LowCoverageException:
                    analysis_counter.counts['low_cov'] += 1
                    continue

                # avoid regions with pathologically high coverage
                if mean_inner > max_coverage or mean_outer > max_coverage:
                    analysis_counter.counts['high_cov'] += 1
                    continue

                if user_limited:
                    # user limited runs get progress updates here rather than
                    # within the read iterator
                    progress.update()

                cumulative_length += seq_len

                # record this observation
                coverage_obs.append(CovInfo(mean_inner, mean_outer, has_junc, ix, _id, seq_len, junc_seq))

                if max_obs is not None and analysis_counter.accepted() >= sampled_obs[n]:
                    break

            if max_obs is not None and analysis_counter.accepted() >= max_obs:
                raise MaxObsLimit

    except MaxObsLimit:
        pass
    finally:
        if progress is not None:
            progress.close()
            progress = None

    if max_obs is not None:
        if analysis_counter.count('accepted') >= max_obs:
            logger.info('Reached requested observation limit [{}]'.format(max_obs))
        else:
            logger.warning('Failed to collect requested number of observations [{}]'.format(max_obs))

    if analysis_counter.fraction('rejected') > 0.9:
        logger.warning('More than 90% of reads were rejected')

    #
    # Reporting
    #
    report = OrderedDict({
        'mode': 'kmer',
        'runtime_info': runtime_info(),
        'input_args': {'kmer_db': kmer_db,
                       'kmer_size': k_size,
                       'read_list': read_files,
                       'enzymes': enzyme_names,
                       'seed': seed,
                       'sample_rate': sample_rate,
                       'max_coverage': max_coverage,
                       'mean_insert': mean_insert,
                       'max_freq_quantile': max_freq_quantile,
                       'max_obs': max_obs,
                       'num_sample': num_sample,
                       'frac_sample': frac_sample
                       },
        'n_parsed_reads': analysis_counter.counts['all'],
        'n_analysed_reads': analysis_counter.analysed(),
        'n_too_short': analysis_counter.count('short'),
        'n_no_flank': analysis_counter.count('flank'),
        'n_ambiguous': analysis_counter.count('ambig'),
        'n_high_cov': analysis_counter.count('high_cov'),
        'n_low_cov': analysis_counter.count('low_cov'),
        'n_zero_cov': analysis_counter.count('zero_cov'),
    })

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
    logger.info('Number of reads filtered [low coverage]: {:,} ({:#.4g}% of analysed)'
                .format(analysis_counter.count('low_cov'),
                        analysis_counter.fraction('low_cov') * 100))
    logger.warning('There were {:,} ({:#.4g}% of analyzed) reads filtered due to gaps in k-mer coverage.'
                   .format(analysis_counter.count('zero_cov'),
                           analysis_counter.fraction('zero_cov') * 100))

    # this might not be needed generally.
    if analysis_counter.fraction('zero_cov') > 0.1:
        logger.warning('Excessive gaps in coverage can lead to low estimates. Consider reducing "--min-quality" '
                       'or check that reads and library match')

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
        # make the boolean explicitly categorical so that a count of zero can be returned below
        all_df['has_junc'] = pandas.Categorical(all_df.has_junc, categories=[True, False])

        # remove any row which had zero coverage in inner or outer region
        z_inner = all_df.mean_inner == 0
        z_outer = all_df.mean_outer == 0
        nz_either = ~(z_inner | z_outer)
        all_df = all_df[nz_either]
        logger.debug('Rows removed due to zero coverage: inner {:,}, outer {:,}, shared {:,}'.format(
            sum(z_inner), sum(z_outer), sum(z_inner & z_outer)))

        # inspect the tally of observations with/without junctions
        n_with_junc, n_without_junc = all_df.groupby('has_junc').size().sort_values(ascending=True)
        report['n_without_junc'] = n_without_junc
        report['n_with_junc'] = n_with_junc
        logger.info('Break down of accepted reads: without junction {:,} with junction {:,}'
                    .format(n_without_junc, n_with_junc))
        if n_without_junc == 0:
            raise InsufficientDataException('No reads without the junction sequence were observed.')
        if n_with_junc == 0:
            raise InsufficientDataException('No reads containing the junction sequence were observed.')

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
                    .format(sum(1 / 4 ** ci.site_len for ci in digest.cutsites.values()) * 100))
        logger.info('Number of accepted reads containing potential junctions: {:,} ({:#.4g}% of accepted)'
                    .format(n_with_junc,
                            n_with_junc / analysis_counter.count('accepted') * 100))
        logger.info('Expected fraction by random chance 50% GC: {:#.4g}%'
                    .format(sum(1 / 4 ** li.junc_len for li in digest.junctions.values()) * 100))

        # the total number of read observations
        n_read_obs = len(all_df)
        # the number of unique sequence ids is assumed to be the number of fragment observations
        n_frag_obs = len(set(all_df.seq))
        logger.info('Number of observations: fragments {}, reads {}'.format(n_frag_obs, n_read_obs))

        if n_read_obs == n_frag_obs:
            logger.warning('Equal number of fragments and reads! Hi-C is expected to be paired.')
        obs_redundancy = n_read_obs / n_frag_obs
        logger.info('Fragment observational redundancy: {:.4g}'.format(obs_redundancy))

        if n_read_obs < 500 or n_frag_obs < 500:
            # Single estimation for a small observation pool
            small_pool = True
            logger.warning('The number of observations was small, bootstrapping will not be performed')
            _est = pvalue_expectation(assign_empirical_pvalues_all(all_df))
            hic_frac = np.array([_est['mean'] - _est['error'], _est['mean'] + _est['error']])
            logger.info('Estimated raw Hi-C read fraction from a single sample via '
                        'p-value sum method: {:#.4g} - {:#.4g} %'
                        .format(*hic_frac * 100))

        else:
            # init same random number state
            random_state = np.random.RandomState(seed)
            # Use bootstrap resampling to estimate confidence interval
            small_pool = False
            _sample_est = []
            for _ in tqdm.tqdm(range(num_sample), desc='Bootstrapping CI'):
                # calculate p-values using a sample of observations
                _smpl = assign_empirical_pvalues_all(extract_sample(all_df, frac_sample, 100, 100000, random_state))
                # estimate hi-c fraction for this sample
                _sample_est.append(pvalue_expectation(_smpl))
            _sample_est = pandas.DataFrame(_sample_est)

            # CI defined by 95% quantile range
            hic_frac = _sample_est['mean'].quantile(q=[0.025, 0.975]).values
            logger.info('Estimated raw Hi-C read fraction from bootstrap resampling '
                        'via p-value sum method: 95% CI [{:#.4g},{:#.4g}] %'
                        .format(*hic_frac * 100))

        report['raw_fraction'] = hic_frac

        # TODO for multi-digests with varying length ligation products, using the longest
        #   will lead to overestimating the observered fraction. This could be calculated as
        #   a weighted sum over abundance of each possible junction.
        obs_frac = observed_fraction(round(mean_read_len), round(mean_insert), 'additive',
                                     k_size, digest.longest_junction())
        if 1 - obs_frac < 0:
            logger.warning('Small fragment size can lead to double-counting due to read-pair overlap')
            logger.warning('Observed fraction exceeds 1 ({:#.4g}), therefore unobserved will be reported as 0'
                           .format(obs_frac))
            report['unobs_fraction'] = 0
        else:
            logger.info('For supplied insert length of {:.0f}nt, estimated unobserved fraction: {:#.4g}'
                        .format(mean_insert, 1 - obs_frac))
            report['unobs_fraction'] = 1 - obs_frac

        report['adj_fraction'] = hic_frac * 1 / obs_frac
        if np.any(report['adj_fraction'] > 1):
            logger.warning('Rejecting nonsensical result for adjusted fraction that exceeded 100%')
            report['adj_fraction'] = None
        else:
            if small_pool:
                logger.info('Estimated Hi-C read fraction adjusted for unobserved: {:#.4g} - {:#.4g} %'
                            .format(*hic_frac * 1/obs_frac * 100))
            else:
                logger.info('Estimated Hi-C read fraction adjusted for unobserved: 95% CI [{:#.4g},{:#.4g}] %'
                            .format(*hic_frac * 1/obs_frac * 100))

        report['digestion'] = {'cutsites': digest.to_dict('cutsite'),
                               'junctions': digest.to_dict('junction')}

        report['junction_frequency'] = {}
        for _e, _counts in digest.gather_tracker(junction_tracker).items():
            for _j, _n in _counts.items():
                logger.info('For {} junction sequence {} found: {}'.format(_e, _j, _n))
                report['junction_frequency']['{} {}'.format(_e, _j)] = _n

        # combine them together
        if output_table is not None:
            import os
            logger.info('Writing observations to gzipped tab-delimited file: {}'.format(output_table))
            with gzip.open(output_table, 'wt') as out_h:
                all_df.to_csv(out_h, sep='\t')

    except InsufficientDataException as e:
        logger.warning(e)

    finally:
        if not no_json:
            write_jsonline(report_path, report, suffix='json', append=False)
        if not no_html:
            write_html_report(report_path, report)
