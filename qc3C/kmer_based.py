import pandas
import tqdm

from collections import OrderedDict
from recordclass import make_dataclass, asdict
from scipy.stats import gmean
from typing import List, Tuple
from qc3C.exceptions import *
from qc3C.ligation import Digest
from qc3C.utils import *
from qc3C._version import runtime_info

try:
    import dna_jellyfish
except ImportError:
    import jellyfish as dna_jellyfish

logger = logging.getLogger(__name__)

CovInfo = make_dataclass('cov_info',
                         [('mean_inner', float), ('mean_outer', float), ('has_junc', bool), ('read_pos', int),
                          ('seq_id', str), ('read_extent', int), ('junc_size', int), ('obs_count', int),
                          ('is_merged', bool)])

OBS_ARRAY_DTYPE = np.dtype([('ratio', np.float64), ('has_junc', np.uint8),
                            ('seq_hash', np.int64), ('read_extent', np.uint64),
                            ('junc_size', np.int64), ('is_merged', np.bool),
                            ('obs_count', np.int64)])


@nb.jit(nopython=True, cache=True)
def obs_fraction(obs, insert_size):
    """
    Calculate the fraction of insert extent that has been observed in this set of observations.
    For merged reads, the extent is considered full observed, while we compare the cumulative extent
    of unmerged pairs against the mean insert size.
    """
    frac = np.zeros(obs.shape[0], np.float64)
    for i in range(obs.shape[0]):
        _oi = obs[i]
        if _oi['is_merged']:
            _mask = make_observable_mask(_oi['read_extent'], _oi['read_extent'], False, 0, _oi['junc_size'])
        else:
            _mask = make_observable_mask(int(_oi['read_extent'] / _oi['obs_count']), insert_size,
                                         _oi['obs_count'] > 1, 0, _oi['junc_size'])

        # binary coverage of the ith fragment
        frac[i] = (_mask > 0).mean()

    # overall coverge summary statistics for this set of observations
    return frac.mean()

@nb.jit(nopython=True, cache=True)
def empirical_fraction(obs, insert_size):
    """
    Calculate an estimate for the empirical fraction from an array of observations. The array contains
    k-mer frequency observations for both those which had junctions and those which did not. Those
    without define our empirical p-values for the null hypothesis.

    :param obs: the table of observations
    :param insert_size: expected insert size
    :return: an scalar estimates for the empirical Hi-C fraction, its variance and the observed fraction
    """
    # separate obs with and without junctions
    wj_ratios = np.sort(obs[obs['has_junc'] == 1]['ratio'])
    nj_ratios = np.sort(obs[obs['has_junc'] == 0]['ratio'])

    # empirical p-values are generated from only unique (no junc) observations
    uniq_ratios = np.unique(nj_ratios)
    max_n = len(uniq_ratios)
    emp_pvals = np.linspace(0, 1, max_n + 3)[2:-1]

    # assign p-values to obs containing junctions, where ratio dictates p-value
    wj_pvals = np.zeros_like(wj_ratios)
    n = 0
    for i in range(len(uniq_ratios)):
        while n < len(wj_ratios) and wj_ratios[n] <= uniq_ratios[i]:
            wj_pvals[n] = emp_pvals[i]
            n += 1
        if n >= len(wj_ratios):
            break

    # for cases when we have additional with-junc obs,
    # these get assigned the largest value we have at our disposal
    for i in range(n, len(wj_ratios)):
        wj_pvals[i] = emp_pvals[-1]

    # n_fragments = len(np.unique(obs['seq_hash']))
    n_fragments = len(obs)

    q = 1 - wj_pvals
    return q.sum() / n_fragments, np.sqrt((q * (1-q)).sum()) / n_fragments, obs_fraction(obs, insert_size)


def dataframe_to_array(_df: pandas.DataFrame):
    """
    Extract a numpy array of the necessary columns for estimating the observed fraction
    :param _df: dataframe of observations
    :return: numpy array
    """
    _obs_arr = _df.loc[:, ['ratio', 'has_junc',
                           'seq_id', 'read_extent',
                           'junc_size', 'is_merged',
                           'obs_count']].values
    _obs_arr = np.array([(_obs_arr[i, 0], _obs_arr[i, 1],
                          hash(_obs_arr[i, 2]), _obs_arr[i, 3],
                          _obs_arr[i, 4], _obs_arr[i, 5],
                          _obs_arr[i, 6])
                         for i in range(_obs_arr.shape[0])], OBS_ARRAY_DTYPE)
    return _obs_arr


class BootstrapFraction(object):

    def __init__(self, _df: pandas.DataFrame, _n_samples: int, _frac: float,
                 _insert_size: int, _min_size: int, _seed: int):
        """
        A bootstrapping analysis of the set of junction k-mer frequency observations.

        :param _df: the table of observations, both with and without junctions
        :param _n_samples: the number of bootstrap samples to use in estimation
        :param _frac: the fraction of the total observations to use per bootstrap
        :param _insert_size: the expected size of inserts
        :param _min_size: the minimum sample size
        :param _seed: a random seed.
        """

        self.obs_arr = dataframe_to_array(_df)
        self.num_obs = len(self.obs_arr)
        self.insert_size = _insert_size
        self.seed = _seed
        self.sample_size = int(self.num_obs * _frac)
        if self.sample_size < _min_size:
            logger.debug('Observation pool is small, limiting sample size to {}'.format(_min_size))
            self.sample_size = _min_size
        self.random_state = np.random.RandomState(_seed)

        # compute the bootstrap estimates, storing in an array
        self.estimates = np.zeros((_n_samples, 3), dtype=np.float)
        self.extents = np.zeros(_n_samples, dtype=np.int)
        for i in tqdm.tqdm(range(_n_samples), desc='Bootstrapping CI'):
            _smpl = self._take_sample()
            self.extents[i] = _smpl['read_extent'].sum()
            self.estimates[i, :] = empirical_fraction(_smpl, self.insert_size)

    def _take_sample(self):
        """
        Take a sample from the pool of observations using the sampling
        characteristics defined at instantiation.
        :return: the sample array
        """
        return self.obs_arr[self.random_state.choice(self.num_obs, size=self.sample_size, replace=True)]

    def summary(self) -> Tuple[float, float, float]:
        """
        Basic summary statistics for the computed set of estimates.
        :return: mean, std, median
        """
        return self.estimates[:, 0].mean(), self.estimates[:, 0].std(), np.median(self.estimates[:, 0])

    def confidence(self, conf_width: float = 0.95) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Confidence interval for the set of computed estimates, returning both raw and adjusted estimations
        :param conf_width: confidence interval width (default 95%)
        :return: low and high quantiles for both raw and adjusted for unobserved extent
        """
        assert 0 < conf_width < 1, 'Confidence interval width must be between 0 and 1.'
        _q = 0.5 * (1 - conf_width)
        _raw = np.quantile(self.estimates[:, 0], q=[_q, 1-_q])
        _obs_frac = np.quantile(self.estimates[:, 2], q=[_q, 1-_q])
        _adj = _raw * 1 / _obs_frac
        return _raw, _adj, _obs_frac

    def obs_frac(self) -> float:
        """
        Calculate the mean observed extent over the samples
        :return: observed fraction
        """
        return self.estimates[:, 2].mean()


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


def get_kmer_frequency_summary(filename: str, q: float = 0.9, threads: int = 1) -> Tuple[int, float]:
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
        kmer_histo = pandas.read_csv(proc.stdout, sep=' ', names=['cov', 'count'])
        max_cov = kmer_histo['cov'].quantile(q=q)
        mean_cov = (kmer_histo['cov'] * kmer_histo['count']).sum() / kmer_histo['count'].sum()
        return round(float(max_cov)), mean_cov
    except Exception as e:
        raise RuntimeError('Encountered a problem while analyzing kmer distribution. [{}]'.format(e))


def next_read(filename: str, searcher, longest_site: int, k_size: int,
              prob_accept: float, progress, counts: Counter,
              tracker: Dict[Digest.DerivedString, int], no_user_limit: bool,
              random_state: np.random.RandomState) -> Tuple[str,
                                                            Optional[int],
                                                            str,
                                                            int,
                                                            bool,
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
    :return: yield tuple of (sequence, site position, read id, seq length, is_merged, junc length, junc seq)
    """
    reader = read_seq(open_input(filename))
    unif = random_state.uniform

    for _id, _desc, _seq, _qual in reader:

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

        # simple test for fastp merged reads
        _is_merged = 'merged_' in _desc

        # TODO this (and previous approach) only matches the first instance of the site within the read.
        #  There is no reason that this would be the correct location.
        _match = searcher(_seq)
        # no site found
        if _match is None:
            yield _seq, None, _id, seq_len, _is_merged, None, None
        else:
            _ix, _end = _match.span()
            _junc_len = _end - _ix
            _junc_seq = _match.group()
            tracker[_junc_seq] += 1
            yield _seq, _ix, _id, seq_len, _is_merged, _junc_len, _junc_seq


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

    def collect_coverage(seq: str, ix: int, site_size: int, k: int, min_cov: float) -> (float, float):
        """
        Collect the k-mer coverage centered around the position ix. From the left, the sliding
        window begins just before the site region and slides right until just after. Means
        are then calculated for the inner (within the junction) and outer (left and right flanks)
        :param seq: the sequence to analyse
        :param ix: the position marking the beginning of a junction or any other arbitrary location if so desired
        :param site_size: the size of the junction site
        :param k: the kmer size
        :param min_cov: minimum accepted coverage for outer
        :return: mean(inner), mean(outer)
        """
        def get_kmer_cov(mer: str) -> float:
            mer = dna_jellyfish.MerDNA(mer)
            mer.canonicalize()
            k_cov = query_jf[mer]
            return k_cov

        r_shift = 0
        l_flank = r_flank = (k - site_size) // 2
        if (k - site_size) % 2 != 0:
            l_flank += 1
            r_shift = 1

        if ix <= k:
            outer = np.zeros(shape=3, dtype=np.int)
            inner = np.zeros(shape=3, dtype=np.int)
            for i in range(-1, 2):
                si = seq[ix + site_size + r_shift + i: ix + site_size + k + r_shift + i]
                outer[i+1] = get_kmer_cov(si)
                si = seq[ix - l_flank + i: ix + site_size + r_flank + i]
                inner[i+1] = get_kmer_cov(si)
                if outer[i+1] == 0 or inner[i+1] == 0:
                    raise ZeroCoverageException()

        elif ix + site_size + r_shift + k -1 >= len(seq):
            outer = np.zeros(shape=3, dtype=np.int)
            inner = np.zeros(shape=3, dtype=np.int)
            for i in range(-1, 2):
                si = seq[ix - k + i: ix + i]
                outer[i+1] = get_kmer_cov(si)
                si = seq[ix - l_flank + i: ix + site_size + r_flank + i]
                inner[i+1] = get_kmer_cov(si)
                if outer[i+1] == 0 or inner[i+1] == 0:
                    raise ZeroCoverageException()

        else:
            outer = np.zeros(shape=6, dtype=np.int)
            inner = np.zeros(shape=3, dtype=np.int)
            for i in range(-1, 2):
                si = seq[ix - k + i: ix + i]
                outer[i+1] = get_kmer_cov(si)
                si = seq[ix + site_size + r_shift + i: ix + site_size + k + r_shift + i]
                outer[i+4] = get_kmer_cov(si)
                si = seq[ix - l_flank + i: ix + site_size + r_flank + i]
                inner[i+1] = get_kmer_cov(si)
                if outer[i+1] == 0 or outer[i+4] == 0 or inner[i+1] == 0:
                    raise ZeroCoverageException()

        cov_inner, cov_outer = gmean(inner), gmean(outer)
        if cov_outer < min_cov:
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

    max_coverage, mean_coverage = get_kmer_frequency_summary(kmer_db, max_freq_quantile, threads)
    logger.info('Maximum k-mer frequency of {} at quantile {} in Jellyfish library {}'
                .format(max_coverage, max_freq_quantile, kmer_db))
    # A rule of thumb min threshold
    min_coverage = mean_coverage / 2.5
    if min_coverage < 1:
        # no less than 1, which is pedantic as later we don't accept zero-coverage k-mers anyhow
        min_coverage = 1
    logger.info('Mean k-mer frequency from Jellyfish library: {:.2f}'.format(mean_coverage))
    logger.info('Minimum acceptable local k-mer frequency: {:.2f}'.format(min_coverage))

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
    coverage_obs = {}

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
                    seq, ix, _id, seq_len, is_merged, junc_len, junc_seq = next(fq_reader)
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

                    # randomly select a region which does not contain ambiguous bases,
                    # but give-up after a few failed attempts
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
                    mean_inner, mean_outer = collect_coverage(seq, ix, junc_len, k_size, min_coverage)
                except ZeroCoverageException:
                    analysis_counter.counts['zero_cov'] += 1
                    continue
                except LowCoverageException:
                    analysis_counter.counts['low_cov'] += 1
                    continue

                # avoid regions with strangely high coverage
                if mean_inner > max_coverage or mean_outer > max_coverage:
                    analysis_counter.counts['high_cov'] += 1
                    continue

                if user_limited:
                    # user limited runs get progress updates here rather than
                    # within the read iterator
                    progress.update()

                cumulative_length += seq_len

                # if this is a new fragment, record the observation
                if _id not in coverage_obs:
                    coverage_obs[_id] = CovInfo(mean_inner, mean_outer, has_junc, ix,
                                                _id, seq_len, junc_len, 1, is_merged)
                # otherwise always supersede preexisting non-junc obs
                # or for junc obs flip a coin to choose which obs to keep
                else:
                    # we must always record the additional read and flank extent
                    _ci = coverage_obs[_id]
                    assert not _ci.is_merged, 'Merged sequence appears more than once!'
                    _ci.read_extent += seq_len
                    _ci.obs_count += 1

                    # when both do or do not contain a junction, flip a coin to accept one obs over the other
                    if ((_ci.has_junc and has_junc) or (not _ci.has_junc and not has_junc)) and randint(2) == 0:
                        # replace obs
                        coverage_obs[_id] = CovInfo(mean_inner, mean_outer, has_junc, ix, _id,
                                                    _ci.read_extent, junc_len, _ci.obs_count, is_merged)

                    # always take a junction observation over non-junc
                    # this emulates what would occur if we could see the entire insert each time.
                    elif has_junc:
                        # replace obs
                        coverage_obs[_id] = CovInfo(mean_inner, mean_outer, has_junc, ix, _id,
                                                    _ci.read_extent, junc_len, _ci.obs_count, is_merged)

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
                       'mean_coverage': mean_coverage,
                       'min_coverage': min_coverage,
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
    logger.warning('Number of reads filtered [zero coverage]: {:,} ({:#.4g}% of analyzed)'
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
        # all_df = pandas.DataFrame(np.hstack(coverage_obs.values()))
        all_df = pandas.DataFrame(map(asdict, coverage_obs.values()))
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
        n_with_junc = (all_df['has_junc'] == True).sum()
        n_without_junc = len(all_df) - n_with_junc
        report['n_without_junc'] = n_without_junc
        report['n_with_junc'] = n_with_junc
        logger.info('Break down of fragment observations: without junction {:,} with junction {:,}'
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
        n_read_obs = all_df.obs_count.sum()
        # the number of unique sequence ids is assumed to be the number of fragment observations
        n_frag_obs = len(all_df)
        logger.info('Number of observations: fragments {}, reads {}'.format(n_frag_obs, n_read_obs))

        # if n_read_obs == n_frag_obs:
        #     logger.warning('Equal number of fragments and reads! Hi-C is expected to be paired.')
        obs_redundancy = n_read_obs / n_frag_obs
        logger.info('Fragment observational redundancy (read obs per fragment): {:.4g}'.format(obs_redundancy))

        # TODO reinstate the small pool code
        if n_frag_obs < 500:
            raise ApplicationException('At least {} fragment observations are required. Found {}'.format(
                500, n_frag_obs))
        # if n_read_obs < 500 or n_frag_obs < 500:
        #     # Single estimation for a small observation pool
        #     small_pool = True
        #     logger.warning('The number of observations was small, bootstrapping will not be performed')
        #     _mean, _err = empirical_fraction(dataframe_to_array(all_df))
        #     raw_frac = np.array([_mean - 2*_err, _mean + 2*_err])
        #     logger.info('Estimated raw Hi-C read fraction from a single sample via '
        #                 'p-value sum method: [{:#.4g},{:#.4g}] %'
        #                 .format(*hic_frac * 100))
        #     obs_frac = all_df.obs_extent.sum() / (mean_insert * n_frag_obs)
        # else:

        small_pool = False

        # Use bootstrap resampling to estimate confidence interval for raw and adjusted
        bs_frac = BootstrapFraction(all_df, num_sample, frac_sample, mean_insert, 100, seed)
        raw_frac, adj_frac, obs_frac = bs_frac.confidence()
        # Also report what the observed fraction was estimated to be
        obs_frac = obs_frac.mean()

        # TODO for multi-digests with varying length ligation products, using the longest
        #   will lead to overestimating the observed fraction. This could be calculated as
        #   a weighted sum over abundance of each possible junction.
        if 1 - obs_frac < 0:
            logger.warning('Small fragment size can lead to double-counting due to read-pair overlap')
            logger.warning('Observed fraction exceeds 1 ({:#.4g}), therefore unobserved will be reported as 0'
                           .format(obs_frac))
            report['unobs_fraction'] = 0
        else:
            logger.info('For supplied insert length of {:.0f}nt, estimated unobserved fraction: {:#.4g}'
                        .format(mean_insert, 1 - obs_frac))
            report['unobs_fraction'] = 1 - obs_frac

        report['raw_fraction'] = raw_frac
        logger.info('Estimated raw Hi-C read fraction from bootstrap resampling '
                    'via p-value sum method: 95% CI [{:#.4g},{:#.4g}] %'
                    .format(*raw_frac * 100))

        report['adj_fraction'] = adj_frac
        if np.any(report['adj_fraction'] > 1):
            logger.warning('Rejecting nonsensical result for adjusted fraction that exceeded 100%')
            report['adj_fraction'] = None
        else:
            if small_pool:
                logger.info('Estimated Hi-C read fraction adjusted for unobserved: {:#.4g} - {:#.4g} %'
                            .format(*adj_frac * 100))
            else:
                logger.info('Estimated Hi-C read fraction adjusted for unobserved: 95% CI [{:#.4g},{:#.4g}] %'
                            .format(*adj_frac * 100))

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
