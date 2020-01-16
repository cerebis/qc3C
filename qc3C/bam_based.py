import logging
import multiprocessing
import numpy as np
import os
import pysam
import subprocess
import tqdm

from typing import Optional, Tuple
from abc import ABC, abstractmethod
from qc3C.exceptions import NameSortingException, UnknownLibraryKitException, MaxObsLimit, InsufficientDataException
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance, LigationInfo, CutSitesDB
from qc3C.utils import init_random_state, warn_if, write_jsonline, observed_fraction
from qc3C._version import runtime_info

logger = logging.getLogger(__name__)

# translation table used for complementation
COMPLEMENT_TABLE = str.maketrans('acgtumrwsykvhdbnACGTUMRWSYKVHDBN',
                                 'TGCAAnnnnnnnnnnnTGCAANNNNNNNNNNN')


def revcomp(seq: str) -> str:
    """
    Reverse complement a string representation of a sequence. This uses string.translate.
    :param seq: input sequence as a string
    :return: revcomp sequence as a string
    """
    return seq.translate(COMPLEMENT_TABLE)[::-1]


def exe_exists(exe_name: str) -> bool:
    """
    Check that a executable exists on the Path.
    :param exe_name: the base executable name
    :return: True, an executable file named exe_name exists and has executable bit set
    """
    p, f = os.path.split(exe_name)
    assert not p, 'include only the base file name, no path specification'

    for pn in os.environ["PATH"].split(':'):
        full_path = os.path.join(pn, exe_name)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return True
    return False


def count_bam_reads(file_name: str, paired: bool = False, mapped: bool = False,
                    mapq: int = None, max_cpu: int = None) -> int:
    """
    Use samtools to quickly count the number of non-header lines in a bam file. This is assumed to equal
    the number of mapped reads.
    :param file_name: a bam file to scan (neither sorted nor an index is required)
    :param paired: when True, count only reads of mapped pairs (primary alignments only)
    :param mapped: when True, count only reads which are mapped
    :param mapq: count only reads with mapping quality greater than this value
    :param max_cpu: set the maximum number of CPUS to use in counting (otherwise all cores)
    :return: estimated number of mapped reads
    """
    assert exe_exists('samtools'), 'required tool samtools was not found on path'
    assert exe_exists('wc'), 'required tool wc was not found on path'
    assert not (paired and mapped), 'Cannot set paired and mapped simultaneously'

    if not os.path.exists(file_name):
        raise IOError('{} does not exist'.format(file_name))
    if not os.path.isfile(file_name):
        raise IOError('{} is not a file'.format(file_name))

    opts = ['samtools', 'view', '-c']
    if max_cpu is None:
        max_cpu = multiprocessing.cpu_count()
    opts.append('-@{}'.format(max_cpu))

    if paired:
        opts.append('-F0xC0C')
    elif mapped:
        opts.append('-F0xC00')

    if mapq is not None:
        assert 0 <= mapq <= 60, 'mapq must be in the range [0,60]'
        opts.append('-q{}'.format(mapq))

    proc = subprocess.Popen(opts + [file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    value_txt = proc.stdout.readline().strip()
    try:
        count = int(value_txt)
        if paired and count % 2 != 0:
            logger.warning('When counting paired reads, the total was not divisible by 2.')
        return count
    except ValueError:
        raise RuntimeError('Encountered a problem determining alignment count. samtools returned [{}]'
                           .format(value_txt))


def greedy_pair_separation(r1: pysam.AlignedSegment, r1ref: int,
                           r2: pysam.AlignedSegment, r2ref: int,
                           min_reflen: int, cutsite_db: CutSitesDB) -> Optional[Tuple[int, int]]:
    """
    Greedily try to estimate a distance between pairs, even when they map to different references. The intervening
    distance is assumed to be the least possible (summing distances of reads from ends of their references)
    :param r1: read 1
    :param r1ref: the length of reference to which read 1 is mapped
    :param r2: read 2
    :param r2ref: the length of reference to which read 2 is mapped
    :param min_reflen: the shortest acceptable reference length with which to be greedy
    :param cutsite_db: an instance of CutSiteDB pertaining to the set of references and enzyme in question.
    :return: an integer distance and number of cutsites or None
    """

    def nearest_edge_details(read, ref_len):
        """
        Determine the nearest contig edge from the read's position and count the number of cutsites
        which exists between the edge and the read's position.
        :param read: the read to assess
        :param ref_len: the length of the mapped reference
        :return: tuple of distance from edge and number of cutsites
        """
        left_d = read.pos
        right_d = ref_len - (read.pos + read.alen)
        if left_d < right_d:
            return left_d, cutsite_db[read.reference_name].count_sites(0, left_d)
        else:
            return right_d, cutsite_db[read.reference_name].count_sites(right_d, ref_len)

    assert r1.reference_id != r2.reference_id, 'references must be on different strands for greedy estimation'

    # only do this for sufficiently long contigs
    if r1ref < min_reflen or r2ref < min_reflen:
        return None

    if r1.reference_start < 1000 or (r1ref - r1.reference_end) < 1000 or \
            r2.reference_start < 1000 or (r2ref - r2.reference_end) < 1000:
        return None

    # estimate a distance by looking at the the placement of reads
    # within their respective contigs. Assume the least possible
    # intervening distance.
    r1_d, r1_ns = nearest_edge_details(r1, r1ref)
    r2_d, r2_ns = nearest_edge_details(r2, r2ref)

    # sum edge distances and cutsites
    d = r1_d + r2_d
    ns = r1_ns + r2_ns

    # for greedy pairs, we require that at least one read map
    # the full distance from its nearest end, otherwise return None
    if r1_d > 10000 or r2_d > 10000:
        return d, ns
    elif r1_d > 5000 or r2_d > 5000:
        return d, ns
    elif r1_d > 1000 or r2_d > 1000:
        return d, ns
    return None


class Counter(ABC):

    def __init__(self, counts: dict):
        self._counts = {'all': 0, 'unmapped': 0, 'sample': 0, 'accepted': 0}
        self._counts.update(counts)

    @abstractmethod
    def accept(self, **kwargs: pysam.AlignedSegment) -> bool:
        pass

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
        else:
            return self._counts[category]

    def fraction(self, category: str, against: str = 'analysed') -> float:
        """
        Return the category's fraction compared to one of (accepted, all, analysed).

        :param category: the category (numerator)
        :param against: the denominator
        :return: a fraction [0,1]
        """
        if against == 'accepted':
            return self.count(category) / self._counts['accepted']
        elif against == 'all':
            return self.count(category) / self._counts['all']
        elif against == 'analysed':
            return self.count(category) / self.analysed()
        else:
            raise RuntimeError('parameter \"against\" must be one of [accepted, analysed or all]')

    def analysed(self) -> int:
        """
        The number of reads/pairs analysed is all items minus those skipped from sub-sampling
        or which were unmapped.

        :return: The number of items analysed
        """
        return self._counts['all'] - self._counts['sample'] - self._counts['unmapped']

    def rejected(self) -> int:
        """
        :return: The number of reads/pairs which were rejected due to specified constraints.
        """
        return self.analysed() - self._counts['accepted']


class ReadFilter(Counter):

    def __init__(self, sample_rate: float, min_mapq: int, min_reflen: int, min_match: int,
                 no_secondary: bool, no_supplementary: bool, no_refterm: bool,
                 ref_lengths, random_state: np.random.RandomState = None):
        """
        Filter reads based on a number of criteria, maintaining counts of the various categorisations.

        :param sample_rate: consider only a fraction of all reads [0..1]
        :param min_mapq: the minimum acceptable map quality score
        :param min_reflen: the shortest acceptable reference sequence length
        :param min_match: match at least this many nt from the start from the read
        :param no_secondary: reject secondary alignments
        :param no_supplementary: reject supplementary alignments
        :param no_refterm: reject alignments which terminate early due to reaching the end of the reference
        :param ref_lengths: a list of all reference lengths
        :param random_state: an optional random state required for sub-sampling
        """
        super().__init__({'mapq': 0, 'sample': 0, 'ref_len': 0, 'secondary': 0,
                          'supplementary': 0, 'weak': 0, 'ref_term': 0})
        self.sample_rate = sample_rate
        self.min_mapq = min_mapq
        self.min_reflen = min_reflen
        self.min_match = min_match
        self.no_secondary = no_secondary
        self.no_supplementary = no_supplementary
        self.no_refterm = no_refterm
        self.ref_lengths = ref_lengths
        if random_state is not None:
            self.unif = random_state.uniform

    def accept(self, r: pysam.AlignedSegment) -> bool:
        """
        Test if a read passes the specified criteria.

        :param r: the read to test
        :return: True if the read is accepted
        """
        self._counts['all'] += 1

        #
        # immediately stop from some conditions
        #

        if r.is_unmapped:
            self._counts['unmapped'] += 1
            return False

        # for sampling, don't continue to analyse, just return
        if self.sample_rate is not None and self.sample_rate < self.unif():
            self._counts['sample'] += 1
            return False

        #
        # test for each condition, otherwise
        #
        _accept = True
        if self.min_mapq is not None and r.mapping_quality < self.min_mapq:
            self._counts['mapq'] += 1
            _accept = False
        if self.min_reflen is not None and self.ref_lengths[r.reference_id] < self.min_reflen:
            self._counts['ref_len'] += 1
            _accept = False
        if self.no_secondary and r.is_secondary:
            self._counts['secondary'] += 1
            _accept = False
        if self.no_supplementary and r.is_supplementary:
            self._counts['supplementary'] += 1
            _accept = False
        if self.min_match is not None:
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            if cig[0] != 0 or cig[1] < self.min_match:
                self._counts['weak'] += 1
                _accept = False
        if self.no_refterm and r.query_alignment_length < r.query_length and \
                (r.reference_end >= self.ref_lengths[r.reference_id] or r.reference_start == 0):
            self._counts['ref_term'] += 1
            _accept = False

        if _accept:
            self._counts['accepted'] += 1

        return _accept


class PairFilter(Counter):

    def __init__(self, no_trans: bool = False):
        """
        Filter paired reads using specified criteria. Keeping track of categorisation counts.

        :param no_trans: reject trans-mapping pairs (reads that map to different references)
        """
        super().__init__({'trans': 0})
        self.no_trans = no_trans

    def count(self, category: str) -> int:
        """
        Return the count of a category

        :param category:
        :return: the count
        """
        if category == 'cis':
            return self._counts['all'] - self._counts['trans']
        else:
            return super().count(category)

    def accept(self, r1: pysam.AlignedSegment, r2: pysam.AlignedSegment) -> bool:
        """
        Test of a pair passes the specified criteria

        :param r1: the first read
        :param r2: the second read
        :return: True if pair is accepted
        """
        _accept = True
        self._counts['all'] += 1
        if r1.reference_id != r2.reference_id:
            # always count trans, even if not emitted
            self._counts['trans'] += 1
            if self.no_trans:
                _accept = False
        if _accept:
            self._counts['accepted'] += 1
        return _accept


class read_pairs(object):

    def __init__(self, bam_path: str, random_state: np.random.RandomState = None, sample_rate: float = None,
                 min_mapq: int = None, min_reflen: int = None, min_match: int = None, no_trans: bool = False,
                 no_secondary: bool = False, no_supplementary: bool = False, no_refterm: bool = False,
                 threads: int = 1, count_reads: bool = False, show_progress: bool = False,
                 is_ipynb: bool = False):
        """
        An iterator over the pairs in a bam file, with support for the 'with' statement.
        Pairs are filtered based on user specified critera. Only pairs, whose read both
        pass these criteria are emitted. Filteration is done both on the read level and
        as pairs.

        :param bam_path: the path to the bam file
        :param random_state: a numpy random state required if sub-sampling
        :param sample_rate: consider only a fraction of all reads [0..1]
        :param min_mapq: the minimum acceptable mapping quality
        :param min_reflen: the minimum acceptable reference length
        :param min_match: the number of nt which must align, beginning from the start of a read
        :param no_trans: reject trans-mapping pairs
        :param no_secondary: reject secondary alignments
        :param no_supplementary: reject supplementary alignments
        :param no_refterm: reject reads which terminate early due to reaching the end of a reference
        :param threads: the number of concurrent threads when accessing the bam file
        :param count_reads: before starting, count the number of alignments in the bam file
        :param show_progress: display progress using tqdm
        :param is_ipynb: use tqdm widget for ipython notepad
        """

        if random_state is None:
            assert sample_rate is None or sample_rate == 1, 'A random state must be supplied when sub-sampling reads'

        self.n_reads = None
        if count_reads:
            logger.info('Counting alignments in {}'.format(bam_path))
            self.n_reads = count_bam_reads(bam_path, max_cpu=threads)
            logger.info('Found {:,} alignments to analyse'.format(self.n_reads))

        self.bam = pysam.AlignmentFile(bam_path, mode='rb', threads=threads)
        _header = self.bam.header['HD']
        if 'SO' not in _header or _header['SO'] != 'queryname':
            raise NameSortingException(bam_path)

        self.random_state = random_state
        self.show_progress = show_progress
        self.progress = None
        self.is_ipynb = is_ipynb
        self.reference_lengths = np.array([li for li in self.bam.lengths])

        # applying sample_rate to reads means that the pass rate of pairs has an upper bound of sample_rate ** 2.
        # therefore take the square, so that the number of pairs is closer to what users expect.
        # this is particularly noticeable for low sample rates (ie passing only 10% of reads means only 1% of pairs)
        pair_sample_rate = None
        if sample_rate is not None:
            pair_sample_rate = sample_rate ** 0.5

        self.read_filter = ReadFilter(pair_sample_rate, min_mapq, min_reflen, min_match,
                                      no_secondary, no_supplementary, no_refterm,
                                      self.reference_lengths, random_state)
        self.pair_filter = PairFilter(no_trans=no_trans)

        self.bam_iter = self.bam.fetch(until_eof=True)
        self.progress = None
        if self.show_progress:
            if self.is_ipynb:
                self.progress = tqdm.tqdm_notebook(total=self.n_reads, desc='Reads')
            else:
                self.progress = tqdm.tqdm(total=self.n_reads, desc='Reads')

    def __enter__(self):
        return self

    def get_reflen(self, r: pysam.AlignedSegment) -> int:
        """
        Return the complete length of the reference sequence on which a read has been aligned. This method
        assumes the read is mapped.
        :param r: the read in question
        :return: the length of the reference
        """
        return self.reference_lengths[r.reference_id]

    def close(self) -> None:
        if self.bam is not None:
            self.bam.close()
            self.bam = None
        if self.progress is not None:
            self.progress.close()
            self.progress = None

    def __exit__(self, *a):
        self.close()

    def __iter__(self) -> (pysam.AlignedSegment, pysam.AlignedSegment):
        """
        Return the next acceptable pair of reads
        :return: r1, r2
        """
        pb_update = None if self.progress is None else self.progress.update
        test_read = self.read_filter.accept
        test_pair = self.pair_filter.accept
        try:

            read_a = next(self.bam_iter)
            pass_a = test_read(read_a)

            if pb_update is not None:
                pb_update()

            while True:

                read_b = next(self.bam_iter)
                pass_b = test_read(read_b)

                if pb_update is not None:
                    pb_update()

                if pass_a and pass_b and read_a.query_name == read_b.query_name:
                    if test_pair(read_a, read_b):
                        yield read_a, read_b

                read_a = read_b
                pass_a = pass_b

        except StopIteration:
            self.close()


def junction_match_length(seq: str, read: pysam.AlignedSegment, lig_info: LigationInfo) -> int:
    """
    For a given read, whose 3' end goes beyond the end of its mapped alignment, check if the
    immediately following nt match sequence the proximity ligation junction. Return the number of
    nts which matched.
    :param seq: the sequence of the read
    :param read: the read mapping object
    :param lig_info: the ligation details for a given enzyme
    :return: the number of matching nt (to PL junction) after the end of the alignment
    """
    if read.is_reverse:
        i = read.query_length - read.query_alignment_start
    else:
        i = read.query_alignment_end

    j = lig_info.vest_len
    while i < len(seq) and j < lig_info.junc_len and seq[i] == lig_info.junction[j]:
        i += 1
        j += 1
    return j - lig_info.vest_len


def analyse(bam_file: str, fasta_file: str, enzyme_name: str,
            seed: int = None, sample_rate: float = None, min_mapq: int = 60, max_pairs: int = None,
            threads: int = 1, report_path: str = None, library_kit: str = 'generic') -> None:
    """
    analyse a bam file which contains Hi-C read-pairs mapped to a set of reference sequences.
    This method attempts to assess details which indicate the overall strength of the
    Hi-C signal. Here signal is the proportion of proximity ligation pairs relative to
    shotgun pairs within the mapped read-sets.

    :param bam_file: the path to the bam file
    :param fasta_file: reference sequences
    :param enzyme_name: the name of the enzyme used in Hi-C library creation
    :param seed: a random integer seed
    :param sample_rate: consider only a portion of all pairs [0..1]
    :param min_mapq: the minimum acceptable mapping quality
    :param max_pairs: the maximum number of accepted pairs to inspect
    :param threads: the number of threads used in accessing the bam file
    :param report_path: append a report in single-line JSON format to the given path.
    :param library_kit: the type of kit used in producing the library (ie. phase, generic)
    """

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

    enzyme = get_enzyme_instance(enzyme_name)
    ligation_info = ligation_junction_seq(enzyme)

    cutsite_db = CutSitesDB(enzyme, fasta_file, use_cache=True, use_tqdm=True)

    # for efficiency, disable sampling if rate is 1
    if sample_rate == 1 or sample_rate is None:
        logger.info('Accepting all usable reads')
        sample_rate = None
    else:
        logger.info('Acceptance threshold: {:#.3g}'.format(sample_rate))

    random_state = init_random_state(seed)

    if max_pairs is not None:
        show_progress = count_reads = False
    else:
        show_progress = count_reads = True

    mp_progress = None
    if max_pairs is not None:
        mp_progress = tqdm.tqdm(total=max_pairs, desc='Pairs')

    cumulative_length = 0
    short_sum = 0
    short_count = 0
    long_counts = np.zeros(3, dtype=np.int)
    long_bins = (1000, 5000, 10000)
    read_counts = {'full_align': 0, 'early_term': 0, 'no_site': 0}
    class_counts = {'dangling': 0, 'self_circle': 0, 'religation': 0, 'ffrr_invalid': 0,
                    'fr_valid': 0, 'rf_valid': 0, 'ffrr_valid': 0, 'valid': 0, 'invalid': 0}
    enzyme_counts = {'cs_full': 0, 'cs_term': 0, 'cs_start': 0, 'read_thru': 0,
                     'is_split': 0, 'partial_readthru': 0}

    n_accepted_obs = 0

    # begin with a guess for insert length
    short_median_est = 300

    pair_parser = read_pairs(bam_file, random_state=random_state, sample_rate=sample_rate,
                             min_mapq=min_mapq, min_match=1, min_reflen=500,
                             no_trans=False, no_secondary=True, no_supplementary=True, no_refterm=True,
                             threads=threads, count_reads=count_reads, show_progress=show_progress)

    logger.info('Beginning analysis...')

    try:

        for r1, r2 in pair_parser:

            # cis-mapping pairs, lets categorise using the standard nomenclature
            if r1.reference_id == r2.reference_id:

                d, ns, flipped = cutsite_db[r1.reference_name].pair_info(r1, r2)

                if flipped:
                    r1, r2 = r2, r1

                # invalid pair
                if ns == 0:

                    if not r1.is_reverse and r2.is_reverse:
                        class_counts['dangling'] += 1
                    elif r1.is_reverse and not r2.is_reverse:
                        class_counts['self_circle'] += 1
                    else:
                        class_counts['ffrr_invalid'] += 1

                    class_counts['invalid'] += 1

                else:

                    # proper orientation
                    if not r1.is_reverse and r2.is_reverse:
                        if ns == 1:
                            class_counts['religation'] += 1
                        else:
                            class_counts['fr_valid'] += 1
                    elif r1.is_reverse and not r2.is_reverse:
                        class_counts['rf_valid'] += 1
                    else:
                        class_counts['ffrr_valid'] += 1

                    class_counts['valid'] += 1

                # when separation has been determined or estimated
                if d is not None and d > 0:

                    # track the number of pairs observed for the requested separation distances
                    if d >= long_bins[2]:
                        long_counts[2] += 1
                        long_counts[1] += 1
                        long_counts[0] += 1
                    elif d >= long_bins[1]:
                        long_counts[1] += 1
                        long_counts[0] += 1
                    elif d >= long_bins[0]:
                        long_counts[0] += 1

                    # for explicitly cis-mapped
                    elif r1.reference_id == r2.reference_id:

                        # keep a count of those pairs which are within the "WGS" region
                        if 50 < d < 1000:
                            short_sum += d
                            short_count += 1

                            # frugal median estimator
                            if short_median_est > d:
                                short_median_est -= 1
                            elif short_median_est < d:
                                short_median_est += 1

            # extract some QC statistics from the reads in both cis and trans pairs,
            for ri in [r1, r2]:

                rlen = ri.query_length
                assert rlen != 0, 'BAM did not contain sequence information for read: {}'.format(ri.query_name)

                cumulative_length += rlen

                # always consider sequences 5'->3'
                if ri.is_reverse:
                    seq = revcomp(ri.seq)
                    aln_seq = seq[ri.query_length - ri.query_alignment_end: ri.query_length - ri.query_alignment_start]
                else:
                    seq = ri.seq
                    aln_seq = seq[ri.query_alignment_start: ri.query_alignment_end]

                # for fully-aligned reads
                if rlen == ri.query_alignment_length:

                    read_counts['full_align'] += 1

                    # check for cut-site
                    _no_site = True
                    if seq.startswith(ligation_info.cut_site):
                        enzyme_counts['cs_start'] += 1

                    if seq.endswith(ligation_info.vestigial):
                        _no_site = False
                        enzyme_counts['cs_full'] += 1
                        enzyme_counts['cs_term'] += 1

                    if _no_site:
                        read_counts['no_site'] += 1

                # for reads terminating early, we can check for the junction duplication
                else:

                    read_counts['early_term'] += 1

                    # check for cut-site and potentially the junction
                    _no_site = True

                    if seq.startswith(ligation_info.cut_site):
                        enzyme_counts['cs_start'] += 1

                    max_match = ligation_info.junc_len - ligation_info.vest_len

                    if aln_seq.endswith(ligation_info.vestigial):
                        _no_site = False

                        enzyme_counts['cs_term'] += 1

                        # look for partial junctions, depending on how much
                        # flanking sequence a read has available
                        if ri.is_reverse:
                            spare = ri.query_alignment_start
                        else:
                            spare = ri.query_length - ri.query_alignment_end

                        m = junction_match_length(seq, ri, ligation_info)
                        if (m < max_match and m == spare) or m == max_match:
                            enzyme_counts['partial_readthru'] += 1

                        # for full matches, keep a separate tally
                        if m == max_match:
                            enzyme_counts['read_thru'] += 1
                            if ri.has_tag('SA'):
                                enzyme_counts['is_split'] += 1

                    if _no_site:
                        read_counts['no_site'] += 1

            if max_pairs is not None:
                n_accepted_obs += 1
                mp_progress.update()
                if n_accepted_obs >= max_pairs:
                    raise MaxObsLimit

    except MaxObsLimit:
        if mp_progress is not None:
            mp_progress.close()
            logger.info('Reached user-defined observation limit [{}]'.format(max_pairs))

    # some shorthand variables which will be used repeatedly
    n_pairs_accepted = pair_parser.pair_filter.count('accepted')
    n_paired_reads = 2 * n_pairs_accepted
    n_cis_pairs = pair_parser.pair_filter.count('cis')

    #
    # Initial values for report
    #

    report = {
        'mode': 'bam',
        'runtime_info': runtime_info(),
        'input_args': {'bam_file': bam_file,
                       'fasta_file': fasta_file,
                       'enzyme': enzyme_name,
                       'seed': seed,
                       'sample_rate': sample_rate,
                       'min_mapq': min_mapq},
        'n_parsed_reads': pair_parser.read_filter.count('all'),
        'n_analysed_reads': pair_parser.read_filter.analysed(),
        'n_accepted_reads': pair_parser.read_filter.count('accepted'),
        'n_unmapped': pair_parser.read_filter.count('unmapped'),
        'n_low_mapq': pair_parser.read_filter.count('mapq'),
        'n_ref_len': pair_parser.read_filter.count('ref_len'),
        'n_secondary': pair_parser.read_filter.count('secondary'),
        'n_supplementary': pair_parser.read_filter.count('supplementary'),
        'n_weak_mapping': pair_parser.read_filter.count('weak'),
        'n_ref_term': pair_parser.read_filter.count('ref_term'),
        'n_analysed_pairs': pair_parser.pair_filter.analysed(),
        'n_accepted_pairs': n_pairs_accepted,
        'n_trans_pairs': pair_parser.pair_filter.count('trans'),
        'n_cis_pairs': n_cis_pairs,
        'n_fully_aligned': read_counts['full_align'],
        'n_align_term': read_counts['early_term'],
        'n_no_site_end': read_counts['no_site'],
        'n_short_inserts': short_count,

        # HiCPro style classifications
        'classification': {

            'informative': {'fr': class_counts['fr_valid'],
                            'rf': class_counts['rf_valid'],
                            'ffrr': class_counts['ffrr_valid']},

            'uninformative': {'religation': class_counts['religation'],
                              'dangling_ends': class_counts['dangling'],
                              'self_circle': class_counts['self_circle'],
                              'ffrr': class_counts['ffrr_invalid']}
        }
    }

    #
    # Log the results
    #

    logger.info('Number of parsed reads: {:,}'
                .format(pair_parser.read_filter.count('all')))
    logger.info('Number of analysed reads: {:,} ({:#.4g}% of all)'
                .format(pair_parser.read_filter.analysed(),
                        pair_parser.read_filter.fraction('analysed', 'all') * 100))
    logger.info('Number of reads filtered [unmapped]: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('unmapped'),
                        pair_parser.read_filter.fraction('unmapped') * 100))
    logger.info('Number of reads filtered [low mapq]: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('mapq'),
                        pair_parser.read_filter.fraction('mapq') * 100))
    logger.info('Number of reads filtered [ref length]: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('ref_len'),
                        pair_parser.read_filter.fraction('ref_len') * 100))
    logger.info('Number of reads filtered [secondary]: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('secondary'),
                        pair_parser.read_filter.fraction('secondary') * 100))
    logger.info('Number of reads filtered [supplementary]: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('supplementary'),
                        pair_parser.read_filter.fraction('supplementary') * 100))
    logger.info('Number of reads filtered [weak mapping]: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('weak'),
                        pair_parser.read_filter.fraction('weak') * 100))
    logger.info('Number of reads filtered [ref terminated]: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('ref_term'),
                        pair_parser.read_filter.fraction('ref_term') * 100))

    # accepted reads
    logger.info('Number of accepted reads: {:,} ({:#.4g}% of analysed)'
                .format(pair_parser.read_filter.count('accepted'),
                        pair_parser.read_filter.fraction('accepted') * 100))

    # for pairs which have accepted read filtration stage
    logger.info('Number of pairs analysed from accepted read pool: {:,}'
                .format(pair_parser.pair_filter.analysed()))
    logger.info('Number of pairs accepted: {:,}'.format(n_pairs_accepted))

    if n_pairs_accepted <= 0:
        raise InsufficientDataException('No read-pairs were accepted during parsing.')

    logger.info('Number of pairs trans-mapping: {:,} ({:#.4g}% of pairs)'
                .format(pair_parser.pair_filter.count('trans'),
                        pair_parser.pair_filter.fraction('trans') * 100))
    logger.info('Number of pairs cis-mapping: {:,} ({:#.4g}% of pairs)'
                .format(pair_parser.pair_filter.count('cis'),
                        pair_parser.pair_filter.fraction('cis') * 100))

    logger.info('Number of paired reads that fully align: {:,} ({:#.4g}% of paired)'
                .format(read_counts['full_align'],
                        read_counts['full_align'] / n_paired_reads * 100))
    logger.info('Number of paired reads whose alignment terminates early: {:,} ({:#.4g}% of paired)'
                .format(read_counts['early_term'],
                        read_counts['early_term'] / n_paired_reads * 100))
    logger.info('Number of paired reads not ending in cut-site remnant: {:,} ({:#.4g}% of paired)'
                .format(read_counts['no_site'],
                        read_counts['no_site'] / n_paired_reads * 100))

    logger.info('Number of short-range cis-mapping pairs: {:,} ({:#.4g}% of cis)'
                .format(short_count,
                        short_count / pair_parser.pair_filter.count('cis') * 100))

    # HiCPro style classifications
    logger.info('Number of pairs [dangling end]: {:,}  ({:#.4g}% of cis)'
                .format(class_counts['dangling'],
                        class_counts['dangling'] / n_cis_pairs * 100))
    logger.info('Number of pairs [self-circle]: {:,}  ({:#.4g}% of cis)'
                .format(class_counts['self_circle'],
                        class_counts['self_circle'] / n_cis_pairs * 100))
    logger.info('Number of pairs [invalid FF/RR]: {:,}  ({:#.4g}% of cis)'
                .format(class_counts['ffrr_invalid'],
                        class_counts['ffrr_invalid'] / n_cis_pairs * 100))
    logger.info('Number of pairs [religation]: {:,}  ({:#.4g}% of cis)'
                .format(class_counts['religation'],
                        class_counts['religation'] / n_cis_pairs * 100))
    logger.info('Number of pairs [valid FF/RR]: {:,}  ({:#.4g}% of cis)'
                .format(class_counts['ffrr_valid'],
                        class_counts['ffrr_valid'] / n_cis_pairs * 100))
    logger.info('Number of pairs [valid RF]: {:,}  ({:#.4g}% of cis)'
                .format(class_counts['rf_valid'],
                        class_counts['rf_valid'] / n_cis_pairs * 100))
    logger.info('Number of pairs [valid FR]: {:,}  ({:#.4g}% of cis)'
                .format(class_counts['fr_valid'],
                        class_counts['fr_valid'] / n_cis_pairs * 100))

    # by default, we assume no fraction has gone unobserved
    if short_count == 0:
        emp_mean = None
        logger.warning('Cannot estimate insert length, no short-range cis-mapping pairs [< 1000nt].')
    else:
        emp_mean = short_sum / short_count
        if short_count > 10 * emp_mean:
            logger.info('Observed short-range mean and median of pair separation: {:.0f}nt {:.0f}nt'
                        .format(emp_mean,
                                short_median_est))
        else:
            logger.warning('Too few short-range pairs for reliable streaming median estimation')
            logger.info('Observed mean of short-range pair separation: {:.0f}nt'.format(emp_mean))
    report['obs_insert_mean'] = emp_mean
    report['obs_insert_median'] = short_median_est

    # predict unobserved fraction using a uniform model of junction location
    # across the observed mean insert length
    mean_read_len = cumulative_length / n_paired_reads
    logger.info('Observed mean read length for paired reads: {:.0f}nt'.format(mean_read_len))
    report['mean_readlen'] = mean_read_len

    # the observed mean could not be estimated
    if emp_mean is None:
        # this will be used later as an indicator of absence
        unobs_frac = None
    else:
        unobs_frac = 1 - observed_fraction(is_phase, int(mean_read_len), int(emp_mean))

        if unobs_frac < 0:
            unobs_frac = 0
            logger.warning('For observed insert length of {:.0f}nt, estimated unobserved fraction '
                           'is invalid (<0). Setting to zero.'
                           .format(emp_mean))
        else:
            logger.info('For observed insert length of {:.0f}nt, estimated unobserved fraction: {:#.4g}'
                        .format(emp_mean,
                                unobs_frac))

        report['unobs_frac'] = unobs_frac

    # enzyme statistics
    enz = ligation_info.enzyme_name
    enz_stats = {'cs_start': enzyme_counts['cs_start'],
                 'cs_term': enzyme_counts['cs_term'],
                 'cs_full': enzyme_counts['cs_full'],
                 'partial_readthru': enzyme_counts['partial_readthru'],
                 'read_thru': enzyme_counts['read_thru'],
                 'is_split':  enzyme_counts['is_split'],
                 'upper_bound': enzyme_counts['cs_term'] - enzyme_counts['cs_full'],
                 'elucidation': ligation_info.elucidation,
                 'vestigial': ligation_info.vestigial}
    report['enzyme_stats'] = {enz: enz_stats}

    logger.info('For {}, number of paired reads which began with complete cut-site: {:,} ({:#.4g}%)'
                .format(enz,
                        enzyme_counts['cs_start'],
                        enzyme_counts['cs_start'] / n_paired_reads * 100))

    logger.info('For {}, cut-site {} has remnant {}'
                .format(enz,
                        ligation_info.elucidation,
                        ligation_info.vestigial))
    logger.info('For {}, the expected fraction by random chance at 50% GC: {:#.3f}%'
                .format(enz,
                        1 / 4 ** ligation_info.vest_len * 100))
    logger.info('For {}, number of paired reads whose alignment ends with cut-site remnant: {:,} ({:#.4g}%)'
                .format(enz,
                        enzyme_counts['cs_term'],
                        enzyme_counts['cs_term'] / n_paired_reads * 100))
    logger.info('For {}, number of paired reads that fully aligned and end with cut-site remnant: {:,} ({:#.4g}%)'
                .format(enz,
                        enzyme_counts['cs_full'],
                        enzyme_counts['cs_full'] / n_paired_reads * 100))

    delta_cs = enzyme_counts['cs_term'] - enzyme_counts['cs_full']
    p_ub = delta_cs / n_paired_reads
    logger.info('For {}, upper bound of read-thru events: {:,} ({:#.4g}%)'
                .format(enz,
                        delta_cs,
                        delta_cs / n_paired_reads * 100))

    logger.info('For {}, number of paired reads whose alignment ends with partial read-thru: {:,} ({:#.4g}%)'
                .format(enz,
                        enzyme_counts['partial_readthru'],
                        enzyme_counts['partial_readthru'] / n_paired_reads * 100))

    p_obs = enzyme_counts['read_thru'] / n_paired_reads
    logger.info('For {}, number of paired reads with observable read-thru: {:,} ({:#.4g}%)'
                .format(enz,
                        enzyme_counts['read_thru'],
                        p_obs*100))
    logger.info('For {}, number of paired reads with read-thru and split alignment: {:,} ({:#.4g}%)'
                .format(enz,
                        enzyme_counts['is_split'],
                        enzyme_counts['is_split'] / n_paired_reads * 100))

    enz_stats['fraction'] = {'read_thru': p_obs, 'upper_bound': p_ub}

    logger.info('For {} the raw estimation of Hi-C fraction: ({:#.4g} - {:#.4g}%)'
                .format(enz,
                        p_obs * 100,
                        p_ub * 100))

    enz_stats['adj_fraction'] = {}
    if unobs_frac is None:
        logger.info('For {} there were not enough observations, no adjustment will be made to Hi-C fraction'
                    .format(enz))
        enz_stats['adj_fraction'].update({'read_thru': None, 'upper_bound': None})
    elif unobs_frac > 0:
        enz_stats['adj_fraction'].update({
            'read_thru': p_obs * (1 / (1 - unobs_frac)),
            'upper_bound': p_ub * (1 / (1 - unobs_frac))
        })

        logger.info('For {} the adjusted estimation of Hi-C fraction: ({:#.4g} - {:#.4g}%)'
                    .format(enz,
                            p_obs * (1 / (1 - unobs_frac)) * 100,
                            p_ub * (1 / (1 - unobs_frac)) * 100))

    # long-range bin counts, with formatting to align fields

    field_width = np.maximum(7, np.log10(np.maximum(long_bins, long_counts)).astype(int) + 1)
    logger.info('Long-range intervals: {:>{w[0]}d}nt, {:>{w[1]}d}nt, {:>{w[2]}d}nt'
                .format(*long_bins, w=field_width))
    logger.info('Number of pairs:       {:>{w[0]},d},   {:>{w[1]},d},   {:>{w[2]},d}'
                .format(*long_counts, w=field_width))
    vs_all_cis = long_counts.astype(np.float) / n_cis_pairs * 100
    logger.info('Relative to all cis:  {:>{w[0]},.4g}%,  {:>{w[1]},.4g}%,  {:>{w[2]},.4g}%'
                .format(*vs_all_cis, w=field_width))
    vs_accepted = long_counts.astype(np.float) / n_pairs_accepted * 100
    logger.info('Relative to accepted: {:>{w[0]},.4g}%,  {:>{w[1]},.4g}%,  {:>{w[2]},.4g}%'
                .format(*vs_accepted, w=field_width))

    report['separation_histogram'] = {'counts': long_counts.tolist(),
                                      'vs_all_cis': vs_all_cis.tolist(),
                                      'vs_accepted': vs_accepted.tolist()}

    # append report to file
    if report_path is not None:
        write_jsonline(report_path, report)
