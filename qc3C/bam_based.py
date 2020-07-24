import logging
import multiprocessing
import numpy as np
import os
import pysam
import subprocess
import tqdm

from collections import OrderedDict
from abc import ABC, abstractmethod
from astropy.stats import sigma_clipped_stats
from typing import Optional, Tuple, List
from qc3C.exceptions import NameSortingException, UnknownLibraryKitException, MaxObsLimit, InsufficientDataException
from qc3C.ligation import CutSitesDB, Digest
from qc3C.utils import init_random_state, write_jsonline, observed_fraction, write_html_report
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

        # only test for the following conditions if the read alignment has passed earlier tests

        if _accept and self.min_match is not None:
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            if cig[0] != 0 or cig[1] < self.min_match:
                self._counts['weak'] += 1
                _accept = False
        if _accept and self.no_refterm and r.query_alignment_length < r.query_length and \
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

                # only remember the previous accepted alignment
                if pass_b:
                    read_a = read_b
                    pass_a = pass_b

        except StopIteration:
            self.close()


def analyse(enzyme_names: List[str], bam_file: str, fasta_file: str,
            seed: int = None, sample_rate: float = None, min_mapq: int = 60, max_obs: int = None,
            threads: int = 1, report_path: str = None, no_json: bool = False, no_html: bool = False) -> None:
    """
    analyse a bam file which contains Hi-C read-pairs mapped to a set of reference sequences.
    This method attempts to assess details which indicate the overall strength of the
    Hi-C signal. Here signal is the proportion of proximity ligation pairs relative to
    shotgun pairs within the mapped read-sets.

    :param enzyme_names: the enzymes used during digestion (max 2)
    :param bam_file: the path to the bam file
    :param fasta_file: reference sequences
    :param seed: a random integer seed
    :param sample_rate: consider only a portion of all pairs [0..1]
    :param min_mapq: the minimum acceptable mapping quality
    :param max_obs: the maximum number of accepted pairs to inspect
    :param threads: the number of threads used in accessing the bam file
    :param report_path: append a report in single-line JSON format to the given path.
    :param no_json: disable json report
    :param no_html: disable html report
    """

    assert 0 < len(enzyme_names) <= 2, 'only 1 or 2 enzymes can be specified'
    digest = Digest(*enzyme_names, no_ambig=False)
    cutsite_db = CutSitesDB(fasta_file, enzyme_names, use_cache=True, use_tqdm=True)

    # for efficiency, disable sampling if rate is 1
    if sample_rate == 1 or sample_rate is None:
        logger.info('Accepting all usable reads')
        sample_rate = None
    else:
        logger.info('Acceptance threshold: {:#.3g}'.format(sample_rate))

    random_state = init_random_state(seed)

    logger.info('Beginning analysis...')

    mp_progress = None
    if max_obs is not None:
        show_progress = count_reads = False
        # special progress bar for max-obs runs
        mp_progress = tqdm.tqdm(total=max_obs, desc='Pairs')
    else:
        show_progress = count_reads = True

    cumulative_length = 0
    short_count = 0
    read_counts = {'full_align': 0, 'early_term': 0, 'no_site': 0}
    class_counts = {'dangling': 0, 'self_circle': 0, 'religation': 0, 'ffrr_invalid': 0,
                    'fr_valid': 0, 'rf_valid': 0, 'ffrr_valid': 0}
    digest_counts = {'cs_full': 0, 'cs_term': 0, 'cs_start': 0, 'read_thru': 0,
                     'is_split': 0}

    n_accepted_obs = 0

    pair_parser = read_pairs(bam_file, random_state=random_state, sample_rate=sample_rate,
                             min_mapq=min_mapq, min_match=1, min_reflen=500,
                             no_trans=False, no_secondary=True, no_supplementary=True, no_refterm=True,
                             threads=threads, count_reads=count_reads, show_progress=show_progress)

    startswith_cutsite = digest.cutsite_searcher('startswith')
    startswith_junction = digest.junction_searcher('startswith')
    endswith_vestigial = digest.vestigial_searcher('endswith')

    junction_tracker = digest.tracker('junction')
    vestigial_tracker = digest.tracker('vestigial')

    pair_separations = []

    try:

        for r1, r2 in pair_parser:

            # cis-mapping pairs, lets categorise using the standard nomenclature
            cis_pair = r1.reference_id == r2.reference_id

            num_sites = 0

            if cis_pair:

                d, num_sites, flipped = cutsite_db[r1.reference_name].pair_info(r1, r2)

                if flipped:
                    r1, r2 = r2, r1

                # invalid pair
                if num_sites == 0:
                    if not r1.is_reverse and r2.is_reverse:
                        class_counts['dangling'] += 1
                    elif r1.is_reverse and not r2.is_reverse:
                        class_counts['self_circle'] += 1
                    else:
                        class_counts['ffrr_invalid'] += 1

                # when separation has been determined or estimated
                if d is not None and d > 0:

                    # track the number of pairs observed for the requested separation distances
                    if d > 50:
                        pair_separations.append(d)

                    if d < 1000:
                        short_count += 1

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
                    if startswith_cutsite(seq) is not None:
                        digest_counts['cs_start'] += 1

                    vest_match = endswith_vestigial(seq)
                    if vest_match is not None:
                        vestigial_tracker[vest_match.group()] += 1
                        _no_site = False
                        digest_counts['cs_full'] += 1
                        digest_counts['cs_term'] += 1

                    if _no_site:
                        read_counts['no_site'] += 1

                # for reads terminating early, we can check for the junction duplication
                else:

                    read_counts['early_term'] += 1

                    # check for cut-site and potentially the junction
                    _no_site = True

                    if startswith_cutsite(seq):
                        digest_counts['cs_start'] += 1

                    vest_match = endswith_vestigial(aln_seq)
                    if vest_match is not None:
                        vestigial_tracker[vest_match.group()] += 1
                        _no_site = False
                        digest_counts['cs_term'] += 1
                        vstart = vest_match.start()
                        # index of vestigial start within full read sequence
                        if ri.is_reverse:
                            i = vstart + ri.query_length - ri.query_alignment_end
                        else:
                            i = vstart + ri.query_alignment_start
                        # now, does the read continue and contain the full junction
                        junc_match = startswith_junction(seq[i:])
                        # for full matches, keep a separate tally
                        if junc_match is not None:
                            junction_tracker[junc_match.group()] += 1
                            digest_counts['read_thru'] += 1
                            if ri.has_tag('SA'):
                                digest_counts['is_split'] += 1

                    if _no_site:
                        read_counts['no_site'] += 1

            if cis_pair and num_sites >= 1:
                if not r1.is_reverse and r2.is_reverse:
                    if num_sites == 1:
                        class_counts['religation'] += 1
                    class_counts['fr_valid'] += 1
                elif r1.is_reverse and not r2.is_reverse:
                    class_counts['rf_valid'] += 1
                else:
                    class_counts['ffrr_valid'] += 1

            if max_obs is not None:
                n_accepted_obs += 1
                mp_progress.update()
                if n_accepted_obs >= max_obs:
                    raise MaxObsLimit

    except MaxObsLimit:
        if mp_progress is not None:
            mp_progress.close()
            logger.info('Reached user-defined observation limit [{}]'.format(max_obs))

    # some shorthand variables which will be used repeatedly
    n_pairs_accepted = pair_parser.pair_filter.count('accepted')
    n_paired_reads = 2 * n_pairs_accepted
    n_cis_pairs = pair_parser.pair_filter.count('cis')
    pair_separations = np.array(pair_separations, dtype=np.int32)

    #
    # Initial values for report
    #

    report = OrderedDict({
        'mode': 'bam',
        'runtime_info': runtime_info(),
        'input_args': {'bam_file': bam_file,
                       'fasta_file': fasta_file,
                       'enzymes': enzyme_names,
                       'seed': seed,
                       'sample_rate': sample_rate,
                       'min_mapq': min_mapq,
                       'max_obs': max_obs},
        'n_parsed_reads': pair_parser.read_filter.count('all'),
        'n_analysed_reads': pair_parser.read_filter.analysed(),
        'n_accepted_reads': pair_parser.read_filter.count('accepted'),
        'n_skipped_reads': pair_parser.read_filter.count('sample'),
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
    })

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

    uninf_pairs = class_counts['dangling'] + class_counts['self_circle'] + \
                  class_counts['ffrr_invalid'] + class_counts['religation']

    logger.info('Number of uninformative pairs: {:,}  ({:#.4g}% of cis)'
                .format(uninf_pairs, uninf_pairs / n_cis_pairs * 100))

    inf_pairs = class_counts['ffrr_valid'] + class_counts['rf_valid'] + class_counts['fr_valid']
    logger.info('Number of informative pairs: {:,}  ({:#.4g}% of cis)'
                .format(inf_pairs, inf_pairs / n_cis_pairs * 100))

    lt_1kb = pair_separations < 1000
    if lt_1kb.sum() < 20:
        emp_mean = None
        emp_median = None
        logger.warning('Insufficient data to estimate insert size due to '
                       'a lack of short-range cis-mapping pairs [< 1000nt].')
    else:
        # use a robust estimator for the mean and median.
        emp_mean, emp_median, emp_sd = sigma_clipped_stats(pair_separations[lt_1kb], sigma=2, maxiters=100)
        logger.info('Estimated insert size mean and median: {:.0f}nt {:.0f}nt'
                    .format(emp_mean,
                            emp_median))
    report['obs_insert_mean'] = emp_mean
    report['obs_insert_median'] = emp_median

    # predict unobserved fraction using a uniform model of junction location
    # across the observed mean insert length
    mean_read_len = cumulative_length / n_paired_reads
    logger.info('Observed mean read length for paired reads: {:.0f}nt'.format(mean_read_len))
    report['mean_readlen'] = mean_read_len

    # the observed mean could not be estimated
    if emp_median is None:
        logger.warning('Unobserved fraction not estimated as insert size was not available')
    else:
        obs_frac = observed_fraction(round(mean_read_len), round(emp_median), 'binary',
                                     junc_size=digest.longest_junction())
        logger.info('For observed insert size of {:.0f}nt, estimated unobserved fraction: {:#.4g}'
                    .format(emp_median, 1 - obs_frac))
        report['unobs_fraction'] = 1 - obs_frac

    # digest statistics
    digest_stats = {'cs_start': digest_counts['cs_start'],
                    'cs_term': digest_counts['cs_term'],
                    'cs_full': digest_counts['cs_full'],
                    'read_thru': digest_counts['read_thru'],
                    'is_split':  digest_counts['is_split'],
                    'vestigial': digest.unique_vestigial()}
    report['digest_stats'] = digest_stats

    for ci in digest.cutsites.values():
        logger.info('The enzyme {} has cut-site {}'.format(ci.name, ci.site))
    for ji in digest.junctions.values():
        logger.info('The enzymatic combination {}/{} produces the {}nt ligation junction {}'
                    .format(ji.enz5p, ji.enz3p, ji.junc_len, ji.junction))

    logger.info('The number of paired reads which began with complete cut-site: {:,} ({:#.4g}%)'
                .format(digest_counts['cs_start'],
                        digest_counts['cs_start'] / n_paired_reads * 100))

    logger.info('The digest contains the following possible remnants {}'
                .format(digest_stats['vestigial']))
    logger.info('The number of paired reads whose alignment ends with cut-site remnant: {:,} ({:#.4g}%)'
                .format(digest_counts['cs_term'],
                        digest_counts['cs_term'] / n_paired_reads * 100))

    p_obs = digest_counts['read_thru'] / n_paired_reads
    logger.info('The number of paired reads with observable read-thru: {:,} ({:#.4g}%)'
                .format(digest_counts['read_thru'],
                        p_obs*100))
    logger.info('The number of paired reads with read-thru and split alignment: {:,} ({:#.4g}%)'
                .format(digest_counts['is_split'],
                        digest_counts['is_split'] / n_paired_reads * 100))

    # Tally coarse bins used by convention in Hi-C field
    long_bins = np.array([1000, 5000, 10000], dtype=np.int)
    long_counts = np.array([(pair_separations >= long_bins[0]).sum(),
                            (pair_separations >= long_bins[1]).sum(),
                            (pair_separations >= long_bins[2]).sum()], dtype=np.int)

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

    report['separation_bins'] = {'bins': long_bins,
                                 'counts': long_counts,
                                 'vs_all_cis': vs_all_cis,
                                 'vs_accepted': vs_accepted}

    report['digestion'] = {'cutsites': digest.to_dict('cutsite'),
                           'junctions': digest.to_dict('junction')}

    report['junction_frequency'] = {}
    for _e, _counts in digest.gather_tracker(junction_tracker).items():
        for _j, _n in _counts.items():
            logger.info(f'For {_e} junction sequence {_j} found: {_n}')
            report['junction_frequency'][f'{_e} {_j}'] = _n

    report['remnant_frequency'] = {}
    for _e, _counts in digest.gather_tracker(vestigial_tracker).items():
        for _v, _n in _counts.items():
            logger.info(f'For {_e} cutsite remnant {_v} found: {_n}')
            report['remnant_frequency'][f'{_e} {_v}'] = _n

    # write all pair separations to a file
    # np.savetxt('pair_separations.txt', np.array(pair_separations, dtype=np.int), fmt='%d')
    hist, edges = np.histogram(pair_separations,
                               bins=np.logspace(np.log10(10), np.log10(1e6), 100), density=True)
    report['separation_histogram'] = {'mid_points': 0.5*(edges[1:] + edges[:-1]), 'counts': hist}

    if not no_json:
        write_jsonline(report_path, report, suffix='json', append=False)

    if not no_html:
        # just drop this large record prior to creating the html version
        del report['separation_histogram']
        write_html_report(report_path, report)
