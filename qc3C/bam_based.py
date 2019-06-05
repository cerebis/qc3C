import logging
import multiprocessing
import numpy as np
import os
import pysam
import subprocess
import tqdm
from typing import Optional
from abc import ABC, abstractmethod

from qc3C.exceptions import NameSortingException
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance, LigationInfo
from qc3C.utils import init_random_state, warn_if

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

    opts = ['samtools', 'view', '-c']
    if max_cpu > 1:
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

    count = int(proc.stdout.readline())
    if paired and count % 2 != 0:
        logger.warning('When counting paired reads, the total was not divisible by 2.')
    return count


def pair_separation(r1: pysam.AlignedSegment, r2: pysam.AlignedSegment) -> Optional[int]:
    """
    Determine the separation between R1 and R2, where the distance is measured from the beginning of each
    read's alignment.
    :param r1: read 1
    :param r2: read 2
    :return: separation d in base-pairs between r1 and r2
    """

    assert r1.reference_id == r2.reference_id, 'R1 and R2 are not mapped to the same reference'

    # Treat FR pairs only
    r1rev, r2rev = r1.is_reverse, r2.is_reverse
    if r1rev and not r2rev or r2rev and not r1rev:
        # standardise order
        if r1rev:
            r1, r2 = r2, r1
        # read starts and length
        x1 = r1.pos
        x2 = r2.pos + r2.alen

        d = abs(x2 - x1)

        # only consider situations where separation of start-points is a positive quantity
        if d > 0:
            return d
        # elif d == 0:
        #     print(r1)
        #     print(r2)
        #     print(x1, x2, r1rev, r2rev)
        #     raise NotImplemented()

    return None


class Counter(ABC):

    def __init__(self, counts: dict):
        self.counts = {'all': 0, 'unmapped': 0, 'sample': 0, 'accepted': 0}
        self.counts.update(counts)

    @abstractmethod
    def accept(self, **kwargs: pysam.AlignedSegment) -> bool:
        pass

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
            return self.count(category) / self.counts['accepted']
        elif against == 'all':
            return self.count(category) / self.counts['all']
        elif against == 'analyzed':
            return self.count(category) / self.analyzed()
        else:
            raise RuntimeError('parameter \"against\" must be one of [accepted, analyzed or all]')

    def analyzed(self) -> int:
        """
        The number of reads/pairs analyzed is all items minus those skipped from sub-sampling
        or which were unmapped.
        :return: The number of items analyzed
        """
        return self.counts['all'] - self.counts['sample'] - self.counts['unmapped']

    def rejected(self) -> int:
        """
        :return: The number of reads/pairs which were rejected due to specified constraints.
        """
        return self.analyzed() - self.counts['accepted']


class ReadFilter(Counter):

    def __init__(self, sample_rate: float, min_mapq: int, min_reflen: int, min_match: int,
                 no_secondary: bool, no_supplementary: bool, no_refterm: bool,
                 ref_lengths, random_state: np.random.RandomState = None):
        """
        Filter reads based on a number of criteria, maintaining counts of the various categorisations
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
        Test if a read passes the specified crtiera.
        :param r: the read to test
        :return: True if the read is accepted
        """
        self.counts['all'] += 1

        #
        # immediately stop from some conditions
        #

        if r.is_unmapped:
            self.counts['unmapped'] += 1
            return False

        # for sampling, don't continue to analyse, just return
        if self.sample_rate is not None and self.sample_rate < self.unif():
            self.counts['sample'] += 1
            return False

        #
        # test for each condition, otherwise
        #
        _accept = True
        if self.min_mapq is not None and r.mapping_quality < self.min_mapq:
            self.counts['mapq'] += 1
            _accept = False
        if self.min_reflen is not None and self.ref_lengths[r.reference_id] < self.min_reflen:
            self.counts['ref_len'] += 1
            _accept = False
        if self.no_secondary and r.is_secondary:
            self.counts['secondary'] += 1
            _accept = False
        if self.no_supplementary and r.is_supplementary:
            self.counts['supplementary'] += 1
            _accept = False
        if self.min_match is not None:
            cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
            if cig[0] != 0 or cig[1] < self.min_match:
                self.counts['weak'] += 1
                _accept = False
        if self.no_refterm and r.query_alignment_length < r.query_length and \
                (r.reference_end >= self.ref_lengths[r.reference_id] or r.reference_start == 0):
            self.counts['ref_term'] += 1
            _accept = False

        if _accept:
            self.counts['accepted'] += 1

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
            return self.cis_count()
        else:
            return super().count(category)

    def cis_count(self) -> int:
        """
        :return: the count of cis-mapping pairs
        """
        return self.counts['all'] - self.counts['trans']

    def accept(self, r1: pysam.AlignedSegment, r2: pysam.AlignedSegment) -> bool:
        """
        Test of a pair passes the specified criteria
        :param r1: the first read
        :param r2: the second read
        :return: True if pair is accepted
        """
        _accept = True
        self.counts['all'] += 1
        if r1.reference_id != r2.reference_id:
            # always count trans, even if not emitted
            self.counts['trans'] += 1
            if self.no_trans:
                _accept = False
        if _accept:
            self.counts['accepted'] += 1
        return _accept


class read_pairs(object):

    def __init__(self, bam_path: str, random_state: np.random.RandomState = None, sample_rate: float = None,
                 min_mapq: int = None, min_reflen: int = None, min_match: int = None, no_trans: bool = False,
                 no_secondary: bool = False, no_supplementary: bool = False, no_refterm: bool = False,
                 threads: int = 1, count_reads: bool = False, show_progress: bool = False, is_ipynb: bool = False):
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

        if count_reads:
            logger.info('Counting alignments in {}'.format(bam_path))
            self.n_reads = count_bam_reads(bam_path, max_cpu=threads)
            logger.info('Found {:,} alignments to analyse'.format(self.n_reads))
        else:
            self.n_reads = None

        self.bam = pysam.AlignmentFile(bam_path, mode='rb', threads=threads)
        _header = self.bam.header['HD']
        if 'SO' not in _header or _header['SO'] != 'queryname':
            raise NameSortingException(bam_path)

        self.random_state = random_state
        self.show_progress = show_progress
        self.progress = None
        self.is_ipynb = is_ipynb

        # applying sample_rate to reads means that the pass rate of pairs has an upper bound of sample_rate ** 2.
        # therefore take the square, so that the number of pairs is closer to what users expect.
        # this is particularly noticeable for low sample rates (ie passing only 10% of reads means only 1% of pairs)
        if sample_rate is not None:
            pair_sample_rate = sample_rate ** 0.5
        else:
            pair_sample_rate = None

        self.read_filter = ReadFilter(pair_sample_rate, min_mapq, min_reflen, min_match, no_secondary, no_supplementary,
                                      no_refterm, np.array([li for li in self.bam.lengths]), random_state)
        self.pair_filter = PairFilter(no_trans=no_trans)

        self.bam_iter = self.bam.fetch(until_eof=True)
        if not self.show_progress:
            self.progress = None
        elif self.is_ipynb:
            self.progress = tqdm.tqdm_notebook(total=self.n_reads)
        else:
            self.progress = tqdm.tqdm(total=self.n_reads)

    def __enter__(self):
        return self

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


def analyze(bam_file: str, enzymes: list, mean_insert: int, seed: int = None,
            sample_rate: float = None, min_mapq: int = 60, threads: int = 1) -> None:
    """
    Analyze a bam file which contains Hi-C read-pairs mapped to a set of reference sequences.
    This method attempts to assess details which indicate the overall strength of the
    Hi-C signal. Here signal is the proportion of proximity ligation pairs relative to
    shotgun pairs within the mapped read-sets.

    :param bam_file: the path to the bam file
    :param enzymes: the list of enzyme names used in Hi-C library creation
    :param mean_insert: the expected mean insert length of the library
    :param seed: a random integer seed
    :param sample_rate: consider only a portion of all pairs [0..1]
    :param min_mapq: the minimum acceptable mapping quality
    :param threads: the number of threads used in accessing the bam file
    """

    ligation_variants = {}
    for enz in enzymes:
        ligation_variants[enz] = ligation_junction_seq(get_enzyme_instance(enz))

    # for efficiency, disable sampling if rate is 1
    if sample_rate == 1 or sample_rate is None:
        logger.info('Accepting all usable reads')
        sample_rate = None
    else:
        logger.info('Acceptance threshold: {:#.3g}'.format(sample_rate))

    random_state = init_random_state(seed)

    with read_pairs(bam_file, random_state=random_state, sample_rate=sample_rate,
                    min_mapq=min_mapq, min_match=1, min_reflen=500,
                    no_trans=False, no_secondary=True, no_supplementary=True, no_refterm=True,
                    threads=threads, count_reads=True, show_progress=True) as pair_parser:

        cumulative_length = 0
        short_sum = 0
        short_count = 0
        long_counts = np.zeros(3, dtype=np.int)
        long_bins = (1000, 5000, 10000)
        pair_counts = {'total': 0, 'full_align': 0, 'early_term': 0, 'no_site': 0}
        enzyme_counts = {enz: {'cs_full': 0, 'cs_term': 0, 'read_thru': 0,
                               'is_split': 0, 'partial_readthru': 0} for enz in enzymes}

        # begin with user supplied median estimate
        short_median_est = mean_insert

        logger.info('Beginning analysis...')

        for r1, r2 in pair_parser:

            # for cis-pairs, determine their separation
            if r1.reference_id == r2.reference_id:

                d = pair_separation(r1, r2)

                # cis-mapping, well mapped pairs
                if d is not None:

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

                    # keep a count of those pairs which are within the "WGS" region
                    elif 50 < d < 1000:
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

                pair_counts['total'] += 1
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

                    pair_counts['full_align'] += 1

                    # check for cut-site
                    _no_site = True
                    for enz, lig_info in ligation_variants.items():
                        if seq.endswith(lig_info.vestigial):
                            _no_site = False
                            enzyme_counts[enz]['cs_full'] += 1
                            enzyme_counts[enz]['cs_term'] += 1
                            break

                    if _no_site:
                        pair_counts['no_site'] += 1

                # for reads terminating early, we can check for the junction duplication
                else:

                    pair_counts['early_term'] += 1

                    # check for cut-site and potentially the junction
                    _no_site = True
                    for enz, lig_info in ligation_variants.items():

                        max_match = lig_info.junc_len - lig_info.vest_len

                        if aln_seq.endswith(lig_info.vestigial):
                            _no_site = False

                            enzyme_counts[enz]['cs_term'] += 1

                            # look for partial junctions, depending on how much
                            # flanking sequence a read has available
                            if ri.is_reverse:
                                spare = ri.query_alignment_start
                            else:
                                spare = ri.query_length - ri.query_alignment_end

                            m = junction_match_length(seq, ri, lig_info)
                            if (m < max_match and m == spare) or m == max_match:
                                enzyme_counts[enz]['partial_readthru'] += 1

                            # for full matches, keep a separate tally
                            if m == max_match:
                                enzyme_counts[enz]['read_thru'] += 1
                                if ri.has_tag('SA'):
                                    enzyme_counts[enz]['is_split'] += 1

                            break

                    if _no_site:
                        pair_counts['no_site'] += 1

        """
        Report the results
        """

        logger.info('Number of parsed reads: {:,}'
                    .format(pair_parser.read_filter.counts['all']))
        logger.info('Number of analysed reads: {:,} ({:#.2f}% of all)'
                    .format(pair_parser.read_filter.analyzed(),
                            pair_parser.read_filter.fraction('analyzed', 'all')*100))
        logger.info('Number of reads filtered [unmapped]: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.counts['unmapped'],
                            pair_parser.read_filter.fraction('unmapped')*100))
        logger.info('Number of reads filtered [low mapq]: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.counts['mapq'],
                            pair_parser.read_filter.fraction('mapq')*100))
        logger.info('Number of reads filtered [ref length]: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.counts['ref_len'],
                            pair_parser.read_filter.fraction('ref_len')*100))
        logger.info('Number of reads filtered [secondary]: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.counts['secondary'],
                            pair_parser.read_filter.fraction('secondary')*100))
        logger.info('Number of reads filtered [supplementary]: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.counts['supplementary'],
                            pair_parser.read_filter.fraction('supplementary')*100))
        logger.info('Number of reads filtered [weak mapping]: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.counts['weak'],
                            pair_parser.read_filter.fraction('weak')*100))
        logger.info('Number of reads filtered [ref terminated]: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.counts['ref_term'],
                            pair_parser.read_filter.fraction('ref_term')*100))

        # accepted reads
        logger.info('Number of accepted reads: {:,} ({:#.2f}% of analyzed)'
                    .format(pair_parser.read_filter.count('accepted'),
                            pair_parser.read_filter.fraction('accepted')*100))

        # for pairs which have accepted read filtration stage
        logger.info('Number of pairs resulting from accepted read pool: {:,}'
                    .format(pair_parser.pair_filter.analyzed()))
        logger.info('Number of pairs trans-mapping: {:,} ({:#.2f}% of pairs)'
                    .format(pair_parser.pair_filter.counts['trans'],
                            pair_parser.pair_filter.fraction('trans')*100))
        logger.info('Number of pairs cis-mapping: {:,} ({:.2f}% of pairs)'
                    .format(pair_parser.pair_filter.cis_count(),
                            pair_parser.pair_filter.fraction('cis')*100))

        logger.info('Number of paired reads that fully align: {:,} ({:#.2f}% of paired)'
                    .format(pair_counts['full_align'],
                            pair_counts['full_align'] / pair_counts['total']*100))
        logger.info('Number of paired reads whose alignment terminates early: {:,} ({:#.2f}% of paired)'
                    .format(pair_counts['early_term'],
                            pair_counts['early_term'] / pair_counts['total']*100))
        # TODO this can probably be dropped.
        logger.info('Number of paired reads not ending in cut-site remnant: {:,} ({:#.2f}% of paired)'
                    .format(pair_counts['no_site'],
                            pair_counts['no_site'] / pair_counts['total']*100))

        logger.info('Number of short-range cis-mapping pairs: {:,} ({:#.2f}% of cis)'
                    .format(short_count,
                            short_count / pair_parser.pair_filter.cis_count()*100))

        # by default, we assume no fraction has gone unobserved
        if short_count == 0:
            emp_mean = None
            logger.warning('Cannot estimate insert length, no short-range cis-mapping pairs [< 1000nt].')
        else:
            emp_mean = short_sum / short_count
            mean_err = (emp_mean - mean_insert) / mean_insert
            logger.log(warn_if(mean_err > 0.1),
                       'Observed short-range mean pair separation differs from supplied insert length by {:.1f}%'
                       .format(mean_err * 100))

            if short_count > 10 * mean_insert:
                logger.info('Observed short-range mean and median of pair separation: {:.0f}nt {:.0f}nt'
                            .format(emp_mean,
                                    short_median_est))
            else:
                logger.warning('Too few short-range pairs for reliable streaming median estimation')
                logger.info('Observed mean of short-range pair separation: {:.0f}nt'.format(emp_mean))

        # predict unobserved fraction using a uniform model of junction location
        # across the observed mean insert length
        mean_read_len = cumulative_length / pair_counts['total']
        logger.info('Observed mean read length for paired reads: {:.0f}nt'.format(mean_read_len))
        unobs_frac = {}
        for _tag, _insert_len in {'supplied': mean_insert, 'observed': emp_mean}.items():

            # no mean was supplied or the observed could not be estimated
            if _insert_len is None:
                # this will be used later as an indicator of absence
                unobs_frac[_tag] = None
                continue

            unobs_frac[_tag] = (_insert_len - mean_read_len * 2) / _insert_len
            if unobs_frac[_tag] < 0:
                unobs_frac[_tag] = 0
                logger.warning('For {} insert length of {:.0f}nt, estimation of the unobserved fraction '
                               'is invalid (<0). Assuming an unobserved fraction: {:#.4g}'
                               .format(_tag,
                                       _insert_len,
                                       unobs_frac[_tag]))
            else:
                logger.info('For {} insert length of {:.0f}nt, estimated unobserved fraction: {:#.4g}'
                            .format(_tag,
                                    _insert_len,
                                    unobs_frac[_tag]))

        # per-enzyme statistics
        for enz, counts in enzyme_counts.items():
            logger.info('For {}, the expected fraction by random chance at 50% GC: {:#.3f}%'
                        .format(enz,
                                1 / 4 ** ligation_variants[enz].vest_len * 100))
            logger.info('For {}, cut-site {} has remnant {}'
                        .format(enz,
                                ligation_variants[enz].elucidation,
                                ligation_variants[enz].vestigial))
            logger.info('For {}, number of paired reads whose alignment ends with cut-site remnant: {:,} ({:#.2f}%)'
                        .format(enz,
                                counts['cs_term'],
                                counts['cs_term'] / pair_counts['total']*100))

            logger.info('For {}, number of paired reads that fully aligned and end with cut-site remnant: {:,} ({:#.2f}%)'
                        .format(enz,
                                counts['cs_full'],
                                counts['cs_full'] / pair_counts['total']*100))

            delta_cs = counts['cs_term'] - counts['cs_full']
            logger.info('For {}, upper bound of read-thru events: {:,} ({:#.2f}%)'
                        .format(enz,
                                delta_cs,
                                delta_cs / pair_counts['total']*100))

            logger.info('For {}, number of paired reads whose alignment ends with partial read-thru: {:,} ({:#.2f}%)'
                        .format(enz,
                                counts['partial_readthru'],
                                counts['partial_readthru'] / pair_counts['total']*100))

            p_obs = counts['read_thru'] / pair_counts['total']
            logger.info('For {}, number of paired reads with observable read-thru: {:,} ({:#.2f}%)'
                        .format(enz,
                                counts['read_thru'],
                                p_obs*100))
            logger.info('For {}, number of paired reads with read-thru and split alignment: {:,} ({:#.2f}%)'
                        .format(enz,
                                counts['is_split'],
                                counts['is_split'] / pair_counts['total']*100))

            for _tag, _frac in unobs_frac.items():
                if _frac is None:
                    logger.info('For {} there was insufficient {} data,     no adjustment could be made to Hi-C fraction'
                                .format(enz, _tag))
                    continue

                p_ub = delta_cs / pair_counts['total']
                logger.info('For {} based on {} data, adjusted estimation of Hi-C fraction: ({:#.2f}-{:#.2f}%)'
                            .format(enz, _tag,
                                    (p_obs + p_obs * _frac) * 100,
                                    (p_ub + p_ub * _frac) * 100))

        # long-range bin counts, with formatting to align fields
        field_width = np.maximum(7, np.log10(np.maximum(long_bins, long_counts)).astype(int) + 1)
        logger.info('Long-range distance intervals:  {:>{w[0]}d}nt, {:>{w[1]}d}nt, {:>{w[2]}d}nt'
                    .format(*long_bins,
                            w=field_width))
        logger.info('Number of cis-mapping pairs:    {:>{w[0]},d},   {:>{w[1]},d},   {:>{w[2]},d}'
                    .format(*long_counts,
                            w=field_width))
        frac = long_counts.astype(np.float) / pair_parser.pair_filter.cis_count()
        logger.info('Relative fraction of all cis:  {:#.4g},   {:#.4g},   {:#.4g}'
                    .format(*frac,
                            w=field_width))
