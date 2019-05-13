import logging
import multiprocessing
import numpy as np
import os
import pandas
import pysam
import re
import subprocess
import tqdm

from collections import namedtuple
from collections.abc import Collection
from recordclass import recordclass

from qc3C.exceptions import NameSortingException
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance
from qc3C.utils import init_random_state

logger = logging.getLogger(__name__)

# immutable type used in pair logging
PairInfo = namedtuple('pair_info', ('name', 'pos', 'length', 'is_reverse', 'cigarstring'))

# mutable tuples used in storing counting statistics on evidence of proximity ligation events.
# information pertaining to the entire read-set
GlobalInfo = recordclass('global_info',
                         ('ref_term', 'no_site', 'full_align', 'early_term', 'low_mapq', 'total'),
                         defaults=(0,) * 6)

# information specific an enzyme's action
CutSiteInfo = recordclass('cutsite_info',
                          ('cs_term', 'cs_full', 'read_thru', 'is_split'),
                          defaults=(0,) * 4)


class QcInfo(object):
    """
    Represents the evidence collected for the action of Hi-C proximity ligation
    over an entire BAM file.
    """

    def __init__(self, enzymes):
        """
        :param enzymes: the enzymes used in creating the Hi-C library
        """
        if not isinstance(enzymes, Collection) or isinstance(enzymes, str):
            enzymes = [enzymes]
        self.global_info = GlobalInfo()
        self.enzyme = {en: CutSiteInfo() for en in enzymes}

    def total(self):
        """
        :return: the total number of reads
        """
        return self.global_info.total

    def accepted(self):
        """
        :return: the total number of mapped reads
        """
        return self.global_info.ref_term + self.global_info.full_align + self.global_info.early_term


# Mapping of cigar characters to code values
CODE2CIGAR = dict((y, x) for x, y in enumerate("MIDNSHP=X"))

# Pattern that finds each unit of a cigar i.e. 10M or 9H
CIGAR_ANY = re.compile(r"(\d+)([MIDNSHP=X])")


def cigar_to_tuple(cigar):
    """
    Convert a CIGAR string into code values
    :param cigar:
    :return:
    """
    return [(CODE2CIGAR[t[1]], int(t[0])) for t in CIGAR_ANY.findall(cigar)]


def parse_secondary_alignment_tag(r):
    """
    Extract the secondary alignment tag (SA) information.
    :param r: the read to inspect
    :return: dictionary of fields extracted from the SA tag, or None if no SA tag exists
    """
    if not r.has_tag('SA'):
        return None

    _tag = r.get_tag('SA')
    for aln_i in _tag.split(';'):
        if not aln_i:
            continue
        ti = aln_i.split(',')
        pos = int(ti[1]) - 1
        is_rev = {'-': True, '+': False}[ti[2]]
        cigtup = cigar_to_tuple(ti[3])
        alen = sum([num for op, num in cigtup if op == 0])
        tot = sum([num for op, num in cigtup])

        return {'ref': ti[0],
                'pos': pos,
                'is_reverse': is_rev,
                'cigar': ti[3],
                'cigartuple': cigtup,
                'mapq': int(ti[4]),
                'nm': int(ti[5]),
                'alen': alen,
                'total': tot}


COMPLEMENT_TABLE = str.maketrans('acgtumrwsykvhdbnACGTUMRWSYKVHDBN',
                                 'TGCAAnnnnnnnnnnnTGCAANNNNNNNNNNN')


def revcomp(seq):
    """
    Reverse complement a string representation of a sequence. This uses string.translate.
    :param seq: input sequence as a string
    :return: revcomp sequence as a string
    """
    return seq.translate(COMPLEMENT_TABLE)[::-1]


def get_forward_strand(r):
    """
    Return the forward/top strand orientation of a read's sequence
    :param r: the read
    :returns: forward (5p-3p) sequence of the read
    """
    seq = r.seq
    if r.is_reverse:
        seq = revcomp(r.seq)
    return seq


def exe_exists(exe_name):
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


def count_bam_reads(file_name: str, paired: bool = False, mapped: bool = False, mapq: int = 1, max_cpu: int = None):
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

    if mapq:
        opts.append('-q{}'.format(mapq))

    proc = subprocess.Popen(opts + [file_name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    count = int(proc.stdout.readline())
    if paired and count % 2 != 0:
        logger.warning('When counting paired reads, the total was not divisible by 2.')
    return count


def print_report(report, sep='\t'):
    """
    Create a tabular report in CSV format
    :param report: the QcInfo instance to report
    :param sep: separator to use in table creation
    """

    total = report.total()
    _cnames = ['enzyme', 'variable', 'read_count', 'vs_total', 'vs_accepted']

    df = pandas.DataFrame([['all', 'total', total, None, None],
                           ['all', 'accepted', report.accepted(), None, None],
                           ['all', 'low_mapq', report.global_info.low_mapq, None, None],
                           ['all', 'fully aligned', report.global_info.full_align, None, None],
                           ['all', 'early termination', report.global_info.early_term, None, None],
                           ['all', 'ref termination', report.global_info.ref_term, None, None],
                           ['all', 'no cutsite', report.global_info.no_site, None, None]],
                          columns=_cnames)

    for enz, inf in report.enzyme.items():
        df_enz = pandas.DataFrame([[enz, 'cs fully aligned', inf.cs_full, None, None],
                                   [enz, 'cs termination', inf.cs_term, None, None],
                                   [enz, 'read-thru', inf.read_thru, None, None],
                                   [enz, 'split alignment', inf.is_split, None, None]],
                                  columns=_cnames)
        df = df.append(df_enz, ignore_index=True)

    df.loc[1:, 'vs_total'] = df.loc[1:, 'read_count'] / df.loc[0, 'read_count']
    df.loc[2:, 'vs_accepted'] = df.loc[2:, 'read_count'] / df.loc[1, 'read_count']
    df[['vs_total', 'vs_accepted']] = df[['vs_total', 'vs_accepted']].apply(pandas.to_numeric)

    print('\nQC details:')
    print(df.to_csv(None, sep=sep, float_format="%#.4g", index=False))


def pair_separation(r1: pysam.AlignedSegment, r2: pysam.AlignedSegment):
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

        d = x2 - x1

        # only consider situations where separation of start-points is a positive quantity
        if d > 0:
            return d
        # elif d == 0:
        #     print(r1)
        #     print(r2)
        #     print(x1, x2, r1rev, r2rev)
        #     raise NotImplemented()

    return None


def strong_match(r, minmatch):
    cig = r.cigartuples[-1] if r.is_reverse else r.cigartuples[0]
    return not r.is_secondary and not r.is_supplementary and cig[0] == 0 and cig[1] >= minmatch


def analyze(bam_file, enzymes, mean_insert, seed=None, sample_rate=None, min_mapq=60, threads=1, sep='\t'):

    report = QcInfo(enzymes)

    ligation_variants = []
    for ename in enzymes:
        ligation_variants.append(ligation_junction_seq(get_enzyme_instance(ename)))

    # for efficiency, disable sampling if rate is 1
    if sample_rate == 1 or sample_rate is None:
        logger.info('Accepting all usable reads')
        sample_rate = None
    else:
        logger.info('Acceptance threshold: {:#.3g}'.format(sample_rate))

    with pysam.AlignmentFile(bam_file, 'rb', threads=threads) as bam:

        if 'SO' not in bam.header['HD'] or bam.header['HD']['SO'] != 'queryname':
            raise NameSortingException(bam_file)

        ref_lengths = [li for li in bam.lengths]

        r_prev = None
        all_cis_count = 0
        short_sum = 0
        short_count = 0
        long_counts = np.zeros(3, dtype=np.int)
        long_bins = (1000, 5000, 10000)

        # begin with user supplied median estimate
        short_median_est = mean_insert

        logger.info('Counting alignments...')
        n_reads = count_bam_reads(bam_file, max_cpu=threads)
        logger.info('There were {} alignments in {}'.format(n_reads, bam_file))

        random_state = init_random_state(seed)
        unif = random_state.uniform

        logger.info('Beginning analysis...')
        with tqdm.tqdm(desc='Progress', total=n_reads) as progress:

            bam_iter = bam.fetch(until_eof=True)

            while True:

                try:
                    r = next(bam_iter)
                    progress.update()
                except StopIteration:
                    break

                # if specified, collect only a sampling
                if sample_rate is not None and sample_rate < unif():
                    continue

                report.global_info.total += 1

                if r.mapping_quality < min_mapq:
                    report.global_info.low_mapq += 1
                    continue

                # mate tracking
                if r_prev is not None:

                    # cis-mapping, well mapped pairs
                    if r.query_name == r_prev.query_name and \
                            r.reference_id == r_prev.reference_id and \
                            strong_match(r, 10) and \
                            strong_match(r_prev, 10):

                        d = pair_separation(r, r_prev)
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
                            if 50 < d < 1000:
                                short_sum += d
                                short_count += 1

                                # frugal median estimator
                                if short_median_est > d:
                                    short_median_est -= 1
                                elif short_median_est < d:
                                    short_median_est += 1

                            all_cis_count += 1

                r_prev = r

                if r.reference_end >= ref_lengths[r.reference_id] or r.reference_start == 0:
                    # reads which align to the ends of references are ignored
                    report.global_info.ref_term += 1
                    continue

                if r.query_length == r.reference_length:
                    report.global_info.full_align += 1
                    # fully aligned reads can't be tested for the junction but
                    # we can still test for the cut-site
                    seq = get_forward_strand(r)
                    for lig in ligation_variants:
                        if seq.endswith(lig.end_match):
                            report.enzyme[lig.enzyme_name].cs_full += 1
                            break
                    continue

                report.global_info.early_term += 1

                # inspect all sequences in as 5'-3'
                seq = get_forward_strand(r)

                # the aligned sequence, which should end with cut-site
                aln_seq = seq[:r.reference_length]

                # check that a cut-site exists on the end
                found_lig = False
                for lig in ligation_variants:
                    if aln_seq.endswith(lig.end_match):
                        found_lig = True
                        report.enzyme[lig.enzyme_name].cs_term += 1

                        # a proximity ligation product should contain a characteristic
                        # sequence which duplicates a portion of the cut-site. Check
                        # and see if the read contains this.
                        # Note: less often, read needs may not have enough remaining seq!
                        jnc_seq = seq[:r.reference_length + (lig.junc_len - lig.site_len)]
                        if jnc_seq.endswith(lig.junction):
                            report.enzyme[lig.enzyme_name].read_thru += 1

                            sa_dict = parse_secondary_alignment_tag(r)
                            if sa_dict is not None:
                                # TODO either use this information to stop parsing.
                                report.enzyme[lig.enzyme_name].is_split += 1

                        break

                if not found_lig:
                    report.global_info.no_site += 1

    logger.info('Number of analysed reads: {}'.format(report.global_info.total))
    logger.info('Number of low mapq reads: {}'.format(report.global_info.low_mapq))
    logger.info('Number of cis-mapping pairs: {}'.format(all_cis_count))

    logger.info('Number of short-range cis-mapping pairs: {}'.format(short_count))
    close_avg = short_sum / short_count
    logger.info('Empirical short-range mean differs from \"mean_insert\" by {:.1f}%'
                .format((close_avg - mean_insert) / mean_insert * 100))
    if short_count > 10 * mean_insert:
        logger.info('Mean and median of short-range pair separation: {:.0f}nt {:.0f}nt'
                    .format(close_avg, short_median_est))
    else:
        logger.warning('Too few short pairs, streaming median estimation would be unreliable')
        logger.info('Mean of short-range pair separation: {:.0f}nt'.format(close_avg))

    # line-up the fields for ease of reading.
    field_width = np.maximum(7, np.log10(np.maximum(long_bins, long_counts)).astype(int) + 1)
    logger.info('Long-range distance intervals: {:>{w[0]}d}nt, {:>{w[1]}d}nt, {:>{w[2]}d}nt'
                .format(*long_bins, w=field_width))
    logger.info('Number of cis-mapping pairs:   {:>{w[0]}d},   {:>{w[1]}d},   {:>{w[2]}d}'
                .format(*long_counts, w=field_width))
    prop = long_counts.astype(np.float) / all_cis_count
    logger.info('Relative proportion:           {:#>{w[0]}.4g},   {:#>{w[1]}.4g},   {:#>{w[2]}.4g}'
                .format(*prop, w=field_width))

    print_report(report, sep)
