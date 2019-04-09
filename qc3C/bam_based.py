import json
import multiprocessing
import os
import re
import subprocess
import pandas
import pysam
import tqdm
import logging

from recordclass import recordclass
from collections import namedtuple
from collections.abc import Collection
from qc3C.ligation import ligation_junction_seq, get_enzyme_instance

logger = logging.getLogger(__name__)

# immutable type used in pair logging
PairInfo = namedtuple('pair_info', ('name', 'pos', 'length', 'is_reverse', 'cigarstring'))

# mutable tuples used in storing counting statistics on evidence of proximity ligation events.
# information pertaining to the entire read-set
GlobalInfo = recordclass('global_info',
                         ('ref_term', 'no_site', 'full_align', 'early_term', 'total'),
                         defaults=(0,) * 5)

# information specific an enzyme's action
CutSiteInfo = recordclass('cutsite_info',
                          ('cs_term', 'cs_full', 'read_thru', 'is_split'),
                          defaults=(0,) * 4)


class QcInfo(object):
    """
    Represents the evidence collected for the action of Hi-C proximity ligation
    over an entire BAM file.
    """

    class EncodeCounter(json.JSONEncoder):
        def default(self, o):
            print(o.__class__, o.__class__.__bases__)
            if isinstance(o, QcInfo):
                d = o._asdict()
                d.update({'enzyme': o.enzyme})
                return d
            if isinstance(o, CutSiteInfo):
                return o._asdict()
            raise TypeError()

    def __init__(self, enzymes):
        """
        :param enzymes: the enzymes used in creating the Hi-C library
        """
        if not isinstance(enzymes, Collection) or isinstance(enzymes, str):
            enzymes = [enzymes]
        self._global = GlobalInfo()
        self.enzyme = {en: CutSiteInfo() for en in enzymes}

    def total(self):
        """
        :return: the total number of reads
        """
        return self._global.total

    def mapped(self):
        """
        :return: the total number of mapped reads
        """
        return self._global.ref_term + self._global.full_align + self._global.early_term


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


def count_bam_reads(file_name, max_cpu=None):
    """
    Use samtools to quickly count the number of non-header lines in a bam file. This is assumed to equal
    the number of mapped reads.
    :param file_name: a bam file to scan (neither sorted nor an index is required)
    :param max_cpu: set the maximum number of CPUS to use in counting (otherwise all cores)
    :return: estimated number of mapped reads
    """
    assert exe_exists('samtools'), 'required tool samtools was not found on path'
    assert exe_exists('wc'), 'required tool wc was not found on path'
    if max_cpu is None:
        max_cpu = multiprocessing.cpu_count()
    proc = subprocess.Popen(['samtools', 'view', '-c', '-@{}'.format(max_cpu), file_name],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return int(proc.stdout.readline())


def print_report(report, output=None, sep='\t'):
    """
    Create a tabular report in CSV format
    :param report: the QcInfo instance to report
    :param output: if None print the table to stdout, otherwise to the specified path
    :param sep: separator to use in table creation
    """

    total = report.mapped()
    _cnames = ['enzyme', 'variable', 'read_count', 'vs_total', 'vs_mapped']

    df = pandas.DataFrame([['all', 'total', total, None, None],
                           ['all', 'mapped', report.mapped(), None, None],
                           ['all', 'fully aligned', report._global.full_align, None, None],
                           ['all', 'early termination', report._global.early_term, None, None],
                           ['all', 'ref termination', report._global.ref_term, None, None],
                           ['all', 'no cutsite', report._global.no_site, None, None]],
                          columns=_cnames)

    for enz, inf in report.enzyme.items():
        df_enz = pandas.DataFrame([[enz, 'cs fully aligned', inf.cs_full, None, None],
                                   [enz, 'cs termination', inf.cs_term, None, None],
                                   [enz, 'read-thru', inf.read_thru, None, None],
                                   [enz, 'split alignment', inf.is_split, None, None]],
                                  columns=_cnames)
        df = df.append(df_enz, ignore_index=True)

    df.loc[1:, 'vs_total'] = df.loc[1:, 'read_count'] / df.loc[0, 'read_count'] * 100
    df.loc[2:, 'vs_mapped'] = df.loc[2:, 'read_count'] / df.loc[1, 'read_count'] * 100
    df[['vs_total', 'vs_mapped']] = df[['vs_total', 'vs_mapped']].apply(pandas.to_numeric)

    if output is None:
        print(df.to_csv(None, sep=sep, float_format="%.2f", index=False))
    else:
        df.to_csv(output, sep=sep, float_format="%.2f", index=False)


def analyze(bam_file, enzymes, output, threads=1, sep='\t'):

    report = QcInfo(enzymes)

    ligation_variants = []
    for ename in enzymes:
        ligation_variants.append(ligation_junction_seq(get_enzyme_instance(ename)))

    with pysam.AlignmentFile(bam_file, 'rb', threads=threads) as bam:

        ref_lengths = [li for li in bam.lengths]
        # ref_names = [ni for ni in bam.references]

        # pair_log = {}

        # rmate = None
        # prev_r = None

        progress = None
        try:

            logger.info('Counting reads in bam file...')
            n_reads = count_bam_reads(bam_file, threads)

            logger.info('Analyzing bam file...')
            progress = tqdm.tqdm(total=n_reads)
            bam_iter = bam.fetch(until_eof=True)
            while True:

                try:
                    r = next(bam_iter)
                    progress.update()

                    # report._global.total += 1
                    # if r.is_unmapped:
                    #     continue

                    # if args.log:
                    #     if prev_r is not None:
                    #         if r.query_name == prev_r.query_name:
                    #             rmate = prev_r
                    #         else:
                    #             rmate = None
                    #     prev_r = r

                    if r.reference_end >= ref_lengths[r.reference_id] or r.reference_start == 0:
                        # reads which align to the ends of references are ignored
                        report._global.ref_term += 1
                        continue

                    if r.query_length == r.reference_length:
                        report._global.full_align += 1
                        # fully aligned reads can't be tested for the junction but
                        # we can still test for the cut-site
                        seq = get_forward_strand(r)
                        for lig in ligation_variants:
                            if seq.endswith(lig.end_match):
                                report.enzyme[lig.enzyme_name].cs_full += 1
                                break
                        continue

                    report._global.early_term += 1

                    # inspect all sequences in as 5'-3'
                    seq = get_forward_strand(r)

                    # the aligned sequence, which should end with cut-site
                    aseq = seq[:r.reference_length]

                    # check that a cut-site exists on the end
                    found_lig = False
                    for lig in ligation_variants:
                        if aseq.endswith(lig.end_match):
                            found_lig = True
                            report.enzyme[lig.enzyme_name].cs_term += 1

                            # if args.log:
                            #     # record those pairs which at least have one cut-site terminated read
                            #     pair_log.setdefault(r.query_name, set()).add(
                            #         pair_info(ref_names[r.reference_id], r.reference_start, r.reference_length,
                            #                   r.is_reverse, r.cigarstring))
                            #     if rmate is not None :
                            #         if rmate.query_name != r.query_name:
                            #             break
                            #         pair_log[r.query_name].add(
                            #             pair_info(ref_names[rmate.reference_id], rmate.reference_start,
                            #                       rmate.reference_length, rmate.is_reverse, rmate.cigarstring))

                            # a proximity ligation product should contain a characteristic
                            # sequence which duplicates a portion of the cut-site. Check
                            # and see if the read contains this.
                            # Note: less often, read needs may not have enough remaining seq!
                            jseq = seq[:r.reference_length + (lig.junc_len - lig.site_len)]
                            if jseq.endswith(lig.junction):
                                report.enzyme[lig.enzyme_name].read_thru += 1

                                sa_dict = parse_secondary_alignment_tag(r)
                                if sa_dict is not None:
                                    # TODO either use this inforation to stop parsing.
                                    report.enzyme[lig.enzyme_name].is_split += 1

                            break

                    if not found_lig:
                        report._global.no_site += 1

                except StopIteration as e:
                    break
        finally:
            if progress:
                progress.close()

    print_report(report, output, sep)

    # if args.log:
    #     logger.info('Writing pairs which contained cut-sites to log...')
    #     json.dump(pair_log, open('pairs_log.toml', 'w'))