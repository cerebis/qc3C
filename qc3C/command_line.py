import argparse
import json
import multiprocessing
import os
import re
import subprocess
from collections import namedtuple
from collections.abc import Collection
from difflib import SequenceMatcher

import pysam
import tqdm
from Bio.Restriction import Restriction
from recordclass import recordclass

__version__ = '0.1'


ligation_info = namedtuple('ligation_info', ('enzyme_name', 'junction', 'end_match', 'junc_len', 'site_len'))

pair_info = namedtuple('pair_info', ('name', 'pos', 'length', 'is_reverse', 'cigarstring'))

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


def parse_tag(_tag):
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


class UnknownEnzymeException(Exception):
    """All sequences were excluded during filtering"""
    def __init__(self, target, similar):
        super(UnknownEnzymeException, self).__init__(
            '{} is undefined, but its similar to: {}'.format(target, ', '.join(similar)))


def get_enzyme_instance(enz_name):
    """
    Fetch an instance of a given restriction enzyme by its name.

    :param enz_name: the case-sensitive name of the enzyme
    :return: RestrictionType the enzyme instance
    """
    try:
        # this has to match exactly
        return getattr(Restriction, enz_name)
    except AttributeError:
        # since we're being strict about enzyme names, be helpful with errors
        similar = []
        for a in dir(Restriction):
            score = SequenceMatcher(None, enz_name.lower(), a.lower()).ratio()
            if score >= 0.8:
                similar.append(a)
        raise UnknownEnzymeException(enz_name, similar)


COMPLEMENT_TABLE = str.maketrans('acgtumrwsykvhdbnACGTUMRWSYKVHDBN',
                                 'TGCAAnnnnnnnnnnnTGCAANNNNNNNNNNN')


def revcomp(seq):
    """
    Reverse complement a string representation of a sequence. This uses string.translate.
    :param seq: input sequence as a string
    :return: revcomp sequence as a string
    """
    return seq.translate(COMPLEMENT_TABLE)[::-1]


def ligation_junction_seq(enz, spacer=''):
    """
    Determine the sequence presented after successful enzymatic cleavage and ligation. Due
    to the use of enzymes which possess non-zero overhang and the subsequent end-fill step
    the sequence intervening the cut points gets duplicated.

    This method returns the full junction sequence, containing the 3' and 5' residuals
    and the intervening duplication.

    end5 - dup - end3

    :params enz: biopython restriction instance
    :params spacer: optional string with which to separate site elements (debugging)
    """
    end5, end3 = '', ''
    site = str(enz.site)
    if enz.ovhg % enz.size != 0:
        end5, end3 = enz.site[:enz.fst5], enz.site[enz.fst3:]
        site = site[:enz.fst3]
    junc = '{0}{3}{1}{3}{1}{3}{2}'.format(end5, enz.ovhgseq, end3, spacer)
    return ligation_info(str(enz), junc, site, len(junc), len(site))


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


global_info = recordclass('global_info', ('ref_term', 'no_site', 'full_align', 'early_term'), defaults=(0,)*4)
cutsite_info = recordclass('cutsite_info', ('cs_term', 'cs_full', 'read_thru', 'is_split'), defaults=(0,)*4)


class qc_info(object):

    def __init__(self, enzymes):
        if not isinstance(enzymes, Collection) or isinstance(enzymes, str):
            enzymes = [enzymes]
        self._global = global_info()
        self.enzyme = {en: cutsite_info() for en in enzymes}

    def total(self):
        return sum(self._global)


class EncodeCounter(json.JSONEncoder):
    def default(self, o):
        print(o.__class__, o.__class__.__bases__)
        if isinstance(o, qc_info):
            d = o._asdict()
            d.update({'enzyme': o.enzyme})
            return d
        if isinstance(o, cutsite_info):
            return o._asdict()
        raise TypeError()


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--log', action='store_true', default=False, help='Keep a log of pairs with cut-sites')
    parser.add_argument('-t', '--threads', metavar='N', type=int, default=None, help='Number of threads')
    parser.add_argument('-e', '--enzyme', metavar='NEB_NAME', required=True, action='append',
                        help='Case-sensitive NEB enzyme name. Use multiple times for multiple enzymes')
    parser.add_argument('BAM', help='Input bam file of Hi-C reads mapped to references')
    args = parser.parse_args()

    report = qc_info(args.enzyme)

    ligation_variants = []
    for ename in args.enzyme:
        ligation_variants.append(ligation_junction_seq(get_enzyme_instance(ename)))

    with pysam.AlignmentFile(args.BAM, 'rb') as bam_file:

        ref_lengths = [li for li in bam_file.lengths]
        ref_names = [ni for ni in bam_file.references]

        pair_log = {}

        rmate = None
        prev_r = None

        progress = None
        try:

            print('Counting reads in bam file...')
            n_reads = count_bam_reads(args.BAM, args.threads)

            print('Analyzing bam file...')
            progress = tqdm.tqdm(total=n_reads)
            bam_iter = bam_file.fetch(until_eof=True)
            while True:

                try:
                    r = next(bam_iter)
                    progress.update()

                    if args.log:
                        if prev_r is not None:
                            if r.query_name == prev_r.query_name:
                                rmate = prev_r
                            else:
                                rmate = None
                        prev_r = r

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

                            if args.log:
                                # record those pairs which at least have one cut-site terminated read
                                pair_log.setdefault(r.query_name, set()).add(
                                    pair_info(ref_names[r.reference_id], r.reference_start, r.reference_length,
                                              r.is_reverse, r.cigarstring))
                                if rmate is not None :
                                    if rmate.query_name != r.query_name:
                                        break
                                    pair_log[r.query_name].add(
                                        pair_info(ref_names[rmate.reference_id], rmate.reference_start,
                                                  rmate.reference_length, rmate.is_reverse, rmate.cigarstring))

                            # a proximity ligation product should contain a characteristic
                            # sequence which duplicates a portion of the cut-site. Check
                            # and see if the read contains this.
                            # Note: less often, read needs may not have enough remaining seq!
                            jseq = seq[:r.reference_length + (lig.junc_len - lig.site_len)]
                            if jseq.endswith(lig.junction):
                                report.enzyme[lig.enzyme_name].read_thru += 1

                                if r.has_tag('SA'):
                                    report.enzyme[lig.enzyme_name].is_split += 1
                                    # sec_align = parse_tag(r.get_tag('SA'))

                            break

                    if not found_lig:
                        report._global.no_site += 1

                except StopIteration as e:
                    break
        finally:
            if progress:
                progress.close()

        total = report.total()

        print()
        print('Total reads analyzed: {}'.format(total))
        print()
        n = report._global.full_align
        print('Fully aligned:     {} {:6.2f}%'.format(n, n / total * 100))
        n = report._global.early_term
        print('Early termination: {} {:6.2f}%'.format(n, n / total * 100))
        n = report._global.ref_term
        print('Ref termination:   {} {:6.2f}%'.format(n, n / total * 100))
        n = report._global.no_site
        print('No cut-site:       {} {:6.2f}%'.format(n, n / total * 100))
        for en in report.enzyme:
            n = report.enzyme[en].cs_full
            print('    {} 3p cut-site: {} {:6.2f}%'.format(en, n, n / total * 100))
        print()
        for en, inf in report.enzyme.items():
            print('Suspected {} ligation products'.format(str(en)))
            print('  fully with 3p:   {} {:6.2f}%'.format(inf.cs_full, inf.cs_full / total * 100))
            print('    3p cut-site:   {} {:6.2f}%'.format(inf.cs_term, inf.cs_term / total * 100))
            print('    3p junction:   {} {:6.2f}%'.format(inf.read_thru, inf.read_thru / total * 100))
            print('    split align:   {} {:6.2f}%'.format(inf.is_split, inf.is_split / total * 100))
        print()

    if args.log:
        print('Writing pairs which contained cut-sites to log...')
        json.dump(pair_log, open('pairs_log.toml', 'w'))
