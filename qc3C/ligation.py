import collections
import contextlib
import gzip
import itertools
import logging
import numpy as np
import os
import pickle
import pysam
import re
import tqdm
import Bio.SeqIO

from Bio.SeqRecord import SeqRecord
from Bio.Restriction import Restriction

from collections.abc import MutableMapping
from leven import levenshtein
from numba import jit
from typing import Dict, List, Optional, Tuple, TypeVar, Union

from .exceptions import *
from .utils import open_input, modification_hash, count_sequences

logger = logging.getLogger(__name__)

EnzymeType = TypeVar('EnzymeType',
                     Restriction.RestrictionType,
                     Restriction.AbstractCut,
                     Restriction.Ov5,
                     Restriction.Ov3,
                     Restriction.Defined)

# Information pertaining to an restriction enzyme
enzyme_info = collections.namedtuple('enzyme_info',
                                     ('name', 'site', 'site_len'))

# Information pertaining to a Hi-C ligation junction
ligation_info = collections.namedtuple('ligation_info',
                                       ('enz5p', 'enz3p', 'junction', 'vestigial', 'junc_len',
                                        'vest_len', 'pattern', 'cross'))

# Information pertaining to digestion / private class for readability
digest_info = collections.namedtuple('digest_info', ('enzyme', 'end5', 'end3', 'overhang'))


class Digest(object):
    """
    A Hi-C digestion. This can support one or two enzyme digestions and ambiguous nucleotides within
    the cut-site such as used by Arima.
    """
    def __init__(self, enzyme_a: str, enzyme_b: str = None, no_ambig: bool = True):

        def vestigial_site(enz: EnzymeType, junc: str) -> str:
            """
            Determine the part of the 5-prime end cut-site that will remain when a ligation junction is
            created. Depending on the enzymes involved, this may be the entire cut-site or a smaller portion
            and begins from the left.
            """
            i = 0
            while i < enz.size and (enz.site[i] == junc[i] or enz.site[i] == 'N'):
                i += 1
            return str(junc[:i]).upper()

        def enzyme_ends(enzyme: EnzymeType) -> Tuple[str, str]:
            """
            Determine the 5` and 3` ends of a restriction endonuclease. Here, "ends"
            are the parts of the recognition site not involved in overhang.
            :param enzyme: Biopython ins
            :return: a pair of strings containing the 5' and 3' ends
            """
            end5, end3 = '', ''
            ovhg_size = abs(enzyme.ovhg)
            if ovhg_size > 0 and ovhg_size != enzyme.size:
                a = abs(enzyme.fst5)
                if a > enzyme.size // 2:
                    a = enzyme.size - a
                end5, end3 = enzyme.site[:a], enzyme.site[-a:]
            return end5.upper(), end3.upper()

        def junction_duplication(enzyme_a: EnzymeType, enzyme_b: Optional[EnzymeType] = None) -> Dict[str, ligation_info]:
            """
            For the enzyme cocktail, generate the set of possible ligation junctions.
            :param enzyme_a: the first enzyme
            :param enzyme_b: the optional second enzyme
            :return: a dictionary of ligation_info objects
            """
            end5, end3 = enzyme_ends(enzyme_a)
            enz_list = [digest_info(enzyme_a, end5, end3, enzyme_a.ovhgseq.upper())]

            if enzyme_b is not None:
                end5, end3 = enzyme_ends(enzyme_b)
                enz_list.append(digest_info(enzyme_b, end5, end3, enzyme_b.ovhgseq.upper()))

            _junctions = {}
            for a, b in itertools.product(enz_list, repeat=2):
                # the actual sequence
                _junc_seq = f'{a.end5}{a.overhang}{b.overhang}{b.end3}'
                if _junc_seq not in _junctions:
                    _vest_seq = vestigial_site(a.enzyme, _junc_seq)
                    _junctions[_junc_seq] = ligation_info(str(a.enzyme), str(b.enzyme),
                                                          _junc_seq, _vest_seq,
                                                          len(_junc_seq), len(_vest_seq),
                                                          re.compile(_junc_seq.replace('N', '.')),
                                                          str(a.enzyme) == str(b.enzyme))
            return _junctions

        # convert the str values to instantiated enzyme objects
        enzyme_a, enzyme_b = get_enzymes(enzyme_a, enzyme_b)

        # test digestion
        assert not enzyme_a.is_blunt(), 'blunt-ended enzymes are not supported'
        if enzyme_b is not None:
            assert not enzyme_b.is_blunt(), 'blunt-ended enzymes are not supported'
            assert enzyme_a != enzyme_b, 'enzymes a and b are either same enzyme or equisoschizomers'

        if no_ambig:
            assert not enzyme_a.is_ambiguous(), 'ambiguous symbols in enzymatic site not supported'
            if enzyme_b is not None:
                assert not enzyme_b.is_ambiguous, 'ambiguous symbols in enzymatic site not supported'

        # target cut-sites
        self.cutsites = {enz.site.upper(): enzyme_info(str(enz), enz.site.upper(), len(enz.site))
                         for enz in [enzyme_a, enzyme_b] if enz is not None}

        # hi-c ligation junctions
        self.junctions = junction_duplication(enzyme_a, enzyme_b)

        # prepare regex based matching methods for the various products from ligation

        # cut-sites for the specified enzymes
        self.any_cutsite = re.compile('({})'.format('|'.join(
            sorted(self.cutsites, key=lambda x: (-len(x), x))).replace('N', '.')))
        self.end_cutsite = re.compile('{}$'.format(self.any_cutsite.pattern))

        # Hi-C ligation junctions, which can be any combination of ends made from enzyme cocktail
        self.any_junction = re.compile('({})'.format('|'.join(
            sorted(self.junctions, key=lambda x: (-len(x), x))).replace('N', '.')))
        self.end_junction = re.compile('{}$'.format(self.any_junction.pattern))

        # Digested ends, which are religated will not necessarily possess the entire
        # cut-site, but will possess a remnant/vestigial sequence
        self.any_vestigial = re.compile('({})'.format('|'.join(
            sorted(set([v.vestigial for v in self.junctions.values()]), key=lambda x: (-len(x), x))).replace('N', '.')))
        self.end_vestigial = re.compile('{}$'.format(self.any_vestigial.pattern))

    def unique_vestigial(self) -> List['Digest.DerivedString']:
        uniq_vest = set([(ji.vestigial, ji.enz5p) for ji in self.junctions.values()])
        return [Digest.DerivedString(_s, _e) for _s, _e in uniq_vest]

    def longest_junction(self) -> int:
        return max(li.junc_len for li in self.junctions.values())

    def unambiguous_junctions(self) -> List['Digest.DerivedString']:
        juncs = []
        for ji in self.junctions.values():
            juncs.append(Digest.DerivedString(ji.junction, '{}/{}'.format(ji.enz5p, ji.enz3p)))
        return Digest._expand_ambiguous(*juncs)

    def unambiguous_vestigial(self) -> List['Digest.DerivedString']:
        return Digest._expand_ambiguous(*self.unique_vestigial())

    def tracker(self, kmer_type: str) -> Dict['Digest.DerivedString', int]:
        """
        A dictionary object ready for tracking the number of observations for
        the possible k-mers produced through the specified Hi-C digestion.
        :param kmer_type: either "junction" or "vestigial" k-mers
        :return:
        """
        if kmer_type == 'junction':
            return collections.OrderedDict({ji: 0 for ji in self.unambiguous_junctions()})
        elif kmer_type == 'vestigial':
            return collections.OrderedDict({vi: 0 for vi in self.unambiguous_vestigial()})
        else:
            raise ValueError(f'unsupported kmer_type {kmer_type}')

    @staticmethod
    def gather_tracker(tracker: Dict['Digest.DerivedString', int]) -> Dict[str, Dict[str, int]]:
        """
        Gather tracker results by enzyme participants
        :param tracker: a tracker dictionary produced from Digest.tracker method
        :return: a dict of dicts containing counts per k-mer organised by enzyme
        """
        d = collections.defaultdict(dict)
        for _s, _n in tracker.items():
            d[_s.source][str(_s)] = _n
        return d

    class DerivedString(str):
        """
        A subclass of str which retains a record of its source.
        """
        def __new__(cls, _str: str, _source: str = None):
            """
            :param _str: the string value
            :param _source: the source str from which this string is derived
            """
            obj = str.__new__(cls, _str)
            obj.source = _source
            return obj

        def to_unambiguous(self) -> List['Digest.DerivedString']:
            """
            Generate the list of unambiguous strings from this instance
            :return: a list of derived strings
            """
            _str = str(self)
            _src = _str if self.source is None else self.source

            n_count = _str.count('N')
            if n_count == 0:
                enum_seq = [Digest.DerivedString(_str, _src)]
            else:
                enum_seq = [_str]
                for i in range(_str.count('N')):
                    enum_seq = [Digest.DerivedString(_s.replace('N', _c, 1), _src) for _c in 'ACGT' for _s in enum_seq]
            return enum_seq

    @staticmethod
    def _expand_ambiguous(*seqs: 'Digest.DerivedString') -> List['Digest.DerivedString']:
        """
        Enumerate all sequences when ambiguous bases (N only) are encountered.
        :param s: a DerivedString representing DNA, possibly with ambiguous N
        :return: a list of DerivedString representing enumerated unambiguous sequences
        """
        return [_ds for _s in seqs for _ds in _s.to_unambiguous()]

    def cutsite_searcher(self, method):
        if method == 'startswith':
            return self.any_cutsite.match
        elif method == 'find':
            return self.any_cutsite.search
        elif method == 'endswith':
            return self.end_cutsite.search

    def junction_searcher(self, method):
        if method == 'startswith':
            return self.any_junction.match
        elif method == 'find':
            return self.any_junction.search
        elif method == 'endswith':
            return self.end_junction.search

    def vestigial_searcher(self, method):
        if method == 'startswith':
            return self.any_vestigial.match
        elif method == 'find':
            return self.any_vestigial.search
        if method == 'endswith':
            return self.end_vestigial.search


def digest_sequence(seq: SeqRecord, enzymes: List[EnzymeType], linear: bool = True) -> np.ndarray:
    """
    Determine the restriction sites for a given enzyme and sequence. Return 0-based
    array of sites.

    :param seq: a Bio.SeqRecord object
    :param enzymes: list of Bio.RestrictionType instances (1 or 2 only)
    :param linear: treat the sequence as being linear (=True) or circular (=False)
    :return: 0-based array of site locations in ascending genomic coordinates
    """
    # digest the sequence for each enzyme and collect all sites in a single array
    sites = np.array([_s for _enz in enzymes for _s in _enz.search(seq.seq, linear=linear)], dtype=np.int32) - 1
    # reoder the sites in ascending order and remove duplicates
    return np.unique(np.sort(sites))


@jit(nopython=True)
def _fast_nearest_upstream(pos: int, sites: np.ndarray) -> Optional[int]:
    """
   Report the nearest 3-prime cutsite or None. (linear always)
   :param pos: genomic coordinate
   :param sites: an ascending array of sites along the subject genome/contig
   :return: the nearest integer location or None if no down-stream site exists before the end of the sequence
   """
    x = np.searchsorted(sites, pos)
    if x == 0:
        return sites[0]
    elif x == sites.shape[0]:
        return None
    return sites[x]


@jit(nopython=True)
def _fast_count_between(a: int, b: int, sites: np.ndarray) -> int:
    """
    Report the number of sites between a pair of genomic coordinates. It
    is expected that a<b, but in cases where it is not, the order will be
    reversed.
    :param a: first position
    :param b: last position
    :param sites: a list of genomic sites
    :return: the number of sites between a and b (inclusive)
    """
    if a > b:
        a, b = b, a
    return np.where((sites >= a) & (sites <= b))[0].shape[0]


@jit(nopython=True)
def _fast_pair_info(r1_fwd: int, r1_start: int, r1_end: int,
                    r2_fwd: int, r2_start: int, r2_end: int) -> Tuple[int, int, int, bool, bool]:
    """
    Report the separation between a read pair, whether they were exchanged (flipped)
    due to coordinate order (r1<r2), and whether their alignments are overlapping.

    :param r1_fwd: is read_1 forward orientation
    :param r1_start: beginning of read_1 alignment on reference
    :param r1_end: end of read_1 alignment on reference
    :param r2_fwd: is read_2 forward orientation
    :param r2_start: beginning of read_2 alignment on reference
    :param r2_end: end of read_2 alignment on reference
    :return: tuple of r1 position,
                      r2 position,
                      separation,
                      whether r1 and r2 were flipped,
                      whether r1 and r2 overlap
    """
    flipped = False
    if r1_start > r2_start:
        flipped = True
        r1_fwd, r1_start, r1_end, r2_fwd, r2_start, r2_end = r2_fwd, r2_start, r2_end, r1_fwd, r1_start, r1_end
    return r1_end, r2_start, r2_end - r1_start, flipped, r1_end > r2_start


class CutSitesDB(MutableMapping):

    def __init__(self, fasta_file: str, enzyme_names: List[str], use_cache: bool = False, use_tqdm: bool = False):

        self.fasta_file = fasta_file
        # alphabetic order regardless of input order
        self.enzyme_names = sorted(enzyme_names)
        self.use_cache = use_cache
        self.use_tqdm = use_tqdm
        self.data = None
        self.file_hash = modification_hash(fasta_file)

        if use_cache:
            if os.path.exists(self._cache_file()):
                cached_db = self.load_cache()
                if self.file_hash == cached_db.file_hash:
                    self.__dict__.clear()
                    self.__dict__.update(cached_db.__dict__)
                else:
                    logger.info('Cached database appears to be out of date')
                    os.remove(self._cache_file())
                    self.build_db()
                    self.save_cache()
            else:
                logger.info('No cut-site database cache found for digest involving {}'.format(
                    ', '.join(self.enzyme_names)))
                self.build_db()
                self.save_cache()
        else:
            self.build_db()

    def build_db(self):

        with contextlib.closing(Bio.SeqIO.parse(open_input(self.fasta_file), 'fasta')) as seq_iter:

            logger.info('Building cut-site database against reference sequence and {}'.format(
                ', '.join(self.enzyme_names)))

            if self.use_tqdm:
                total_seq = count_sequences(self.fasta_file, 'fasta', 4)
                seq_iter = tqdm.tqdm(seq_iter, total=total_seq)

            enzymes = get_enzymes(*self.enzyme_names)

            cutsite_db = {}
            for seq in seq_iter:
                cutsite_db[seq.id] = CutSites(seq, enzymes)

            self.data = cutsite_db
            self.file_hash = modification_hash(self.fasta_file)

    def load_cache(self):
        with gzip.open(self._cache_file(), 'rb') as cache_h:
            logger.info('Loading cached cut-site database')
            return pickle.load(cache_h)

    def save_cache(self):
        with gzip.open(self._cache_file(), 'wb') as cache_h:
            logger.info('Saving cut-site database for reuse')
            pickle.dump(self, cache_h, pickle.HIGHEST_PROTOCOL)

    def _cache_file(self):
        return '{}_sites_{}'.format(self.fasta_file, '-'.join(str(_en) for _en in self.enzyme_names).lower())

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, key, value):
        self.data[key] = CutSites(value, get_enzymes(*self.enzyme_names))

    def __delitem__(self, key):
        del self.data[key]

    def __iter__(self):
        return self.data.__iter__()

    def __len__(self):
        return self.data.__len__()


class CutSites(object):

    def __init__(self, seq: SeqRecord, enzymes: List[EnzymeType]):
        """
        Instances of this class contain information on the restriction
        digest of the supplied sequence. This can be queried for nearby
        sites, or information on sites related to the given pair.

        :param seq: a SeqRecord object
        :param enzymes: a list of Bio.Restriction enzyme instances (1 or 2 only)
        """
        self.fwd_sites = digest_sequence(seq, enzymes)
        self.rev_sites = -1 * self.fwd_sites[::-1]
        self.n_sites = self.fwd_sites.shape[0]

    def nearest_upstream(self, read):
        """
        Report the nearest cut-site to a reads 3-prime end.

        :param read: a AlignedSegment
        :return: an integer genomic position
        """
        if read.is_reverse:
            x = _fast_nearest_upstream(-read.reference_start, self.rev_sites)
            if x is None:
                return None
            else:
                return -x
        else:
            return _fast_nearest_upstream(read.reference_end, self.fwd_sites)

    def count_sites(self, a, b):
        """
        Count the number of sites between genomic coordinates a and b, inclusive.
        :param a: first position
        :param b: second position
        :return: number of sites
        """
        return _fast_count_between(a, b, self.fwd_sites)

    def pair_info(self, read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> Tuple[int, int, bool]:
        """
        Report information relating to an aligned pair.

        :param read1: first read
        :param read2: second read
        :return: a tuple containing
                   genomic separation,
                   the number of intervening cutsites.
                   whether the read order was flipped,
        """
        a, b, d, flipped, overlapped = _fast_pair_info(not read1.is_reverse, read1.reference_start, read1.reference_end,
                                                       not read2.is_reverse, read2.reference_start, read2.reference_end)
        if overlapped:
            ns = 0
        else:
            ns = _fast_count_between(a, b, self.fwd_sites)
        return d, ns, flipped


def leven_ratio(a: str, b: str) -> float:
    """
    Compute the levenshtein distance between two strings and return
    the ratio relative to the combined string length.

    :param a: the first string to compare with b
    :param b: the second string to compare with a
    :return: a value between 0 and 1, 1 being a perfect match (0 distance)
    """
    d = levenshtein(a, b)
    lsum = len(a) + len(b)
    return (lsum - d) / lsum


def get_enzymes(*enzyme_names: str) -> Union[EnzymeType, List[EnzymeType]]:
    """
    Return multiple instances of Bio.Restriction enzyme classes from
    the string name.

    :param enzyme_names: one or more enzyme names
    :return: a list of instantiated objects
    """

    def get_enzyme_instance(enz_name: str) -> EnzymeType:
        """
        Fetch an instance of a given restriction enzyme by its name.

        :param enz_name: the case-sensitive name of the enzyme
        :return: RestrictionType the enzyme instance
        """
        # try an exact match
        try:
            return getattr(Restriction, enz_name)

        # otherwise, supply a helpful error message
        except AttributeError:
            # check that the name uses only valid characters
            if re.search(r'[^0-9A-Za-z]', enz_name) is not None:
                raise InvalidEnzymeException(enz_name)

            # now look for similar enzyme names to suggest
            lower_ename = enz_name.lower()
            similar = []
            for a in dir(Restriction):
                if a[0].isupper() and a[-1].isupper():
                    similar.append((a, leven_ratio(lower_ename, a.lower())))

            similar = np.array(similar, dtype=np.dtype([('name', 'S20'), ('score', np.float)]))
            top = similar[np.argsort(similar['score'])[-3:]][::-1]

            ix = top['score'] > 0.9
            # if there are near-perfect matches, only report those
            if ix.sum() > 0:
                top = top['name'][ix]
            # otherwise, try and suggest a few maybes
            else:
                top = top['name'][top['score'] > 0.7]
            raise UnknownEnzymeException(enz_name, [s.decode('UTF-8') for s in top])

    return [_en if _en is None else get_enzyme_instance(_en) for _en in enzyme_names]


def restriction_enzyme_by_site(size: int = None) -> dict:
    """
    Create a dictionary of restriction enzyme names known to BioPython, indexed by their recognition site.

    :return: a dict (k=site, v=list of enzyme names)
    """
    by_site = {}
    for cl_name in dir(Restriction):
        # avoid some additional non-enzyme classes within the Restriction module
        if cl_name[0].isupper() and cl_name[-1].isupper():
            try:
                clazz = getattr(Restriction, cl_name)
                if size is not None and clazz.size == size:
                    by_site.setdefault(clazz.site, []).append(str(clazz))
                else:
                    by_site.setdefault(clazz.site, []).append(str(clazz))
            except AttributeError:
                pass
    return by_site
