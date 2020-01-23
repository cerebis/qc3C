import contextlib
import gzip
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
from collections import namedtuple
from collections.abc import MutableMapping
from leven import levenshtein
from numba import jit
from typing import TypeVar, Tuple, Optional

from .exceptions import *
from .utils import open_input, modification_hash, count_sequences

logger = logging.getLogger(__name__)

EnzymeType = TypeVar('EnzymeType',
                     Restriction.RestrictionType,
                     Restriction.AbstractCut,
                     Restriction.Ov5,
                     Restriction.Ov3,
                     Restriction.Defined)

# immutable type used in storing information about enzymatic byproducts in proximity ligation
LigationInfo = namedtuple('ligation_info', ('enzyme_name', 'junction', 'cut_site', 'vestigial',
                                            'junc_len', 'site_len', 'vest_len', 'elucidation'))


def digest_sequence(enzyme: EnzymeType, seq: SeqRecord, linear: bool = True) -> np.ndarray:
    """
    Determine the restriction sites for a given enzyme and sequence. Return 0-based
    array of sites.

    :param enzyme: the case-sensitive enzyme name or Bio.RestrictionType instance
    :param seq: a Bio.SeqRecord object
    :param linear: treat the sequence as being linear (=True) or circular (=False)
    :return: 0-based array of site locations in ascending genomic coordinates
    """
    if isinstance(enzyme, str):
        enzyme = get_enzyme_instance(enzyme)
    sites = np.array(enzyme.search(seq.seq, linear=linear), dtype=np.int32) - 1
    # although biopython should return sites in ascending order, lets be pedantic
    return np.sort(sites)


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

    def build_db(self):

        with contextlib.closing(Bio.SeqIO.parse(open_input(self.fasta_file), 'fasta')) as seq_iter:

            logger.info('Building cut-site database from supplied enzyme and reference sequences')

            if self.use_tqdm:
                total_seq = count_sequences(self.fasta_file, 'fasta', 4)
                seq_iter = tqdm.tqdm(seq_iter, total=total_seq)

            cutsite_db = {}
            for seq in seq_iter:
                cutsite_db[seq.id] = CutSites(seq, self.enzyme)
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
        return '{}_sites_{}'.format(self.fasta_file, str(self.enzyme).lower())

    def __init__(self, enzyme, fasta_file, use_cache=False, use_tqdm=False):

        self.enzyme = enzyme
        self.fasta_file = fasta_file
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
                logger.info('No cut-site database cache found')
                self.build_db()
                self.save_cache()
        else:
            self.build_db()

    def __getitem__(self, item):
        return self.data[item]

    def __setitem__(self, key, value):
        self.data[key] = CutSites(value, self.enzyme)

    def __delitem__(self, key):
        del self.data[key]

    def __iter__(self):
        return self.data.__iter__()

    def __len__(self):
        return self.data.__len__()


class CutSites(object):

    def __init__(self, seq: SeqRecord, enzyme: EnzymeType):
        """
        Instances of this class contain information on the restriction
        digest of the supplied sequence. This can be queried for nearby
        sites, or information on sites related to the given pair.

        :param seq: a SeqRecord object
        :param enzyme: a Bio.Restriction enzyme or NEB-based enzyme name of such
        """
        self.fwd_sites = digest_sequence(enzyme, seq)
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


def get_enzyme_instance(enz_name: str) -> Restriction.RestrictionType:
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


def ligation_junction_seq(enz: EnzymeType, spacer: str = '', no_ambig: bool = True) -> LigationInfo:
    """
    Determine the sequence presented after successful enzymatic cleavage and ligation. Due
    to the use of enzymes which possess non-zero overhang and the subsequent end-fill step
    the sequence intervening the cut points gets duplicated.

    This method returns the full junction sequence, containing the 3' and 5' residuals
    and the intervening duplication.

    end5 - dup - end3

    :params enz: biopython restriction instance
    :params spacer: optional string with which to separate site elements (debugging)
    :params no_ambig: throws an exception when true and an enzyme contains ambiguous symbols in its recognition site
    """
    def vestigial_site():
        """
        Determine the part of the cut-site that will remain when a ligation junction is
        created. This can be the entire cut-site or a smaller portion, beginning from
        the left.
        """
        i = 0
        while i < enz.size and enz.site[i] == junc[i]:
            i += 1
        return str(enz.site[:i])

    assert not enz.is_blunt(), 'enzyme produces blunt ends'
    if no_ambig:
        assert not enz.is_ambiguous(), 'ambiguous symbols in enzymatic site not supported'

    end5, end3 = '', ''

    ovhg_size = abs(enz.ovhg)
    if ovhg_size > 0 and ovhg_size != enz.size:
        a = abs(enz.fst5)
        if a > enz.size // 2:
            a = enz.size - a
        end5, end3 = enz.site[:a], enz.site[-a:]
    junc = '{0}{3}{1}{3}{1}{3}{2}'.format(end5, enz.ovhgseq, end3, spacer)
    vest = vestigial_site()
    return LigationInfo(str(enz), junc.upper(), enz.site.upper(), vest.upper(),
                        len(junc), enz.size, len(vest), enz.elucidate())


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
