import re
import numpy as np

from leven import levenshtein
from collections import namedtuple
from Bio.Restriction import Restriction
from qc3C.exceptions import *

# immutable type used in storing information about enzymatic byproducts in proximity ligation
LigationInfo = namedtuple('ligation_info', ('enzyme_name', 'junction', 'cut_site', 'junc_len', 'site_len'))


def leven_ratio(a: str, b: str):
    """
    Compute the levenshtein distance between two strings and return
    the ratio relative to the combined string length
    :param a: the first string to compare with b
    :param b: the second string to compare with a
    :return: a value between 0 and 1, 1 being a perfect match (0 distance)
    """
    d = levenshtein(a, b)
    lsum = len(a) + len(b)
    return (lsum - d) / lsum


def get_enzyme_instance(enz_name):
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
    assert not enz.is_ambiguous(), 'ambiguous symbols in enzymatic site not supported'

    end5, end3 = '', ''
    site = str(enz.site)

    ovhg_size = abs(enz.ovhg)
    if ovhg_size > 0 and ovhg_size != enz.size:
        a = abs(enz.fst5)
        if a > enz.size // 2:
            a = enz.size - a
        end5, end3 = enz.site[:a], enz.site[-a:]
        site = site[:-a]
    junc = '{0}{3}{1}{3}{1}{3}{2}'.format(end5, enz.ovhgseq, end3, spacer)
    return LigationInfo(str(enz), junc.upper(), site.upper(), len(junc), len(site))
