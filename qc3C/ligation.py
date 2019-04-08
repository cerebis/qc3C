from collections import namedtuple
from difflib import SequenceMatcher
from Bio.Restriction import Restriction
from qc3C.exceptions import *

# immutable type used in storing information about enzymatic byproducts in proximity ligation
LigationInfo = namedtuple('ligation_info', ('enzyme_name', 'junction', 'end_match', 'junc_len', 'site_len'))


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
    return LigationInfo(str(enz), junc, site, len(junc), len(site))
