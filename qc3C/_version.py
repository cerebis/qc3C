__version__ = '0.2.6.1'

__copyright__ = """Copyright (C) 2019 Matthew Z DeMaere
This is free software.  You may redistribute copies of it under the terms of
the GNU Affero General Public License <https://www.gnu.org/licenses/agpl.html>.
There is NO WARRANTY, to the extent permitted by law.
"""


def version_stamp(full: bool = True) -> str:
    """
    Create a string indicating the version and possibly extended details such as copyright
    :param full: when True add extended details (multi-line)
    :return: a version stamp string
    """
    if full:
        return 'qc3C {}\n{}'.format(__version__, __copyright__)
    else:
        return 'qc3C {}'.format(__version__)
