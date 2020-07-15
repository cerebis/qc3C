__version__ = '0.5rc1'

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


def date_stamp() -> str:
    """
    :return: Return a datetime stamp in the format YYYY-mm-dd hh:mm:ss.f
    """
    from datetime import datetime
    _now = datetime.now()
    return _now.strftime('%Y-%m-%d %H:%M:%S.%f')


def runtime_info() -> dict:
    """
    :return: Return runtime info (version and date stamps) as a dict
    """
    return {'qc3C_version': version_stamp(False), 'run_timestamp': date_stamp()}
