import numpy as np
import logging
import os

logger = logging.getLogger(__name__)


def warn_if(_test: bool):
    """
    Conditional help for logging warnings.
    :param _test: if true, level is that of WARNING
    :return: return WARNING if true, otherwise INFO
    """
    if _test:
        return logging.WARNING
    else:
        return logging.INFO


def test_for_exe(prog_name):
    """
    Test whether a program exists on the system. This can be either a full path
    or just the executable name. This is case sensitive.
    :param prog_name: the program name with/without path
    :return: the full path or None
    """

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(prog_name)
    if fpath:
        if is_exe(prog_name):
            return prog_name
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, prog_name)
            if is_exe(exe_file):
                return exe_file

    return None


def init_random_state(seed: int = None):
    """
    Helper function to initialise a random state using a given seed, or
    from a random seed. The used seed value is sent to log.
    :param seed: an integer or None. If none, a randomly selected seed is used.
    :return: an initialised numpy.random.RandomState
    """
    # set up random number generation
    if seed is None:
        seed = np.random.randint(1000000, 99999999)
        logger.info('Random seed was not set, using {}'.format(seed))
    else:
        logger.info('Random seed was {}'.format(seed))
    return np.random.RandomState(seed)
