import numpy as np
import logging

logger = logging.getLogger(__name__)


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
