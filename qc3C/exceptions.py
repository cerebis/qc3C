class UnknownEnzymeException(Exception):
    """All sequences were excluded during filtering"""
    def __init__(self, target, similar):
        super(UnknownEnzymeException, self).__init__(
            '{} is undefined, but its similar to: {}'.format(target, ', '.join(similar)))
