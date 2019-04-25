class UnknownEnzymeException(Exception):
    """All sequences were excluded during filtering"""
    def __init__(self, target, similar):
        super(UnknownEnzymeException, self).__init__(
            '{} is undefined, but its similar to: {}'.format(target, ', '.join(similar)))


class NameSortingException(Exception):
    """Bam is not sorted by read name"""
    def __init__(self, filename):
        super(NameSortingException, self).__init__(
            '{} must be sorted by name'.format(filename))
