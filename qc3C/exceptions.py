class UnknownEnzymeException(Exception):
    """Supplied enzyme name was not found in biopython"""
    def __init__(self, target, similar):
        super(UnknownEnzymeException, self).__init__(
            '{} is undefined, but its similar to: {}'.format(target, ', '.join(similar)))


class InvalidEnzymeException(Exception):
    """Supplied enzyme name contains invalid characters"""
    def __init__(self, target):
        super(InvalidEnzymeException, self).__init__(
            'The enzyme name \"{}\" contains invalid characters'.format(target))


class NameSortingException(Exception):
    """Bam is not sorted by read name"""
    def __init__(self, filename):
        super(NameSortingException, self).__init__(
            '{} must be sorted by name'.format(filename))
