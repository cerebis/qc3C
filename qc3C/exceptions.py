class ApplicationException(Exception):
    """Base exception class"""
    def __init__(self, message):
        super(ApplicationException, self).__init__(message)


class MaxObsLimit(StopIteration):
    """Maximum observation limit reached"""
    pass


class ZeroCoverageException(ApplicationException):
    """Queried observations of K-mer coverage summed to zero"""
    def __init__(self):
        super(ZeroCoverageException, self).__init__(
            'Queried observations of K-mer coverage summed to zero')


class UnknownLibraryKitException(ApplicationException):
    """The library kit is unknown"""
    def __init__(self, library_kit: str):
        super(UnknownLibraryKitException, self).__init__(
            'The library kit {} is unknown'.format(library_kit))


class UnknownEnzymeException(ApplicationException):
    """Supplied enzyme name was not found in biopython"""
    def __init__(self, target: str, similar: list):
        super(UnknownEnzymeException, self).__init__(
            '{} is undefined, but its similar to: {}'.format(target, ', '.join(similar)))


class InvalidEnzymeException(ApplicationException):
    """Supplied enzyme name contains invalid characters"""
    def __init__(self, target: str):
        super(InvalidEnzymeException, self).__init__(
            'The enzyme name \"{}\" contains invalid characters'.format(target))


class NameSortingException(ApplicationException):
    """Bam is not sorted by read name"""
    def __init__(self, filename: str):
        super(NameSortingException, self).__init__(
            '{} must be sorted by name'.format(filename))


class InsufficientDataException(ApplicationException):
    """Not enough observational data was obtained for subsequent analysis"""
    def __init__(self, reason: str):
        super(InsufficientDataException, self).__init__(reason)
