import logging
import qc3C.bam_based as bam
import qc3C.kmer_based as kmer

__version__ = '0.2.3'
__log_name__ = 'qc3C.log'
__copyright__ = """Copyright (C) 2019 Matthew Z DeMaere
This is free software.  You may redistribute copies of it under the terms of
the GNU Affero General Public License <https://www.gnu.org/licenses/agpl.html>.
There is NO WARRANTY, to the extent permitted by law.
"""


def mk_version():
    return 'qc3C {}\n{}'.format(__version__, __copyright__)


def main():
    import argparse
    import sys

    """
    Shared CLI arguments
    """
    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('-p', '--sample-rate', default=None, type=float,
                               help='Sample only a proportion of all read-pairs [None]')
    global_parser.add_argument('-s', '--seed', type=int,
                               help='Random seed used in sampling the read-set')
    global_parser.add_argument('-t', '--threads', metavar='N', type=int, default=1, help='Number of threads')
    global_parser.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append', required=True,
                               help='One or more case-sensitive NEB enzyme names '
                                    '(Use multiple times for multiple files enzymes)')
    global_parser.add_argument('-m', '--mean-insert', type=int, required=True,
                               help='Mean fragment length to use in estimating the unobserved junction rate')

    parser = argparse.ArgumentParser(description='qc3C: Hi-C quality control')
    parser.add_argument('-V', '--version', default=False, action='store_true', help='Version')
    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')
    subparsers.required = False
    cmd_bam = subparsers.add_parser('bam', parents=[global_parser], description='Alignment-based analysis.')
    cmd_kmer = subparsers.add_parser('kmer', parents=[global_parser], description='Kmer-based analysis.')

    """
    CLI for BAM based analysis
    """
    cmd_bam.add_argument('-b', '--bam', required=True,
                         help='Input name-sorted bam file of Hi-C reads mapped to references')

    """
    CLI for Kmer based analysis
    """
    cmd_kmer.add_argument('--save-cov', default=False, action='store_true',
                          help='Save the collected coverage information to file')
    cmd_kmer.add_argument('-x', '--max-coverage', default=500, type=int,
                          help='Ignore regions with more than this coverage [500]')
    cmd_kmer.add_argument('-k','--kmer-size', type=int, required=True,
                          help='Kmer size used in database')
    cmd_kmer.add_argument('-l', '--lib', metavar='KMER_LIB', required=True,
                          help='Jellyfish kmer database')
    cmd_kmer.add_argument('-r', '--reads', metavar='FASTQ_FILE', action='append', required=True,
                          help='FastQ format reads used in making the kmer database '
                               '(Use multiple times for multiple files)')

    args = parser.parse_args()

    if args.version:
        print(mk_version())
        sys.exit(0)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    logging.captureWarnings(True)
    logger = logging.getLogger('main')

    # root log listens to everything
    root = logging.getLogger('')
    root.setLevel(logging.DEBUG)

    # log message format
    formatter = logging.Formatter(fmt='%(levelname)-8s | %(asctime)s | %(name)7s | %(message)s')

    # Runtime console listens to INFO by default
    ch = logging.StreamHandler()
    if args.verbose:
        ch.setLevel(logging.DEBUG)
    else:
        ch.setLevel(logging.INFO)
    ch.setFormatter(formatter)
    root.addHandler(ch)

    fh = logging.FileHandler(__log_name__, mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(mk_version())
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    try:

        if args.threads < 1:
            parser.error('The number of threads must be greater than 1')

        if args.sample_rate is not None and not (0 < args.sample_rate <= 1):
            parser.error('Sample rate must be within the range (0,1]')

        # BAM based analysis
        if args.command == 'bam':

            bam.analyze(args.bam, args.enzyme, args.mean_insert, seed=args.seed,
                        sample_rate=args.sample_rate, threads=args.threads)

        # Kmer based analysis
        elif args.command == 'kmer':

            if len(set(args.reads)) != len(args.reads):
                parser.error('Some supplied input read-sets are the same file')

            if len(args.enzyme) > 1:
                parser.error('kmer analysis mode only supports a single enyzme')

            kmer.analyze(args.kmer_size, args.enzyme[0], args.lib, args.reads, args.mean_insert,
                         sample_rate=args.sample_rate, seed=args.seed, max_coverage=args.max_coverage,
                         threads=args.threads, save_cov=args.save_cov)

    except Exception as ex:
        logger.exception(ex)
        sys.exit(1)
