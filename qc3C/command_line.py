import logging
import qc3C.bam_based as bam
import qc3C.kmer_based as kmer
from qc3C.exceptions import ApplicationException
from qc3C._version import version_stamp

__log_name__ = 'qc3C.log'


def main():
    import argparse
    import sys

    class UniqueStore(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if getattr(namespace, self.dest, self.default) is not None:
                parser.error('duplicate use of the option: {}'.format(option_string))
            setattr(namespace, self.dest, values)

    """
    Shared CLI arguments
    """
    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('-p', '--sample-rate', default=None, type=float, action=UniqueStore,
                               help='Sample only a proportion of all read-pairs [None]')
    global_parser.add_argument('-s', '--seed', type=int, action=UniqueStore,
                               help='Random seed used in sampling the read-set [None]')
    global_parser.add_argument('-t', '--threads', metavar='N', type=int, default=1, help='Number of threads [1]')
    # global_parser.add_argument('-n', '--num-obs', type=int, action=UniqueStore,
    #                            help='Stop parsing after collecting N observations [None]')
    global_parser.add_argument('-m', '--mean-insert', type=int, required=True, action=UniqueStore,
                               help='Mean fragment length to use in estimating the unobserved junction rate')
    global_parser.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append', required=True,
                               help='One or more case-sensitive NEB enzyme names'
                                    '(Use multiple times for multiple files enzymes)')

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
    cmd_bam.add_argument('-q', '--min-mapq', default=60, type=int, action=UniqueStore,
                         help='Minimum acceptable mapping quality [60]')
    cmd_bam.add_argument('-b', '--bam', required=True, action=UniqueStore,
                         help='Input name-sorted bam file of Hi-C reads mapped to references')

    """
    CLI for Kmer based analysis
    """
    cmd_kmer.add_argument('--output-table', default=None, action=UniqueStore,
                          help='Save the collected per-read statistics table to a file')
    cmd_kmer.add_argument('-x', '--max-coverage', default=500, type=int, action=UniqueStore,
                          help='Ignore regions with more than this coverage [500]')
    cmd_kmer.add_argument('-l', '--lib', metavar='KMER_LIB', required=True, action=UniqueStore,
                          help='Jellyfish kmer database')
    cmd_kmer.add_argument('-r', '--reads', metavar='FASTQ_FILE', action='append', required=True,
                          help='FastQ format reads used in making the kmer database '
                               '(Use multiple times for multiple files)')

    args = parser.parse_args()

    if args.version:
        print(version_stamp())
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
    logger.debug(version_stamp(False))
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
                        sample_rate=args.sample_rate, threads=args.threads,
                        min_mapq=args.min_mapq, num_obs=None)

        # Kmer based analysis
        elif args.command == 'kmer':

            if len(set(args.reads)) != len(args.reads):
                parser.error('Some supplied input read-sets are the same file')

            if len(args.enzyme) > 1:
                parser.error('kmer analysis mode only supports a single enyzme')

            kmer.analyze(args.enzyme[0], args.lib, args.reads, args.mean_insert,
                         sample_rate=args.sample_rate, seed=args.seed, max_coverage=args.max_coverage,
                         threads=args.threads, output_table=args.output_table, num_obs=None)

    except ApplicationException as ex:
        logger.error(str(ex))
        sys.exit(1)

    except Exception as ex:
        logger.exception(ex)
        sys.exit(1)
