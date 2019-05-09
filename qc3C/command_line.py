import logging
import qc3C.bam_based as bam
import qc3C.kmer_based as kmer

__version__ = '0.2'
__log_name__ = 'qc3C.log'


def mk_version():
    return 'qc3C v{}'.format(__version__)


def main():
    import argparse
    import sys

    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('--sep', default='\t', help='Delimiter to use in report table')
    global_parser.add_argument('-e', '--enzyme', metavar='NEB_NAME', required=True, action='append',
                               help='Case-sensitive NEB enzyme name. Use multiple times for multiple enzymes')

    parser = argparse.ArgumentParser(description='qc3C: Hi-C quality control')
    parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    parser.add_argument('-V', '--version', default=False, action='store_true', help='Version')
    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')
    subparsers.required = True
    cmd_bam = subparsers.add_parser('bam', parents=[global_parser],
                                      description='Alignment-based analysis.')
    cmd_kmer = subparsers.add_parser('kmer', parents=[global_parser],
                                        description='Kmer-based analysis.')
    """
    CLI for BAM based analysis
    """
    cmd_bam.add_argument('-t', '--threads', metavar='N', type=int, default=1, help='Number of threads')
    cmd_bam.add_argument('-m', '--mean-insert', type=int, required=True,
                          help='Mean fragment length to use in estimating the unobserved junction rate')
    cmd_bam.add_argument('BAM', help='Input name-sorted bam file of Hi-C reads mapped to references')

    """
    CLI for Kmer based analysis
    """
    cmd_kmer.add_argument('-s', '--seed', type=int,
                          help='Random seed used in sampling the read-set')
    cmd_kmer.add_argument('-n', '--max-reads', default=500000, type=int,
                          help='Stop after collecting N sample reads')
    cmd_kmer.add_argument('-a', '--accept-all', default=False, action='store_true',
                          help='Override acceptance rate and accept all useable reads')
    cmd_kmer.add_argument('-m', '--mean-insert', type=int,
                          help='Mean fragment length to use in estimating the unobserved junction rate')
    cmd_kmer.add_argument('-x', '--max-coverage', default=500, type=int,
                          help='Ignore regions with more than this coverage')
    cmd_kmer.add_argument('-N', '--pool-size', type=int,
                          help='The total number of reads which are being provided for consideration')
    cmd_kmer.add_argument('KMER_SIZE', type=int, help='Kmer size used in database')
    cmd_kmer.add_argument('FASTQ', help='FastQ file used in making the kmer database')
    cmd_kmer.add_argument('KMER_DB', help='Jellyfish kmer database')

    args = parser.parse_args()

    if args.version:
        print(mk_version())
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

        # BAM based analysis
        if args.command == 'bam':
            bam.analyze(args.BAM, args.enzyme, args.mean_insert, threads=args.threads, sep=args.sep)

        # Kmer based analysis
        elif args.command == 'kmer':
            assert len(args.enzyme) == 1, 'Kmer-based approach currently supports only a single enzyme'
            kmer.analyze(args.KMER_SIZE, args.enzyme[0], args.KMER_DB, args.FASTQ, args.mean_insert,
                         pool_size=args.pool_size, max_reads=args.max_reads, seed=args.seed,
                         max_coverage=args.max_coverage, accept_all=args.accept_all)

    except Exception as ex:
        logger.error(ex)
        sys.exit(1)
