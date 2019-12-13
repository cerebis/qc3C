import os
import logging
import qc3C.bam_based as bam
import qc3C.kmer_based as kmer

from qc3C.exceptions import ApplicationException
from qc3C._version import version_stamp


__log_name__ = 'qc3C.log'
__report_name__ = 'report.jsonl'
__table_name__ = 'obs_table.tsv.gz'


def main():
    import argparse
    import sys

    class UniqueStore(argparse.Action):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.count = 0

        def __call__(self, parser, namespace, values, option_string=None):
            if self.count >= 1:
                parser.error('duplicate use of the option: {}'.format(option_string))
            setattr(namespace, self.dest, values)
            self.count += 1

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
    global_parser.add_argument('--output-path', metavar='PATH', default='.',
                               help='Write output files to this folder [.]')
    global_parser.add_argument('--write-report', default=False, action='store_true',
                               help='Create a result report in JSONLines format')
    global_parser.add_argument('-k', '--library-kit', choices=['phase', 'generic'], default='generic',
                               help='The library kit type [generic]')
    global_parser.add_argument('-m', '--mean-insert', type=int, required=True, action=UniqueStore,
                               help='Mean fragment length to use in estimating the unobserved junction rate')
    global_parser.add_argument('-e', '--enzyme', metavar='NEB_NAME', action=UniqueStore, required=True,
                               help='A case-sensitive NEB enzyme name')

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
    cmd_bam.add_argument('-f', '--fasta', required=True, action=UniqueStore,
                         help='Reference sequences')
    cmd_bam.add_argument('-b', '--bam', required=True, action=UniqueStore,
                         help='Input name-sorted bam file of Hi-C reads mapped to references')

    """
    CLI for Kmer based analysis
    """
    cmd_kmer.add_argument('--write-table', default=False, action='store_true',
                          help='Save the collected observations to a file')
    cmd_kmer.add_argument('-x', '--max-freq-quantile', default=0.9, type=float, action=UniqueStore,
                          help='Ignore regions possessing k-mer frequencies above this quantile [0.9]')
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

    # create the output folder if it did not exist
    if not os.path.exists(args.output_path):
        os.mkdir(args.output_path)

    fh = logging.FileHandler(os.path.join(args.output_path, __log_name__), mode='a')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    root.addHandler(fh)

    # Add some environmental details
    logger.debug(version_stamp(False))
    logger.debug(sys.version.replace('\n', ' '))
    logger.debug('Command line: {}'.format(' '.join(sys.argv)))

    report_path = None
    if args.write_report:
        report_path = os.path.join(args.output_path, __report_name__)

    try:

        if args.threads < 1:
            parser.error('The number of threads must be greater than 1')

        if args.sample_rate is not None and not (0 < args.sample_rate <= 1):
            parser.error('Sample rate must be within the range (0,1]')

        # BAM based analysis
        if args.command == 'bam':

            bam.analyse(args.bam, args.fasta, args.enzyme, args.mean_insert,
                        seed=args.seed, sample_rate=args.sample_rate, threads=args.threads,
                        min_mapq=args.min_mapq, report_path=report_path,
                        library_kit=args.library_kit)

        # Kmer based analysis
        elif args.command == 'kmer':

            table_path = None
            if args.write_table:
                table_path = os.path.join(args.output_path, __table_name__)

            if len(set(args.reads)) != len(args.reads):
                parser.error('Some supplied input read-sets are the same file')

            kmer.analyse(args.enzyme, args.lib, args.reads, args.mean_insert,
                         sample_rate=args.sample_rate, seed=args.seed, max_freq_quantile=args.max_freq_quantile,
                         threads=args.threads, output_table=table_path, report_path=report_path,
                         library_kit=args.library_kit)

    except ApplicationException as ex:
        logger.error(str(ex))
        sys.exit(1)

    except Exception as ex:
        logger.exception(ex)
        sys.exit(1)
