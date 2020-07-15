import os
import logging
import pathlib
import qc3C.bam_based as bam
import qc3C.kmer_based as kmer
from qc3C.jellyfish import mk_database

from qc3C.exceptions import ApplicationException
from qc3C._version import version_stamp


__log_name__ = 'qc3C.log'
__report_name__ = 'report.qc3C'
__table_name__ = 'obs_table.tsv.gz'


def main():
    import argparse
    import sys

    class UniqueStore(argparse.Action):
        """
        Restrict an argparse argument to only one occurrence
        """
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.count = 0

        def __call__(self, parser, namespace, values, option_string=None):
            if self.count >= 1:
                parser.error('duplicate use of the option: {}'.format(option_string))
            setattr(namespace, self.dest, values)
            self.count += 1

    class Range(object):
        """
        Define a range of acceptable values for argparse argument
        """
        def __init__(self, start, end, include_start: bool = True, include_end: bool = True):
            self.start = start
            self.end = end
            self.include_start = include_start
            self.include_end = include_end

        def __eq__(self, other):
            if not self.include_start and not self.include_end:
                return self.start < other < self.end
            elif not self.include_start:
                return self.start < other <= self.end
            elif not self.include_end:
                return self.start <= other < self.end
            else:
                return self.start <= other <= self.end

        def __contains__(self, item):
            return self.__eq__(item)

        def __iter__(self):
            yield self

        def __repr__(self):
            if not self.include_start and not self.include_end:
                return 'range ({0},{1})'.format(self.start, self.end)
            elif not self.include_start:
                return 'range ({0},{1}]'.format(self.start, self.end)
            elif not self.include_end:
                return 'range [{0},{1})'.format(self.start, self.end)
            else:
                return 'range [{0},{1}]'.format(self.start, self.end)

    """
    Shared CLI arguments
    """
    global_parser = argparse.ArgumentParser(add_help=False)
    global_parser.add_argument('-d', '--debug', default=False, action='store_true', help='Enable debug output')
    global_parser.add_argument('-v', '--verbose', default=False, action='store_true', help='Verbose output')
    global_parser.add_argument('-t', '--threads', metavar='N', type=int, default=1, help='Number of threads [1]')
    global_parser.add_argument('-o', '--output-path', metavar='PATH', default='.', action=UniqueStore,
                               help='Write output files to this folder [.]')

    analysis_parser = argparse.ArgumentParser(add_help=False)
    analysis_parser.add_argument('-p', '--sample-rate', default=None, type=float,
                                 choices=Range(0, 1, include_start=False), action=UniqueStore,
                                 help='Sample only a proportion of all read-pairs [None]')
    analysis_parser.add_argument('-s', '--seed', type=int, action=UniqueStore,
                                 help='Random seed used in sampling the read-set [None]')
    analysis_parser.add_argument('-M', '--max-obs', type=int, action=UniqueStore,
                                 help='Stop after collecting this many observations')
    analysis_parser.add_argument('--no-json', default=False, action='store_true',
                                 help='Do not write a JSON report')
    analysis_parser.add_argument('--no-html', default=False, action='store_true',
                                 help='Do not write an HTML report')

    digestion_group = analysis_parser.add_mutually_exclusive_group(required=True)
    digestion_group.add_argument('-k', '--library-kit', choices=['phase', 'arima'], default=None, action=UniqueStore,
                                 help='Define digest by the commercial library kit used')
    digestion_group.add_argument('-e', '--enzyme', metavar='NEB_NAME', action='append',
                                 help='Define digest by explicitly naming up to two case-sensitive NEB enzyme names')

    parser = argparse.ArgumentParser(description='qc3C: Hi-C quality control')
    parser.add_argument('-V', '--version', default=False, action='store_true', help='Version')

    subparsers = parser.add_subparsers(title='commands', dest='command', description='Valid commands',
                                       help='choose an analysis stage for further options')
    subparsers.required = False
    cmd_bam = subparsers.add_parser('bam', parents=[global_parser, analysis_parser],
                                    description='Alignment-based analysis.')
    cmd_kmer = subparsers.add_parser('kmer', parents=[global_parser, analysis_parser],
                                     description='Kmer-based analysis.')
    cmd_mkdb = subparsers.add_parser('mkdb', parents=[global_parser],
                                     description='Create kmer database.')

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
    cmd_kmer.add_argument('--merged-reads', default=False, action='store_true',
                          help='Input reads are merged pairs')
    cmd_kmer.add_argument('--write-table', default=False, action='store_true',
                          help='Save the collected observations to a file')
    cmd_kmer.add_argument('--num-sample', type=int, default=50, action=UniqueStore,
                          help='Number of samples to use in bootstrapping confidence interval [50]')
    cmd_kmer.add_argument('--frac-sample', type=float, default=0.33,
                          choices=Range(0, 1, include_start=False), action=UniqueStore,
                          help='Fraction of observations to use per-bootstrap iteration [1/3]')
    cmd_kmer.add_argument('-x', '--max-freq-quantile', default=0.9, type=float,
                          choices=Range(0, 1, include_start=False), action=UniqueStore,
                          help='Ignore k-mers possessing frequencies above this quantile [0.9]')
    cmd_kmer.add_argument('-m', '--mean-insert', type=int, required=True, action=UniqueStore,
                          help='Mean fragment length to use in estimating the unobserved junction rate')
    cmd_kmer.add_argument('-l', '--lib', metavar='KMER_LIB', required=True, action=UniqueStore,
                          help='Jellyfish kmer database')
    cmd_kmer.add_argument('-r', '--reads', metavar='FASTQ_FILE', action='append', required=True,
                          help='FastQ format reads that were used to generate the k-mer library '
                               '(use multiple times for multiple files)')

    """
    CLI for creating kmer database
    """
    cmd_mkdb.add_argument('-s', '--hash-size', default='10M', action=UniqueStore,
                          help='Initial hash size (eg. 10M, 2G) [10M]')
    cmd_mkdb.add_argument('-k', '--kmer-size', default=24, type=int, action=UniqueStore,
                          help='Library kmer size [24]')
    cmd_mkdb.add_argument('-r', '--reads', metavar='FASTQ_FILE', action='append', required=True,
                          help='FastQ format reads to use in generating the k-mer library '
                               '(use multiple times for multiple files)')
    cmd_mkdb.add_argument('-l', '--lib', metavar='KMER_LIB', required=True,
                          help='Output Jellyfish k-mer library base name')

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
    logging.getLogger('numba').setLevel(logging.INFO)

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
        pathlib.Path(args.output_path).mkdir(parents=True, exist_ok=True)

    fh = logging.FileHandler(os.path.join(args.output_path, __log_name__), mode='a')
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

        report_path = os.path.join(args.output_path, __report_name__)

        # check if the user has employed the library kit option to declare enzymes
        if args.command in {'kmer', 'bam'} and args.library_kit is not None:
            # commercial kit definitions
            kit_choices = {'phase': ['Sau3AI', 'MluCI'],
                           'arima': ['DpnII', 'HinfI']}
            args.enzyme = kit_choices[args.library_kit]
            logger.info('Library kit {} declares enzymes {}'.format(args.library_kit, args.enzyme))

        if args.command == 'mkdb':

            if os.path.basename(args.lib) != args.lib:
                logger.error('The supplied library name should not include any path')
                sys.exit(1)

            if not args.lib.endswith('.jf'):
                args.lib = '{}.jf'.format(args.lib)

            klib_path = os.path.join(args.output_path, args.lib)

            if os.path.exists(klib_path):
                logger.error('Output path already exists: {}'.format(klib_path))
                sys.exit(1)

            mk_database(klib_path, args.reads, args.kmer_size, args.hash_size, args.threads)

        # BAM based analysis
        elif args.command == 'bam':

            bam.analyse(args.enzyme, args.bam, args.fasta,
                        seed=args.seed, sample_rate=args.sample_rate, threads=args.threads,
                        min_mapq=args.min_mapq, max_obs=args.max_obs, report_path=report_path,
                        no_json=args.no_json, no_html=args.no_html)

        # Kmer based analysis
        elif args.command == 'kmer':

            table_path = None if not args.write_table else os.path.join(args.output_path, __table_name__)

            if len(set(args.reads)) != len(args.reads):
                parser.error('Some supplied input read-sets are the same file')

            if args.merged_reads:
                logger.info('Merged reads specified, insert size will be ignored')

            kmer.analyse(args.enzyme, args.lib, args.reads, args.mean_insert,
                         sample_rate=args.sample_rate, seed=args.seed, max_freq_quantile=args.max_freq_quantile,
                         threads=args.threads, output_table=table_path, report_path=report_path,
                         no_json=args.no_json, no_html=args.no_html, max_obs=args.max_obs,
                         num_sample=args.num_sample, frac_sample=args.frac_sample)

    except ApplicationException as ex:
        logger.error(str(ex))
        if args.debug:
            logger.exception(ex)
        sys.exit(1)

    except Exception as ex:
        logger.error(str(ex))
        if args.debug:
            logger.exception(ex)
        sys.exit(1)
