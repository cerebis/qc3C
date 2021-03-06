#!/usr/bin/env python3

import os
import argparse
import subprocess
import tempfile

# from qc3C.jellyfish import mk_database

if __name__ == '__main__':

    parser = argparse.ArgumentParser('Make a Jellyfish kmer database')
    parser.add_argument('-t', '--threads', default=1, type=int,
                        help='Number of threads [1]')
    parser.add_argument('-k', '--kmer-size', default=24, type=int,
                        help='Library kmer size [24]')
    parser.add_argument('-o', '--output', required=True, 
                        help='Output database name')
    parser.add_argument('FASTQ', nargs='+', help='FastQ read file')
    args = parser.parse_args()

    if not args.output.endswith('.jf'):
        args.output = '{}.jf'.format(args.output)

    with tempfile.NamedTemporaryFile(mode='wt') as gen_h:
        for fn in args.FASTQ:
            gen_h.write('zcat -f {}\n'.format(fn))

        out_dir = os.path.dirname(os.path.realpath(args.output))
        if not out_dir:
            raise IOError('Failed to determine containing output directory')

        print('Beginning library creation')
        with open(os.path.join(out_dir, 'mk_kmer_db.log'), 'w+') as stdout:
            stdout.write('Beginning library creation\n')
            stdout.write('Requested kmer size is: {}\n'.format(args.kmer_size))
            stdout.write('Input FastQ files: {}\n'.format(' '.join(args.FASTQ)))
            stdout.write('Output library file: {}\n'.format(args.output))
            subprocess.check_call(['jellyfish', 'count', '-t', str(args.threads), '-m', str(args.kmer_size),
                                   '-s', '100M', '-C', '-o', args.output, '-g', gen_h.name],
                                  stdout=stdout, stderr=subprocess.STDOUT)
            stdout.write('Finished\n')
        print('Finished')

