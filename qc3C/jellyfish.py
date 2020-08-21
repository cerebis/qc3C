import logging
import re
import subprocess
import tempfile
import psutil

from typing import List
from qc3C.exceptions import ApplicationException
from qc3C.utils import guess_quality_encoding

logger = logging.getLogger(__name__)


def mk_database(db_path: str, fastq_files: List[str], kmer_size: int, hash_size: int,
                ascii_base: int = None, min_qual: int = 5, threads: int = 1):
    """
    Create a Jellyfish kmer database from the supplied fasta files.

    :param db_path: output database path
    :param fastq_files: input FastQ file list
    :param kmer_size: k-mer size
    :param ascii_base: Ascii-base for quality score encoding  (33 new encoding, 64 old encoding)
    :param min_qual: Minimum acceptable base quality or position converted to N
    :param hash_size: starting hash size used by jellyfish
    :param threads: number of concurrent threads
    """

    if min_qual is None:
        min_qual = 0

    if re.fullmatch(r'[0-9]+[mMgG]', hash_size) is None:
        raise ValueError('invalid hash_size format supplied')

    # guess the quality encoding if not specified
    if ascii_base is None:
        ascii_base = guess_quality_encoding(fastq_files[0])
        logger.info('Inferred FastQ quality encoding from sampled reads: {}'.format(ascii_base))
    else:
        logger.info('User specified FastQ quality encoding is: {}'.format(ascii_base))

    with tempfile.NamedTemporaryFile(mode='wt') as gen_h:

        # Satisfy jellyfish's strange generator file input method for
        # multiple fasta file (R1/R2) support.
        for fn in fastq_files:
            gen_h.write('zcat -f {}\n'.format(fn))
        # make sure the buffer is flushed to disk
        gen_h.flush()

        try:
            logger.info('Beginning library creation')
            logger.info('Requested minimum quality: {}'.format(min_qual))
            logger.info('Requested kmer size: {}'.format(kmer_size))
            logger.info('Input FastQ files: {}'.format(' '.join(fastq_files)))
            logger.info('Creating library: {}'.format(db_path))
            p = subprocess.Popen(['jellyfish', 'count',
                                  '--quality-start', str(ascii_base),
                                  '--min-quality', str(min_qual),
                                  '-C',
                                  '-t', str(threads),
                                  '-G', str(threads),
                                  '-m', str(kmer_size),
                                  '-s', hash_size,
                                  '-g', gen_h.name,
                                  '-o', db_path],
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)

            return_value = None
            while True:
                try:
                    return_value = p.poll()
                    if return_value is not None:
                        break
                    out, err = p.communicate(timeout=3)
                except subprocess.TimeoutExpired:
                    proc = psutil.Process(p.pid)
                    if proc.status() == psutil.STATUS_ZOMBIE:
                        logger.error('Jellyfish aborted, orphaned subprocesses will require termination')
                        break
                    elif not proc.is_running():
                        # jellyfish ended normally
                        break

            if return_value is None or return_value > 0 :
                logger.warning('Jellyfish subprocess did not return 0 (ok). Return value was: {}'.format(return_value))
                if err:
                    logger.warning('Jellyfish stderr: {}'.format(err.decode()))
            else:
                if out:
                    logger.debug('Jellyfish stdout: {}'.format(out.decode()))

        except subprocess.CalledProcessError as e:
            logger.error('An exception occurred during kmer database creation')
            raise ApplicationException(e.stdout)
