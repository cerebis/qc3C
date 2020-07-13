import logging
import re
import subprocess
import tempfile
import psutil

from qc3C.exceptions import ApplicationException

logger = logging.getLogger(__name__)


def mk_database(db_path, fasta_files, kmer_size, hash_size, threads=1):
    """
    Create a Jellyfish kmer database from the supplied fasta files.

    :param db_path: output database path
    :param fasta_files: input fasta file list
    :param kmer_size: k-mer size
    :param hash_size: starting hash size used by jellyfish
    :param threads: number of concurrent threads
    """

    if re.fullmatch(r'[0-9]+[mMgG]', hash_size) is None:
        raise ValueError('invalid hash_size format supplied')

    with tempfile.NamedTemporaryFile(mode='wt') as gen_h:

        # Satisfy jellyfish's strange generator file input method for
        # multiple fasta file (R1/R2) support.
        for fn in fasta_files:
            gen_h.write('zcat -f {}\n'.format(fn))
        # make sure the buffer is flushed to disk
        gen_h.flush()

        try:
            logger.info('Beginning library creation')
            logger.info('Requested kmer size: {}'.format(kmer_size))
            logger.info('Input FastQ files: {}'.format(' '.join(fasta_files)))
            logger.info('Creating library: {}'.format(db_path))
            p = subprocess.Popen(['jellyfish', 'count',
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
