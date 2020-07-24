from qc3C.utils import *
import pytest


def test_simple_observed_fraction():
    # Setup
    obs_frac = 0.5
    mean_frag_size = 300
    n_fragments = 1000
    obs_extent = 150000
    # Exercise
    result = simple_observed_fraction(obs_extent, mean_frag_size, n_fragments)
    # Verify
    assert result == obs_frac
    # Cleanup - not required


@pytest.mark.parametrize('exc_class,file_path,fmt,ncpu,wait_timeout',
                         [(FileNotFoundError, 'foobar', 'fasta', 1, None),
                          (ApplicationException, 'data/10seq.fa', 'foobar', 1, None),
                          (ApplicationException, 'data/10seq.fa', 'fasta', -1, None),
                          (ApplicationException, 'data/10seq.fa', 'fasta', 1, 0)])
def test_count_sequences_exceptions(exc_class, file_path, fmt, ncpu, wait_timeout):
    # Exercise and verify
    with pytest.raises(exc_class):
        count_sequences(file_path, fmt, ncpu, wait_timeout)


@pytest.mark.parametrize("nseq,file_path,fmt,ncpu",
                         [(10, 'data/10seq.fa', 'fasta', 1),
                          (10, 'data/10seq.fq', 'fastq', 1),
                          (10, 'data/10seq.fa.gz', 'fasta', 1),
                          (10, 'data/10seq.fq.gz', 'fastq', 1),
                          (10, 'data/10seq.fa', 'fasta', 2),
                          (10, 'data/10seq.fa.gz', 'fasta', 2)])
def test_count_sequences(nseq, file_path, fmt, ncpu):
    # Exercise and verify
    assert nseq == count_sequences(file_path, fmt, ncpu)
