from qc3C.utils import *
import pytest


@pytest.mark.parametrize('obs_frac, mean_frag_size',
                         [(0.25, 600),
                          (0.5, 300),
                          (1.0, 150)])
def test_simple_observed_fraction(obs_frac, mean_frag_size):
    # Setup
    n_fragments = 1000
    obs_extent = 150000
    # Exercise and verify
    assert simple_observed_fraction(obs_extent, mean_frag_size, n_fragments) == obs_frac


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
    assert count_sequences(file_path, fmt, ncpu) == nseq


@pytest.mark.parametrize('obs_mask, read_len, insert_len, kmer_size, junc_size',
                         [([1] * 20, 10, 20, 0, 0),
                          ([0] + [1] * 8 + [0, 0] + [1] * 8 + [0], 10, 20, 1, 0),
                          ([0] * 2 + [1] * 6 + [0] * 4 + [1] * 6 + [0] * 2, 10, 20, 2, 0),
                          ([1] * 9 + [0, 0] + [1] * 9, 10, 20, 0, 1),
                          ([0] + [1] * 7 + [0] * 4 + [1] * 7 + [0], 10, 20, 1, 1)])
def test_make_observable_mask(obs_mask, read_len, insert_len, kmer_size, junc_size):
    # setup
    obs_mask = np.array(obs_mask, dtype=np.float)
    # exercise and verify
    result = make_observable_mask(read_len, insert_len, kmer_size, junc_size)
    assert np.all(result == obs_mask)


@pytest.mark.parametrize('obs_frac, read_len, insert_len, method, kmer_size, junc_size',
                         [(1.0, 10, 20, 'binary', 0, 0),
                          (0.5, 5, 20, 'binary', 0, 0),
                          (1.0, 15, 20, 'binary', 0, 0),
                          (0.9, 15, 20, 'binary', 1, 1),
                          (1.5, 15, 20, 'additive', 0, 0),
                          (1.0, 10, 20, 'additive', 0, 0),
                          (1.2, 15, 20, 'additive', 1, 1)])
def test_observed_fraction(obs_frac, read_len, insert_len, method, kmer_size, junc_size):
    assert observed_fraction(read_len, insert_len, method, kmer_size, junc_size) == obs_frac


@pytest.mark.parametrize('exc_class, method', [(ApplicationException, 'foobar')])
def test_observed_fraction_exceptions(exc_class, method):
    with pytest.raises(exc_class):
        observed_fraction(10, 20, method, 0, 0)
