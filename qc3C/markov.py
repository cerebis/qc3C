import itertools
import numpy as np
import tqdm

from collections import Counter as collectionsCounter
from typing import List, Dict, Optional


def choice(options, probs, rs):
    """
    A faster method for non-uniform selection, when compared
    to numpy.random.choice.

    NOTE This could be further doubled in speed, if numba is used,
    however an initialised random state in numba is only maintained
    within the invoking code block and therefore a larger method would
    have to encapsulate this method.

    :param options: the list of options from which to select
    :param probs: the probability of each choice
    :param rs: numpy random state
    :return: a randomly chosen option
    """
    x = rs.rand()
    cum = 0
    i = 0
    for i, p in enumerate(probs):
        cum += p
        if x < cum:
            break
    return options[i]


class ReadMarkovModeller:
    """
    Use a k-th order Markov model to compute the expected number of junction
    sequences that would occur in a set of sequences with the given composition
    """

    def __init__(self, kmer_size: int, seed: int, reject_prob: float, targets: List[str]):
        """
        :param kmer_size:
        :param seed:
        :param reject_prob:
        :param targets:
        """
        self.kmer_size = kmer_size
        self.counts = collectionsCounter()
        self.counts_plus_one = collectionsCounter()
        self.trans_prob = None
        self.kmers_prob = None
        self.kmers = None
        self.reject_prob = reject_prob
        self.targets = [ti.encode() for ti in targets]
        self.random_state = np.random.RandomState(seed)
        self.total_obs = 0
        self.total_extent = 0

    def average_read_length(self) -> int:
        """
        :return: Average read length accumulated by the model instance
        """
        return self.total_extent // self.total_obs

    @staticmethod
    def filter_n(kmer_dict: collectionsCounter) -> Dict[bytes, int]:
        """
        Remove records whose k-mer contains one or more Ns
        """
        return {k: v for k, v in kmer_dict.items() if b'N' not in k}

    def contains_target(self, s):
        return any(ti in s for ti in self.targets)

    def prepare_model(self):
        """
        Once some training data has been added, this method will prepare the transition
        probabilty data structure for use in simulating reads.
        """
        # filter k-mers containing ambiguous bases
        clean_counts = ReadMarkovModeller.filter_n(self.counts)

        # rather than rely on kmers from counts, we will assume there
        # are missing observations and instead explicitly iterate
        # over all possibilities.
        self.kmers = [bytes(ki) for ki in itertools.product(b'ACGT', repeat=self.kmer_size)]
        # assign probabilities for each kmer
        self.kmers_prob = [clean_counts[kmer] if kmer in clean_counts else 1 for kmer in self.kmers]
        # normalise these counts to probs
        self.kmers_prob = np.array(self.kmers_prob, dtype=np.float)
        self.kmers_prob /= self.kmers_prob.sum()

        uniform_prob = np.array([0.25] * 4)
        self.trans_prob = {}
        for kmer in self.kmers:
            if kmer not in clean_counts:
                # if we didn't see the kmer, assume uniform probability
                self.trans_prob[kmer] = uniform_prob
            else:
                # generate the four posssible k+1-mers
                plus_ones = [kmer + nt for nt in [b'A', b'C', b'G', b'T']]
                # get the counts for these possibilities
                next_counts = np.array([self.counts_plus_one[po] if po else 0 for po in plus_ones])
                # transition probabilities for this k-mer to the next 4 possibilities
                self.trans_prob[kmer] = next_counts / next_counts.sum()

    def find_target_mid(self, seq: bytes) -> Optional[int]:
        """
        Find target sequences within a read and return the modpoint of the
        matched region. Only the first match is found.

        :param seq: the sequence to search
        """
        self.random_state.shuffle(self.targets)
        for ti in self.targets:
            ix = seq.find(ti)
            if ix != -1:
                return ix + len(ti) // 2
        return None

    def __iadd__(self, other: str):
        """
        Introduce another sequence for training model. In the case of reads which contain
        a target sequence, split the reads at the target sequence midpoint prior to extracting
        k-mers from each end. This aims to suppress k-mer bias in high quality libraries which
        contain a high proportion of target sequences.

        :param other: (sequence is assumed to already be upper case)
        """
        other = other.upper().encode()

        # split reads containing junction at midpoint at given rate
        # this surpresses Hi-C signal bias which is a problem for high quality libraries
        if self.contains_target(other) and self.random_state.uniform() < self.reject_prob:
            break_point = self.find_target_mid(other)
            for oi in [other[:break_point], other[break_point:]]:
                # extract k and k+1 mers
                self.counts.update(oi[i: i+self.kmer_size] for i in range(len(oi) - self.kmer_size + 1))
                self.counts_plus_one.update(oi[i: i+self.kmer_size+1] for i in range(len(oi) - self.kmer_size))
        else:
            # extract k and k+1 mers
            self.counts.update(other[i: i+self.kmer_size] for i in range(len(other) - self.kmer_size + 1))
            self.counts_plus_one.update(other[i: i+self.kmer_size+1] for i in range(len(other) - self.kmer_size))

        self.total_obs += 1
        self.total_extent += len(other)
        return self

    def estimate_freq(self, sample_size: int) -> float:
        """
        Compute the expected frequency at which any of the target sequences appear in reads
        simulated from the model.

        :param sample_size:
        :return:
        """

        read_len = int(self.average_read_length())

        # simulate sequences at the target read length and count the fraction
        # that contain the junction sequence(s)
        acgt = bytearray(b'ACGT')
        # convert the list of junctions to a set of bytes
        containing_reads = 0
        # generate an initial starting point for each read
        initial_kmers = self.random_state.choice(self.kmers, size=sample_size, p=self.kmers_prob)
        for kmer in tqdm.tqdm(initial_kmers):
            # choose a kmer based on observed prob
            # start the simulated sequence from this
            sim_seq = bytearray(kmer + (b'N' * (read_len - self.kmer_size)))
            for i in range(self.kmer_size, len(sim_seq)):
                sim_seq[i] = choice(acgt, self.trans_prob[bytes(sim_seq[i-self.kmer_size: i])], self.random_state)
            # does our simulated sequence contain the junction?
            if self.contains_target(sim_seq):
                containing_reads += 1

        return containing_reads / sample_size