"""Tests for presence-absence matrix construction."""

from skipalign.matrix import build_pa_matrix, filter_conserved


def test_build_pa_matrix():
    kmer_sets = {
        "genome1": {"AAA", "BBB", "CCC"},
        "genome2": {"AAA", "BBB", "DDD"},
        "genome3": {"AAA", "EEE"},
    }
    matrix, kmer_index, genome_names = build_pa_matrix(kmer_sets)
    assert matrix.shape == (5, 3)  # 5 unique k-mers, 3 genomes
    assert genome_names == ["genome1", "genome2", "genome3"]


def test_filter_conserved():
    kmer_sets = {
        "g1": {"AAA", "BBB", "CCC"},
        "g2": {"AAA", "BBB"},
        "g3": {"AAA"},
    }
    matrix, kmer_index, _ = build_pa_matrix(kmer_sets)
    conserved, counts = filter_conserved(matrix, kmer_index, min_genomes=2)
    assert "AAA" in conserved
    assert "BBB" in conserved
    assert "CCC" not in conserved
