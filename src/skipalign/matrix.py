"""Presence-absence matrix construction from k-mer sets."""

from __future__ import annotations

import numpy as np
from scipy.sparse import csr_matrix


def build_pa_matrix(
    kmer_sets: dict[str, set[str]],
) -> tuple[csr_matrix, list[str], list[str]]:
    """Build a binary presence-absence matrix (k-mers × genomes).

    Returns:
        matrix: sparse binary matrix
        kmer_index: ordered list of k-mer strings (row labels)
        genome_names: ordered list of genome names (column labels)
    """
    genome_names = sorted(kmer_sets.keys())
    all_kmers = set()
    for kmers in kmer_sets.values():
        all_kmers.update(kmers)
    kmer_index = sorted(all_kmers)
    kmer_to_row = {km: i for i, km in enumerate(kmer_index)}

    rows, cols = [], []
    for col_idx, name in enumerate(genome_names):
        for km in kmer_sets[name]:
            rows.append(kmer_to_row[km])
            cols.append(col_idx)

    data = np.ones(len(rows), dtype=np.uint8)
    matrix = csr_matrix((data, (rows, cols)), shape=(len(kmer_index), len(genome_names)))
    return matrix, kmer_index, genome_names


def filter_conserved(
    matrix: csr_matrix,
    kmer_index: list[str],
    min_genomes: int,
) -> tuple[list[str], np.ndarray]:
    """Filter k-mers present in at least min_genomes genomes.

    Returns:
        conserved_kmers: list of conserved k-mer strings
        genome_counts: array of genome counts per conserved k-mer
    """
    counts = np.asarray(matrix.sum(axis=1)).flatten()
    mask = counts >= min_genomes
    conserved_kmers = [kmer_index[i] for i in range(len(kmer_index)) if mask[i]]
    genome_counts = counts[mask]
    return conserved_kmers, genome_counts
