"""Tests for strand orientation inference and normalization."""

from pathlib import Path

from Bio.Seq import Seq

from skipalign.io import load_genomes
from skipalign.kmer import count_kmers
from skipalign.mapper import find_exact_matches
from skipalign.matrix import build_pa_matrix, filter_conserved
from skipalign.strand import infer_orientations, normalize_genomes
from skipalign.unitig import extract_unitigs

GENOMES_DIR = Path(__file__).parent / "data" / "genomes"


def test_rc_genomes_detected():
    """RC genomes should be detected and marked with orientation '-'."""
    genomes = load_genomes(GENOMES_DIR)
    k = 19
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    matrix, kmer_index, _ = build_pa_matrix(kmer_sets)
    conserved, _ = filter_conserved(matrix, kmer_index, min_genomes=3)
    unitigs = extract_unitigs(conserved)
    hits = find_exact_matches(unitigs, genomes)

    ref, calls = infer_orientations(genomes, hits)

    # RC genomes should be detected
    rc_genomes = [g for g, c in calls.items() if c.orientation == "-"]
    assert "TestVirus_G_RC" in rc_genomes
    assert "TestVirus_H_RC" in rc_genomes

    # Forward genomes should stay "+"
    for name in ["TestVirus_A", "TestVirus_B", "TestVirus_C"]:
        assert calls[name].orientation == "+"


def test_normalization_produces_same_orientation():
    """After normalization, RC genomes should match their forward originals."""
    genomes = load_genomes(GENOMES_DIR)
    k = 19
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    matrix, kmer_index, _ = build_pa_matrix(kmer_sets)
    conserved, _ = filter_conserved(matrix, kmer_index, min_genomes=3)
    unitigs = extract_unitigs(conserved)
    hits = find_exact_matches(unitigs, genomes)

    _, calls = infer_orientations(genomes, hits)
    normalized = normalize_genomes(genomes, calls)

    # TestVirus_G_RC normalized should equal TestVirus_G original
    assert normalized["TestVirus_G_RC"] == genomes["TestVirus_G"]
    assert normalized["TestVirus_H_RC"] == genomes["TestVirus_H"]

    # Forward genomes should be unchanged
    assert normalized["TestVirus_A"] == genomes["TestVirus_A"]
