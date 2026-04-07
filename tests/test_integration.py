"""End-to-end integration test using synthetic test genomes.

The test data has 8 genomes with a planted conserved region at positions
1000-1300 (300bp, ~97% identity). The pipeline should discover this region.
"""

from pathlib import Path

from skipalign.io import load_genomes
from skipalign.kmer import count_kmers
from skipalign.matrix import build_pa_matrix, filter_conserved
from skipalign.unitig import extract_unitigs
from skipalign.mapper import find_exact_matches, parse_gff3, annotate_hits
from skipalign.scorer import score_windows

TEST_DATA = Path(__file__).parent / "data"
GENOMES_DIR = TEST_DATA / "genomes"
ANNOTATIONS_DIR = TEST_DATA / "annotations"

# Known planted conserved region
CONSERVED_START = 1000
CONSERVED_END = 1300
TOLERANCE = 50  # allow some positional tolerance


def test_load_genomes():
    genomes = load_genomes(GENOMES_DIR)
    assert len(genomes) == 8
    assert all(len(seq) == 2000 for seq in genomes.values())


def test_kmer_counting():
    genomes = load_genomes(GENOMES_DIR)
    k = 19
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    assert len(kmer_sets) == 8
    # Each 2000bp genome should have ~1982 possible 19-mers (minus non-ACGT)
    for name, kmers in kmer_sets.items():
        assert len(kmers) > 1000, f"{name} has too few k-mers: {len(kmers)}"


def test_pa_matrix_and_conservation():
    genomes = load_genomes(GENOMES_DIR)
    k = 19
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    matrix, kmer_index, genome_names = build_pa_matrix(kmer_sets)

    assert matrix.shape[1] == 8  # 8 genomes
    assert matrix.shape[0] > 1000  # many unique k-mers

    # Filter conserved k-mers (present in ≥6 of 8 genomes = 75%)
    conserved, counts = filter_conserved(matrix, kmer_index, min_genomes=6)
    # The conserved region should produce shared k-mers
    assert len(conserved) > 0, "No conserved k-mers found!"
    assert len(conserved) < len(kmer_index), "Too many k-mers are conserved"


def test_unitig_extraction():
    genomes = load_genomes(GENOMES_DIR)
    k = 19
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    matrix, kmer_index, _ = build_pa_matrix(kmer_sets)
    conserved, _ = filter_conserved(matrix, kmer_index, min_genomes=6)

    unitigs = extract_unitigs(conserved)
    assert len(unitigs) > 0, "No unitigs extracted"
    # Unitigs should be at least k-mer length
    assert all(len(u) >= k for u in unitigs)


def test_mapping_finds_conserved_region():
    genomes = load_genomes(GENOMES_DIR)
    k = 19
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    matrix, kmer_index, _ = build_pa_matrix(kmer_sets)
    conserved, _ = filter_conserved(matrix, kmer_index, min_genomes=6)
    unitigs = extract_unitigs(conserved)

    hits = find_exact_matches(unitigs, genomes)
    assert len(hits) > 0, "No mapping hits found"

    # Hits should cluster around the conserved region (1000-1300)
    hit_positions = [h.start for h in hits]
    hits_in_conserved = [p for p in hit_positions
                         if CONSERVED_START - TOLERANCE <= p <= CONSERVED_END + TOLERANCE]
    assert len(hits_in_conserved) > 0, (
        f"No hits in conserved region {CONSERVED_START}-{CONSERVED_END}. "
        f"Hit positions: {sorted(set(hit_positions))[:20]}"
    )


def test_window_scoring_finds_conserved_region():
    """Full pipeline: the top-scoring window should overlap the planted conserved region."""
    genomes = load_genomes(GENOMES_DIR)
    k = 19
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    matrix, kmer_index, _ = build_pa_matrix(kmer_sets)
    conserved, _ = filter_conserved(matrix, kmer_index, min_genomes=6)
    unitigs = extract_unitigs(conserved)
    hits = find_exact_matches(unitigs, genomes)

    genome_lengths = {name: len(seq) for name, seq in genomes.items()}
    windows = score_windows(
        hits, genome_lengths,
        window_size=300,
        min_genome_count=3,  # lower threshold for 8 genomes
        step=10,
    )

    assert len(windows) > 0, "No scored windows found"

    # The top window should overlap the planted conserved region
    top = windows[0]
    assert top.start < CONSERVED_END and top.end > CONSERVED_START, (
        f"Top window ({top.start}-{top.end}) does not overlap "
        f"conserved region ({CONSERVED_START}-{CONSERVED_END})"
    )


def test_gff3_annotation():
    gff_path = ANNOTATIONS_DIR / "TestVirus_A.gff3"
    features = parse_gff3(str(gff_path))
    assert len(features) >= 1
    ns5 = [f for f in features if f.attributes.get("Name") == "NS5"]
    assert len(ns5) >= 1
    assert ns5[0].start == 800  # 0-based (GFF3 says 801)
    assert ns5[0].end == 1600
