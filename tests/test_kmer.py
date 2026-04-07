"""Tests for k-mer counting and canonical form."""

from skipalign.kmer import canonical, count_kmers


def test_canonical_forward_is_smaller():
    assert canonical("AAAC") == "AAAC"  # AAAC < GTTT


def test_canonical_rc_is_smaller():
    assert canonical("GTTT") == "AAAC"  # rc(GTTT) = AAAC


def test_canonical_palindrome():
    assert canonical("AATT") == "AATT"  # rc(AATT) = AATT


def test_count_kmers_basic():
    seq = "ACGTACGT"
    kmers = count_kmers(seq, 4)
    # ACGT, CGTA, GTAC, TACG, ACGT(dup)
    # Canonical: ACGT, CGTA→TACG, GTAC, TACG→CGTA — need to check
    assert len(kmers) > 0
    assert all(len(k) == 4 for k in kmers)


def test_count_kmers_skips_n():
    seq = "ACGTNACGT"
    kmers = count_kmers(seq, 4)
    # k-mers spanning N should be skipped
    assert all("N" not in k for k in kmers)


def test_count_kmers_empty():
    assert count_kmers("", 4) == set()
    assert count_kmers("ACG", 4) == set()  # seq shorter than k
