"""Tests for De Bruijn graph and unitig extraction."""

from skipalign.unitig import build_debruijn, extract_unitigs


def test_build_debruijn_simple():
    # k=3: ABC -> prefix=AB, suffix=BC
    kmers = ["ABC", "BCD", "CDE"]
    graph = build_debruijn(kmers)
    assert "AB" in graph
    assert "BC" in graph["AB"]


def test_extract_unitigs_linear_path():
    # Linear chain: ABCDE as overlapping 3-mers
    kmers = ["ABC", "BCD", "CDE"]
    unitigs = extract_unitigs(kmers)
    assert len(unitigs) == 1
    assert unitigs[0] == "ABCDE"


def test_extract_unitigs_branching():
    # Branch at B: ABC->BCD and ABC->BCX
    kmers = ["ABC", "BCD", "BCX"]
    unitigs = extract_unitigs(kmers)
    assert len(unitigs) >= 2


def test_extract_unitigs_empty():
    assert extract_unitigs([]) == []
