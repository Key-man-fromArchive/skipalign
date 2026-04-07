"""Tests for primer-probe validation rules."""

from skipalign.primer import degeneracy, validate_forward_primer, validate_probe


def test_degeneracy_no_ambiguity():
    assert degeneracy("ACGT") == 1


def test_degeneracy_with_r():
    assert degeneracy("R") == 2  # A or G


def test_degeneracy_with_n():
    assert degeneracy("N") == 4


def test_degeneracy_mixed():
    assert degeneracy("ARN") == 1 * 2 * 4  # 8


def test_validate_forward_primer_clean():
    issues = validate_forward_primer("ACGTACGTACGTACGTACGTACG")
    assert issues == []


def test_validate_forward_primer_g_run():
    issues = validate_forward_primer("ACGTGGGGGACGT")
    assert any(">4 consecutive G" in i for i in issues)


def test_validate_probe_5prime_g():
    issues = validate_probe("GACGTACGTACGT")
    assert any("5' G" in i for i in issues)


def test_validate_probe_a_run():
    issues = validate_probe("TACAAAAAAGT")
    assert any(">=6 consecutive A" in i for i in issues)
