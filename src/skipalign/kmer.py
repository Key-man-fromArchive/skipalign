"""k-mer counting and canonical form utilities."""


def canonical(kmer: str) -> str:
    """Return the lexicographically smaller of a k-mer and its reverse complement."""
    complement = str.maketrans("ACGT", "TGCA")
    rc = kmer.translate(complement)[::-1]
    return min(kmer, rc)


def count_kmers(sequence: str, k: int) -> set[str]:
    """Extract all canonical k-mers from a sequence, skipping non-ACGT characters."""
    valid = set("ACGT")
    kmers = set()
    seq = sequence.upper()
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if all(c in valid for c in kmer):
            kmers.add(canonical(kmer))
    return kmers
