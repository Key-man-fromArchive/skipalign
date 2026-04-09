"""Strand orientation inference and genome normalization."""

from __future__ import annotations

from dataclasses import dataclass

from Bio.Seq import Seq

from skipalign.mapper import Hit


@dataclass
class OrientationCall:
    genome: str
    orientation: str  # "+" or "-"
    same_votes: float
    reverse_votes: float
    anchors_used: int
    confidence: float
    status: str  # "ok" | "ambiguous"


def _is_palindrome(seq: str) -> bool:
    return seq.upper() == str(Seq(seq).reverse_complement()).upper()


def _rc(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def choose_reference(
    genomes: dict[str, str],
    hits: list[Hit],
    annotations: dict[str, str] | None = None,
    user_ref: str | None = None,
) -> str:
    """Choose a reference genome for orientation."""
    if user_ref and user_ref in genomes:
        return user_ref
    if annotations:
        annotated = [g for g in genomes if g in annotations]
        if annotated:
            return max(annotated, key=lambda g: sum(1 for h in hits if h.genome == g))
    return max(genomes, key=lambda g: sum(1 for h in hits if h.genome == g))


def _index_unique_hits(hits: list[Hit]) -> dict[tuple[str, str], Hit]:
    """Index hits that are unique per (genome, unitig) pair."""
    by_key: dict[tuple[str, str], list[Hit]] = {}
    for h in hits:
        by_key.setdefault((h.genome, h.unitig), []).append(h)
    return {k: v[0] for k, v in by_key.items() if len(v) == 1}


def infer_orientations(
    genomes: dict[str, str],
    hits: list[Hit],
    annotations: dict[str, str] | None = None,
    user_ref: str | None = None,
    min_anchors: int = 3,
    min_confidence: float = 0.60,
) -> tuple[str, dict[str, OrientationCall]]:
    """Infer strand orientation of each genome relative to a reference.

    Returns:
        ref: name of the reference genome
        calls: dict mapping genome name to OrientationCall
    """
    ref = choose_reference(genomes, hits, annotations, user_ref)
    unique = _index_unique_hits(hits)

    # Reference is always "+"
    calls: dict[str, OrientationCall] = {
        ref: OrientationCall(ref, "+", 0, 0, 0, 1.0, "ok")
    }

    # Get reference's unique non-palindromic unitig hits
    ref_hits: dict[str, Hit] = {}
    for (genome, unitig), hit in unique.items():
        if genome == ref and not _is_palindrome(unitig):
            ref_hits[unitig] = hit

    for genome in genomes:
        if genome == ref:
            continue

        same_votes = 0.0
        reverse_votes = 0.0
        anchors_used = 0

        for unitig, ref_hit in ref_hits.items():
            tgt_hit = unique.get((genome, unitig))
            if not tgt_hit:
                continue

            weight = len(unitig)
            if ref_hit.strand == tgt_hit.strand:
                same_votes += weight
            else:
                reverse_votes += weight
            anchors_used += 1

        total = same_votes + reverse_votes
        if anchors_used < min_anchors or total == 0:
            calls[genome] = OrientationCall(genome, "+", same_votes, reverse_votes, anchors_used, 0.0, "ambiguous")
            continue

        if reverse_votes > same_votes:
            orientation = "-"
            conf = reverse_votes / total
        else:
            orientation = "+"
            conf = same_votes / total

        status = "ok" if conf >= min_confidence else "ambiguous"
        calls[genome] = OrientationCall(genome, orientation, same_votes, reverse_votes, anchors_used, conf, status)

    return ref, calls


def normalize_genomes(
    genomes: dict[str, str],
    calls: dict[str, OrientationCall],
) -> dict[str, str]:
    """Reverse-complement genomes that are in reverse orientation."""
    normalized = {}
    for genome, seq in genomes.items():
        call = calls.get(genome)
        if call and call.orientation == "-":
            normalized[genome] = _rc(seq)
        else:
            normalized[genome] = seq
    return normalized
