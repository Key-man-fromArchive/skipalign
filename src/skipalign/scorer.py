"""Sliding window scoring for conserved region discovery."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from skipalign.mapper import Hit


@dataclass
class ScoredWindow:
    feature: str | None
    start: int
    end: int
    coverage_score: int
    genome_count: int


def score_windows(
    hits: list[Hit],
    genome_lengths: dict[str, int],
    window_size: int = 300,
    min_genome_count: int = 15,
) -> list[ScoredWindow]:
    """Score sliding windows by unitig coverage and genome count.

    For each genome, slides a window across the sequence and scores by:
      1. Total unitig bp covered within the window
      2. Number of distinct genomes with unitig hits in the window

    Returns windows passing the min_genome_count threshold, sorted by score.
    """
    if not hits:
        return []

    # Group hits by genome
    hits_by_genome: dict[str, list[Hit]] = {}
    for hit in hits:
        hits_by_genome.setdefault(hit.genome, []).append(hit)

    # Use a reference genome to define window positions
    # (pick the one with the most hits for robustness)
    ref_genome = max(hits_by_genome, key=lambda g: len(hits_by_genome[g]))
    ref_length = genome_lengths[ref_genome]

    scored: list[ScoredWindow] = []

    for win_start in range(0, ref_length - window_size + 1, 1):
        win_end = win_start + window_size
        total_coverage = 0
        genomes_with_hits: set[str] = set()

        for genome_name, genome_hits in hits_by_genome.items():
            for hit in genome_hits:
                # Check overlap with window
                overlap_start = max(hit.start, win_start)
                overlap_end = min(hit.end, win_end)
                if overlap_start < overlap_end:
                    total_coverage += overlap_end - overlap_start
                    genomes_with_hits.add(genome_name)

        if len(genomes_with_hits) >= min_genome_count:
            feature = None
            for hit in hits_by_genome.get(ref_genome, []):
                if hit.start < win_end and hit.end > win_start and hit.feature:
                    feature = hit.feature
                    break

            scored.append(ScoredWindow(
                feature=feature,
                start=win_start,
                end=win_end,
                coverage_score=total_coverage,
                genome_count=len(genomes_with_hits),
            ))

    # Sort by genome_count (primary) then coverage_score (secondary), descending
    scored.sort(key=lambda w: (w.genome_count, w.coverage_score), reverse=True)
    return scored
