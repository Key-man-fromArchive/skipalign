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
    step: int = 1,
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

    for win_start in range(0, ref_length - window_size + 1, step):
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


@dataclass
class RegionCandidate:
    rank: int
    cluster_start: int
    cluster_end: int
    summit_start: int
    summit_end: int
    peak_score: int
    peak_genome_count: int
    area_score: int
    design_start: int
    design_end: int
    feature: str | None
    window_count: int


def cluster_windows(
    windows: list[ScoredWindow],
    genome_length: int,
    design_region: int = 600,
    merge_gap: int | None = None,
    window_size: int = 300,
    top_n: int = 5,
    max_overlap: float = 0.5,
) -> list[RegionCandidate]:
    """Cluster adjacent high-scoring windows into distinct conserved regions.

    Args:
        windows: scored windows (from score_windows)
        genome_length: reference genome length
        design_region: extraction region size (bp)
        merge_gap: max gap to merge windows (default: window_size // 3)
        top_n: number of top regions to return
        max_overlap: max design-interval overlap fraction before suppression
    """
    if not windows:
        return []

    if merge_gap is None:
        merge_gap = min(150, max(50, window_size // 3))

    # Sort by position
    ws = sorted(windows, key=lambda w: (w.start, w.end))

    # Merge into clusters
    clusters: list[list[ScoredWindow]] = []
    current = [ws[0]]
    current_end = ws[0].end

    for w in ws[1:]:
        gap = w.start - current_end
        if gap <= merge_gap:
            current.append(w)
            current_end = max(current_end, w.end)
        else:
            clusters.append(current)
            current = [w]
            current_end = w.end
    clusters.append(current)

    # Summarize each cluster
    regions = []
    for cluster in clusters:
        peak = max(cluster, key=lambda w: (w.genome_count, w.coverage_score))

        # Core windows: within 85% of peak score and within 1 of peak genome count
        core = [
            w for w in cluster
            if w.coverage_score >= peak.coverage_score * 0.85
            and w.genome_count >= peak.genome_count - 1
        ]
        if not core:
            core = cluster

        cluster_start = min(w.start for w in core)
        cluster_end = max(w.end for w in core)

        # Weighted summit center
        weight_sum = sum(w.coverage_score for w in core)
        if weight_sum > 0:
            summit_center = int(sum(((w.start + w.end) / 2) * w.coverage_score for w in core) / weight_sum)
        else:
            summit_center = (cluster_start + cluster_end) // 2

        design_start = max(0, summit_center - design_region // 2)
        design_end = min(genome_length, design_start + design_region)
        design_start = max(0, design_end - design_region)

        # Dominant feature
        features = [w.feature for w in core if w.feature]
        feature = max(set(features), key=features.count) if features else None

        regions.append(RegionCandidate(
            rank=0,
            cluster_start=cluster_start,
            cluster_end=cluster_end,
            summit_start=peak.start,
            summit_end=peak.end,
            peak_score=peak.coverage_score,
            peak_genome_count=peak.genome_count,
            area_score=sum(w.coverage_score for w in core),
            design_start=design_start,
            design_end=design_end,
            feature=feature,
            window_count=len(cluster),
        ))

    # Rank by genome support → peak score → area score → narrower
    regions.sort(
        key=lambda r: (r.peak_genome_count, r.peak_score, r.area_score, -(r.cluster_end - r.cluster_start)),
        reverse=True,
    )

    # Non-maximum suppression on design intervals
    chosen: list[RegionCandidate] = []
    for r in regions:
        overlap_ok = all(
            _overlap_fraction(r.design_start, r.design_end, c.design_start, c.design_end) < max_overlap
            for c in chosen
        )
        if overlap_ok:
            r.rank = len(chosen) + 1
            chosen.append(r)
        if len(chosen) >= top_n:
            break

    return chosen


def _overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    ov = max(0, min(a_end, b_end) - max(a_start, b_start))
    shorter = min(a_end - a_start, b_end - b_start)
    return ov / shorter if shorter > 0 else 0.0
