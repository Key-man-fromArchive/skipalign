"""Full pipeline orchestration — connects all modules end-to-end."""

from __future__ import annotations

import json
import time
from pathlib import Path

from rich.console import Console

from skipalign import __version__
from skipalign.io import load_all_gff3, load_genomes
from skipalign.kmer import count_kmers
from skipalign.mapper import Hit, annotate_hits, find_exact_matches, parse_gff3
from skipalign.matrix import build_pa_matrix, filter_conserved
from skipalign.primer import (
    PrimerProbeSet,
    check_mafft,
    compute_consensus,
    compute_conservation_scores,
    design_primers_from_consensus,
    extract_conserved_region,
    run_mafft,
)
from skipalign.reporter import generate_html_report
from skipalign.scorer import RegionCandidate, ScoredWindow, cluster_windows, score_windows
from skipalign.strand import infer_orientations, normalize_genomes
from skipalign.unitig import extract_unitigs
from skipalign.validator import (
    check_mfeprimer as check_mfeprimer_bin,
    format_validation_summary,
    run_validation,
)

console = Console()


def _localize_region_per_genome(
    genomes: dict[str, str],
    hits: list[Hit],
    ref_start: int,
    ref_end: int,
    design_region: int,
) -> dict[str, tuple[int, int]]:
    """Localize each genome's extraction interval using anchor hits.

    Instead of blindly extracting the same X:Y from every genome,
    find the median anchor position within the reference region for
    each genome and center the extraction around that.
    """
    # Collect unitigs that overlap the reference region
    region_unitigs: set[str] = set()
    for h in hits:
        if h.start < ref_end and h.end > ref_start:
            region_unitigs.add(h.unitig)

    intervals: dict[str, tuple[int, int]] = {}
    for genome, seq in genomes.items():
        genome_hits = [h for h in hits if h.genome == genome and h.unitig in region_unitigs]
        if not genome_hits:
            # Fallback: use reference coordinates
            start = max(0, min(ref_start, len(seq) - design_region))
            intervals[genome] = (start, start + design_region)
            continue

        centers = sorted((h.start + h.end) // 2 for h in genome_hits)
        median_center = centers[len(centers) // 2]
        start = max(0, median_center - design_region // 2)
        end = min(len(seq), start + design_region)
        start = max(0, end - design_region)
        intervals[genome] = (start, end)

    return intervals


def _extract_region_per_genome(
    genomes: dict[str, str],
    intervals: dict[str, tuple[int, int]],
    output_fasta: Path,
) -> None:
    """Extract genome-specific regions to a multi-FASTA file."""
    with open(output_fasta, "w") as f:
        for name in sorted(genomes.keys()):
            seq = genomes[name]
            start, end = intervals.get(name, (0, min(600, len(seq))))
            region = seq[start:end]
            f.write(f">{name}:{start}-{end}\n{region}\n")


def run_pipeline(
    input_dir: str,
    output_dir: str,
    annotations_dir: str | None = None,
    k: int = 19,
    min_genomes: int = 3,
    window_size: int = 300,
    design_region: int = 600,
    top_n: int = 5,
    top_regions: int = 3,
    step: int = 10,
    no_report: bool = False,
    verbose: bool = False,
) -> dict:
    """Run the full alignment-free primer design pipeline."""
    start_time = time.time()
    out_path = Path(output_dir)
    work_path = out_path / "work"
    out_path.mkdir(parents=True, exist_ok=True)
    work_path.mkdir(parents=True, exist_ok=True)

    console.print(f"\n[bold]skipalign v{__version__}[/bold] — Alignment-free primer design\n")

    # ── Step 1: Load genomes ──
    console.print("  [1/7] Loading genomes ...", end=" ")
    raw_genomes = load_genomes(input_dir)
    console.print(f"[green]✓[/green]  {len(raw_genomes)} genomes loaded")

    # ── Step 2: k-mer discovery (strand-agnostic) ──
    console.print("  [2/7] Counting k-mers & building PA matrix ...", end=" ")
    kmer_sets = {name: count_kmers(seq, k) for name, seq in raw_genomes.items()}
    matrix, kmer_index, genome_names = build_pa_matrix(kmer_sets)
    conserved_kmers, genome_counts = filter_conserved(matrix, kmer_index, min_genomes=min_genomes)
    total_unique = len(kmer_index)
    console.print(f"[green]✓[/green]  {total_unique:,} unique {k}-mers, {len(conserved_kmers):,} conserved")

    if not conserved_kmers:
        console.print("[red]  ✗ No conserved k-mers found. Try lowering --min-genomes or adjusting --k[/red]")
        return {"error": "no_conserved_kmers"}

    # ── Step 3: Unitigs ──
    console.print("  [3/7] Extracting unitigs ...", end=" ")
    unitigs = extract_unitigs(conserved_kmers)
    console.print(f"[green]✓[/green]  {len(unitigs)} unitigs ({min(len(u) for u in unitigs)}-{max(len(u) for u in unitigs)}bp)")

    if not unitigs:
        console.print("[red]  ✗ No unitigs extracted[/red]")
        return {"error": "no_unitigs"}

    # ── Step 4: Strand normalization (two-pass) ──
    console.print("  [4/7] Normalizing strand orientation ...", end=" ")
    raw_hits = find_exact_matches(unitigs, raw_genomes)

    gff_files = load_all_gff3(annotations_dir) if annotations_dir else {}
    ref_genome, orientation_calls = infer_orientations(
        raw_genomes, raw_hits, annotations=gff_files if gff_files else None,
    )
    genomes = normalize_genomes(raw_genomes, orientation_calls)

    rc_count = sum(1 for c in orientation_calls.values() if c.orientation == "-")
    ambig_count = sum(1 for c in orientation_calls.values() if c.status == "ambiguous")
    console.print(f"[green]✓[/green]  ref={ref_genome}, {rc_count} RC'd, {ambig_count} ambiguous")

    # Re-map on normalized genomes
    hits = find_exact_matches(unitigs, genomes)
    if gff_files:
        for gname, gff_path in gff_files.items():
            features = parse_gff3(gff_path)
            genome_hits = [h for h in hits if h.genome == gname]
            annotate_hits(genome_hits, features)

    # ── Step 5: Window scoring + clustering ──
    console.print("  [5/7] Scoring & clustering windows ...", end=" ")
    genome_lengths = {name: len(seq) for name, seq in genomes.items()}
    windows = score_windows(hits, genome_lengths, window_size=window_size, min_genome_count=min_genomes, step=step)

    if not windows:
        console.print("[red]  ✗ No conserved windows found. Try lowering --min-genomes[/red]")
        return {"error": "no_windows"}

    ref_length = genome_lengths[ref_genome]
    regions = cluster_windows(
        windows, ref_length,
        design_region=design_region,
        window_size=window_size,
        top_n=top_regions,
    )
    console.print(f"[green]✓[/green]  {len(windows)} windows → {len(regions)} candidate regions")

    # ── Step 6: Template extraction + primer suggestions per region ──
    console.print("  [6/7] Extracting templates & designing primers ...", end=" ")
    all_primers: list[tuple[RegionCandidate, list[PrimerProbeSet], list[float]]] = []

    for region in regions:
        # Per-genome localized extraction
        intervals = _localize_region_per_genome(
            genomes, hits, region.design_start, region.design_end, design_region,
        )
        region_fasta = work_path / f"region_{region.rank}.fasta"
        _extract_region_per_genome(genomes, intervals, region_fasta)

        primers: list[PrimerProbeSet] = []
        conservation: list[float] = []

        if check_mafft():
            msa_path = work_path / f"msa_region_{region.rank}.fasta"
            run_mafft(region_fasta, msa_path)
            consensus = compute_consensus(msa_path)
            conservation = compute_conservation_scores(msa_path)
            primers = design_primers_from_consensus(consensus, conservation, top_n=top_n)

            # Copy to output
            import shutil
            shutil.copy(region_fasta, out_path / f"conserved_region_{region.rank}.fasta")
            shutil.copy(msa_path, out_path / f"msa_region_{region.rank}.fasta")

        all_primers.append((region, primers, conservation))

    total_primer_sets = sum(len(p) for _, p, _ in all_primers)
    if check_mafft():
        console.print(f"[green]✓[/green]  {len(regions)} regions, {total_primer_sets} primer sets total")
    else:
        console.print("[yellow]⚠[/yellow]  MAFFT not found — templates extracted, primers skipped")

    # ── Step 7: MFEprimer validation ──
    validation_results = []
    if total_primer_sets > 0 and check_mfeprimer_bin():
        console.print("  [7/7] Validating with MFEprimer ...", end=" ")
        db_fasta = work_path / "all_genomes.fasta"
        with open(db_fasta, "w") as f:
            for name, seq in sorted(genomes.items()):
                f.write(f">{name}\n{seq}\n")

        # Validate primers from the top region
        top_region, top_primers, _ = all_primers[0]
        if top_primers:
            max_amp = design_region + 300
            validation_results = run_validation(top_primers, db_fasta, work_path / "validation", max_amplicon_size=max_amp)
            coverages = [f"{v.coverage_percent}%" for v in validation_results]
            console.print(f"[green]✓[/green]  coverage: {', '.join(coverages)}")
            (out_path / "validation_summary.txt").write_text(format_validation_summary(validation_results))
        else:
            console.print("[yellow]⚠[/yellow]  No primers to validate")
    elif total_primer_sets > 0:
        console.print("  [7/7] [yellow]⚠[/yellow]  MFEprimer not found — skipping validation")

    elapsed = time.time() - start_time

    # ── Write outputs ──
    _write_regions_tsv(out_path / "candidate_regions.tsv", regions)
    if all_primers:
        _, top_primers_list, _ = all_primers[0]
        _write_primers_tsv(out_path / "primers.tsv", top_primers_list)
    summary = _write_summary(
        out_path / "pipeline_summary.json",
        genomes=genomes,
        k=k,
        window_size=window_size,
        min_genomes=min_genomes,
        total_unique_kmers=total_unique,
        conserved_kmers_count=len(conserved_kmers),
        unitigs_count=len(unitigs),
        hits_count=len(hits),
        windows_count=len(windows),
        regions=regions,
        all_primers=all_primers,
        elapsed=elapsed,
        validation_results=validation_results,
        ref_genome=ref_genome,
        rc_count=rc_count,
    )

    # HTML report (uses top region)
    if not no_report and all_primers:
        top_region, top_primers_list, top_conservation = all_primers[0]
        generate_html_report(
            output_path=out_path / "report.html",
            summary=summary,
            windows=windows,
            primers=top_primers_list,
            conservation_scores=top_conservation,
            region_start=top_region.design_start,
            region_end=top_region.design_end,
            validation_results=validation_results,
        )

    console.print(f"\n  [green]✓ Pipeline complete in {elapsed:.1f}s[/green]")
    console.print(f"  Results: {out_path}/")

    return summary


def _write_regions_tsv(path: Path, regions: list[RegionCandidate]) -> None:
    with open(path, "w") as f:
        f.write("rank\tcluster_start\tcluster_end\tdesign_start\tdesign_end\t"
                "peak_score\tpeak_genome_count\tarea_score\tfeature\twindow_count\n")
        for r in regions:
            f.write(f"{r.rank}\t{r.cluster_start}\t{r.cluster_end}\t"
                    f"{r.design_start}\t{r.design_end}\t{r.peak_score}\t"
                    f"{r.peak_genome_count}\t{r.area_score}\t{r.feature or 'N/A'}\t{r.window_count}\n")


def _write_primers_tsv(path: Path, primers: list[PrimerProbeSet]) -> None:
    with open(path, "w") as f:
        f.write("set_id\ttype\tsequence\tlength\ttm\tgc_percent\tdegeneracy\tissues\n")
        for i, pps in enumerate(primers, 1):
            for oligo_type, seq, tm, gc, deg, issues in [
                ("forward", pps.forward, pps.forward_tm, pps.forward_gc, pps.forward_degeneracy, pps.forward_issues),
                ("reverse", pps.reverse, pps.reverse_tm, pps.reverse_gc, pps.reverse_degeneracy, pps.reverse_issues),
                ("probe", pps.probe, pps.probe_tm, pps.probe_gc, pps.probe_degeneracy, pps.probe_issues),
            ]:
                issues_str = "; ".join(issues) if issues else ""
                f.write(f"PS_{i:03d}\t{oligo_type}\t{seq}\t{len(seq)}\t{tm}\t{gc}\t{deg}\t{issues_str}\n")


def _write_summary(
    path: Path,
    genomes: dict,
    k: int,
    window_size: int,
    min_genomes: int,
    total_unique_kmers: int,
    conserved_kmers_count: int,
    unitigs_count: int,
    hits_count: int,
    windows_count: int,
    regions: list[RegionCandidate],
    all_primers: list,
    elapsed: float,
    validation_results: list | None = None,
    ref_genome: str = "",
    rc_count: int = 0,
) -> dict:
    summary = {
        "version": __version__,
        "params": {"k": k, "window": window_size, "min_genomes": min_genomes},
        "input": {
            "genomes": len(genomes),
            "total_bp": sum(len(s) for s in genomes.values()),
        },
        "strand_normalization": {
            "reference": ref_genome,
            "reversed_genomes": rc_count,
        },
        "results": {
            "unique_kmers": total_unique_kmers,
            "conserved_kmers": conserved_kmers_count,
            "unitigs": unitigs_count,
            "mapping_hits": hits_count,
            "scored_windows": windows_count,
            "candidate_regions": [
                {
                    "rank": r.rank,
                    "design_start": r.design_start,
                    "design_end": r.design_end,
                    "feature": r.feature,
                    "genome_count": r.peak_genome_count,
                    "peak_score": r.peak_score,
                }
                for r in regions
            ],
            "primer_candidates": sum(len(p) for _, p, _ in all_primers),
        },
        "runtime_seconds": round(elapsed, 2),
    }
    if validation_results:
        summary["validation"] = [
            {
                "primer_set": v.primer_set_id,
                "hit_sequences": v.hit_sequences,
                "total_sequences": v.total_db_sequences,
                "coverage_percent": v.coverage_percent,
            }
            for v in validation_results
        ]
    path.write_text(json.dumps(summary, indent=2))
    return summary
