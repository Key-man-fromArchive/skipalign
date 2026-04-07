"""Full pipeline orchestration — connects all modules end-to-end."""

from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path

from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn

from skipalign import __version__
from skipalign.io import load_genomes, load_all_gff3
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
from skipalign.scorer import ScoredWindow, score_windows
from skipalign.unitig import extract_unitigs

console = Console()


def run_pipeline(
    input_dir: str,
    output_dir: str,
    annotations_dir: str | None = None,
    k: int = 19,
    min_genomes: int = 3,
    window_size: int = 300,
    design_region: int = 600,
    top_n: int = 5,
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

    # Step 1: Load genomes
    console.print("  [1/6] Loading genomes ...", end=" ")
    genomes = load_genomes(input_dir)
    console.print(f"[green]✓[/green]  {len(genomes)} genomes loaded")

    # Step 2: Count k-mers and build PA matrix
    console.print("  [2/6] Counting k-mers & building PA matrix ...", end=" ")
    kmer_sets = {name: count_kmers(seq, k) for name, seq in genomes.items()}
    matrix, kmer_index, genome_names = build_pa_matrix(kmer_sets)
    conserved_kmers, genome_counts = filter_conserved(matrix, kmer_index, min_genomes=min_genomes)
    total_unique = len(kmer_index)
    console.print(f"[green]✓[/green]  {total_unique:,} unique {k}-mers, {len(conserved_kmers):,} conserved")

    if not conserved_kmers:
        console.print("[red]  ✗ No conserved k-mers found. Try lowering --min-genomes or adjusting --k[/red]")
        return {"error": "no_conserved_kmers"}

    # Step 3: Extract unitigs
    console.print("  [3/6] Extracting unitigs ...", end=" ")
    unitigs = extract_unitigs(conserved_kmers)
    console.print(f"[green]✓[/green]  {len(unitigs)} unitigs ({min(len(u) for u in unitigs)}-{max(len(u) for u in unitigs)}bp)")

    if not unitigs:
        console.print("[red]  ✗ No unitigs extracted[/red]")
        return {"error": "no_unitigs"}

    # Step 4: Map unitigs to genomes + annotate
    console.print("  [4/6] Mapping to genomes ...", end=" ")
    hits = find_exact_matches(unitigs, genomes)

    # Annotate with GFF3 if available
    if annotations_dir:
        gff_files = load_all_gff3(annotations_dir)
        for genome_name, gff_path in gff_files.items():
            features = parse_gff3(gff_path)
            genome_hits = [h for h in hits if h.genome == genome_name]
            annotate_hits(genome_hits, features)

    console.print(f"[green]✓[/green]  {len(hits):,} hits across {len(set(h.genome for h in hits))} genomes")

    # Step 5: Score windows
    console.print("  [5/6] Scoring windows ...", end=" ")
    genome_lengths = {name: len(seq) for name, seq in genomes.items()}
    windows = score_windows(hits, genome_lengths, window_size=window_size, min_genome_count=min_genomes, step=step)
    console.print(f"[green]✓[/green]  {len(windows)} windows above threshold")

    if not windows:
        console.print("[red]  ✗ No conserved windows found. Try lowering --min-genomes[/red]")
        return {"error": "no_windows"}

    # Determine design region from top windows
    top_window = windows[0]
    region_center = (top_window.start + top_window.end) // 2
    region_start = max(0, region_center - design_region // 2)
    region_end = region_start + design_region

    # Step 6: Primer design
    console.print("  [6/6] Designing primers ...", end=" ")

    primers: list[PrimerProbeSet] = []
    msa_path = work_path / "msa_alignment.fasta"
    region_fasta = work_path / "conserved_region.fasta"

    extract_conserved_region(genomes, region_start, region_end, region_fasta)

    if check_mafft():
        run_mafft(region_fasta, msa_path)
        consensus = compute_consensus(msa_path)
        conservation = compute_conservation_scores(msa_path)
        primers = design_primers_from_consensus(consensus, conservation, top_n=top_n)

        # Copy MSA to output
        (out_path / "msa_alignment.fasta").write_text(msa_path.read_text())
        console.print(f"[green]✓[/green]  {len(primers)} candidate sets")
    else:
        console.print("[yellow]⚠[/yellow]  MAFFT not found — skipping primer design")
        console.print("         Install via: [bold]conda install -c bioconda mafft[/bold]")
        consensus = ""
        conservation = []

    elapsed = time.time() - start_time

    # Write outputs
    _write_conserved_region(out_path, region_fasta)
    _write_primers_tsv(out_path / "primers.tsv", primers)
    summary = _write_summary(
        out_path / "pipeline_summary.json",
        genomes=genomes,
        k=k,
        window_size=window_size,
        min_genomes=min_genomes,
        total_unique_kmers=total_unique,
        conserved_kmers=len(conserved_kmers),
        unitigs_count=len(unitigs),
        hits_count=len(hits),
        windows_count=len(windows),
        top_window=top_window,
        region_start=region_start,
        region_end=region_end,
        primer_count=len(primers),
        elapsed=elapsed,
    )

    console.print(f"\n  [green]✓ Pipeline complete in {elapsed:.1f}s[/green]")
    console.print(f"  Results: {out_path}/")

    return summary


def _write_conserved_region(out_path: Path, region_fasta: Path) -> None:
    import shutil
    shutil.copy(region_fasta, out_path / "conserved_region.fasta")


def _write_primers_tsv(path: Path, primers: list[PrimerProbeSet]) -> None:
    header = "set_id\ttype\tsequence\tlength\ttm\tgc_percent\tdegeneracy\tissues\n"
    with open(path, "w") as f:
        f.write(header)
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
    conserved_kmers: int,
    unitigs_count: int,
    hits_count: int,
    windows_count: int,
    top_window: ScoredWindow,
    region_start: int,
    region_end: int,
    primer_count: int,
    elapsed: float,
) -> dict:
    summary = {
        "version": __version__,
        "params": {"k": k, "window": window_size, "min_genomes": min_genomes},
        "input": {
            "genomes": len(genomes),
            "total_bp": sum(len(s) for s in genomes.values()),
        },
        "results": {
            "unique_kmers": total_unique_kmers,
            "conserved_kmers": conserved_kmers,
            "unitigs": unitigs_count,
            "mapping_hits": hits_count,
            "scored_windows": windows_count,
            "conserved_region": {
                "start": region_start,
                "end": region_end,
                "feature": top_window.feature,
                "genome_count": top_window.genome_count,
                "coverage_score": top_window.coverage_score,
            },
            "primer_candidates": primer_count,
        },
        "runtime_seconds": round(elapsed, 2),
    }
    path.write_text(json.dumps(summary, indent=2))
    return summary
