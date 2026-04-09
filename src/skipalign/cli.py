"""CLI entry point for skipalign."""

from __future__ import annotations

import typer

from skipalign import __version__

app = typer.Typer(
    name="skipalign",
    help="Alignment-free conserved region discovery and RT-qPCR primer-probe design.",
    no_args_is_help=True,
)


@app.command()
def run(
    input_dir: str = typer.Option(..., "--input", "-i", help="Directory containing genome FASTA files"),
    annotations_dir: str = typer.Option(None, "--annotations", "-g", help="Directory containing GFF3 annotation files"),
    k: int = typer.Option(19, "--k", "-k", help="k-mer length"),
    min_genomes: int = typer.Option(3, "--min-genomes", help="Minimum genomes for high-confidence unitigs"),
    window: int = typer.Option(300, "--window", "-w", help="Scoring window size in bp"),
    design_region: int = typer.Option(600, "--design-region", help="Design region extraction size in bp"),
    output_dir: str = typer.Option("results", "--output", "-o", help="Output directory"),
    top: int = typer.Option(5, "--top", help="Number of top primer candidates per region"),
    top_regions: int = typer.Option(3, "--top-regions", help="Number of candidate conserved regions"),
    step: int = typer.Option(10, "--step", help="Window scoring step size in bp"),
    no_report: bool = typer.Option(False, "--no-report", help="Skip HTML report generation"),
) -> None:
    """Run the full alignment-free primer design pipeline."""
    from skipalign.pipeline import run_pipeline

    run_pipeline(
        input_dir=input_dir,
        output_dir=output_dir,
        annotations_dir=annotations_dir,
        k=k,
        min_genomes=min_genomes,
        window_size=window,
        design_region=design_region,
        top_n=top,
        top_regions=top_regions,
        step=step,
        no_report=no_report,
    )


@app.command()
def find_k(
    input_dir: str = typer.Option(..., "--input", "-i", help="Directory containing genome FASTA files"),
    k_min: int = typer.Option(9, "--k-min", help="Minimum k to test"),
    k_max: int = typer.Option(51, "--k-max", help="Maximum k to test"),
    step: int = typer.Option(2, "--step", help="k increment"),
    output: str = typer.Option("acf_result.tsv", "--output", "-o", help="Output ACF table"),
) -> None:
    """Discover optimal k-mer length via ACF analysis."""
    from skipalign.io import load_genomes
    from skipalign.kmer import count_kmers
    from skipalign.matrix import build_pa_matrix, filter_conserved

    typer.echo(f"skipalign v{__version__} — ACF sweep k={k_min}..{k_max}")
    genomes = load_genomes(input_dir)
    typer.echo(f"Loaded {len(genomes)} genomes")

    with open(output, "w") as f:
        f.write("k\ttotal_kmers\tshared_kmers\tshared_fraction\n")
        for k_val in range(k_min, k_max + 1, step):
            kmer_sets = {name: count_kmers(seq, k_val) for name, seq in genomes.items()}
            matrix, kmer_index, _ = build_pa_matrix(kmer_sets)
            conserved, _ = filter_conserved(matrix, kmer_index, min_genomes=len(genomes))
            total = len(kmer_index)
            shared = len(conserved)
            frac = shared / total if total > 0 else 0
            f.write(f"{k_val}\t{total}\t{shared}\t{frac:.6f}\n")
            typer.echo(f"  k={k_val}: {total:,} total, {shared:,} shared ({frac:.4%})")

    typer.echo(f"\nResults written to {output}")


@app.callback(invoke_without_command=True)
def version(
    version: bool = typer.Option(False, "--version", "-V", help="Show version"),
) -> None:
    if version:
        typer.echo(f"skipalign {__version__}")
        raise typer.Exit()
