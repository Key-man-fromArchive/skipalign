"""CLI entry point for skipalign."""

import typer

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
    output_dir: str = typer.Option("results", "--output", "-o", help="Output directory"),
) -> None:
    """Run the full alignment-free primer design pipeline."""
    typer.echo(f"skipalign v{__version__} — full pipeline")
    typer.echo(f"Input: {input_dir}, k={k}, window={window}bp")
    # TODO: implement full pipeline orchestration


@app.command()
def find_k(
    input_dir: str = typer.Option(..., "--input", "-i", help="Directory containing genome FASTA files"),
    k_min: int = typer.Option(9, "--k-min", help="Minimum k to test"),
    k_max: int = typer.Option(51, "--k-max", help="Maximum k to test"),
    output: str = typer.Option("acf_result.tsv", "--output", "-o", help="Output ACF table"),
) -> None:
    """Discover optimal k-mer length via ACF analysis."""
    typer.echo(f"ACF sweep k={k_min}..{k_max}")
    # TODO: implement ACF sweep


from skipalign import __version__  # noqa: E402
