"""Gradio web interface for skipalign."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path

import gradio as gr

from skipalign import __version__


def run_web_pipeline(
    files: list[str],
    gff_files: list[str] | None,
    k: int,
    min_genomes: int,
    window_size: int,
    design_region: int,
    top_n: int,
    top_regions: int,
) -> tuple[str, str | None, str | None, str | None, str | None, str | None]:
    """Run the pipeline from uploaded files and return results."""
    from skipalign.pipeline import run_pipeline

    if not files:
        return "No FASTA files uploaded.", None, None, None, None, None

    # Create temp directories for input and output
    with tempfile.TemporaryDirectory() as tmpdir:
        input_dir = Path(tmpdir) / "genomes"
        input_dir.mkdir()
        output_dir = Path(tmpdir) / "results"

        # Copy uploaded FASTA files
        for f in files:
            src = Path(f)
            dst = input_dir / src.name
            shutil.copy(src, dst)

        # Copy GFF3 files if provided
        ann_dir = None
        if gff_files:
            ann_dir_path = Path(tmpdir) / "annotations"
            ann_dir_path.mkdir()
            for f in gff_files:
                src = Path(f)
                shutil.copy(src, ann_dir_path / src.name)
            ann_dir = str(ann_dir_path)

        # Run pipeline
        try:
            summary = run_pipeline(
                input_dir=str(input_dir),
                output_dir=str(output_dir),
                annotations_dir=ann_dir,
                k=k,
                min_genomes=min_genomes,
                window_size=window_size,
                design_region=design_region,
                top_n=top_n,
                top_regions=top_regions,
                step=10,
            )
        except Exception as e:
            return f"Pipeline failed: {e}", None, None, None, None, None

        if "error" in summary:
            return f"Pipeline stopped: {summary['error']}", None, None, None, None, None

        # Read results
        results_text = _format_summary(summary)

        # Collect output files for download
        report_path = output_dir / "report.html"
        primers_path = output_dir / "primers.tsv"
        regions_path = output_dir / "candidate_regions.tsv"
        msa_files = sorted(output_dir.glob("msa_region_*.fasta"))
        conserved_files = sorted(output_dir.glob("conserved_region_*.fasta"))

        # Copy to persistent temp files (gradio needs them to persist)
        downloads = {}
        for name, src in [
            ("report", report_path),
            ("primers", primers_path),
            ("regions", regions_path),
        ]:
            if src.exists():
                dst = Path(tempfile.mktemp(suffix=f"_{src.name}"))
                shutil.copy(src, dst)
                downloads[name] = str(dst)

        # Bundle MSA + conserved region files
        if msa_files:
            dst = Path(tempfile.mktemp(suffix=f"_{msa_files[0].name}"))
            shutil.copy(msa_files[0], dst)
            downloads["msa"] = str(dst)

        if conserved_files:
            dst = Path(tempfile.mktemp(suffix=f"_{conserved_files[0].name}"))
            shutil.copy(conserved_files[0], dst)
            downloads["conserved"] = str(dst)

        return (
            results_text,
            downloads.get("report"),
            downloads.get("primers"),
            downloads.get("msa"),
            downloads.get("conserved"),
            downloads.get("regions"),
        )


def _format_summary(summary: dict) -> str:
    r = summary.get("results", {})
    s = summary.get("strand_normalization", {})
    lines = [
        f"skipalign v{summary.get('version', '?')}",
        f"Runtime: {summary.get('runtime_seconds', '?')}s",
        "",
        f"Input: {summary['input']['genomes']} genomes, {summary['input']['total_bp']:,} bp",
        f"Reference: {s.get('reference', '?')}, {s.get('reversed_genomes', 0)} RC'd",
        "",
        f"Unique {summary['params']['k']}-mers: {r.get('unique_kmers', 0):,}",
        f"Conserved k-mers: {r.get('conserved_kmers', 0):,}",
        f"Unitigs: {r.get('unitigs', 0)}",
        f"Mapping hits: {r.get('mapping_hits', 0):,}",
        f"Scored windows: {r.get('scored_windows', 0)}",
        "",
        "Candidate Regions:",
    ]
    for region in r.get("candidate_regions", []):
        lines.append(
            f"  #{region['rank']}: {region['design_start']}-{region['design_end']}bp "
            f"({region['feature'] or 'N/A'}), {region['genome_count']} genomes, "
            f"score={region['peak_score']}"
        )
    lines.append("")
    lines.append(f"Primer candidates: {r.get('primer_candidates', 0)}")

    val = summary.get("validation", [])
    if val:
        lines.append("")
        lines.append("MFEprimer Validation:")
        for v in val:
            lines.append(f"  {v['primer_set']}: {v['coverage_percent']}% coverage "
                         f"({v['hit_sequences']}/{v['total_sequences']})")

    return "\n".join(lines)


def create_app() -> gr.Blocks:
    with gr.Blocks(title="skipalign") as app:
        gr.Markdown(f"""
        # skipalign v{__version__}
        **Alignment-free conserved region discovery for RT-qPCR primer-probe design**

        Upload genome FASTA files to discover conserved regions across highly divergent viral genera.
        """)

        with gr.Row():
            with gr.Column(scale=2):
                fasta_input = gr.File(
                    label="Genome FASTA files",
                    file_count="multiple",
                    file_types=[".fasta", ".fa", ".fna"],
                )
                gff_input = gr.File(
                    label="GFF3 annotations (optional)",
                    file_count="multiple",
                    file_types=[".gff3", ".gff"],
                )

            with gr.Column(scale=1):
                k_input = gr.Slider(9, 51, value=19, step=2, label="k-mer length")
                min_genomes_input = gr.Slider(2, 50, value=3, step=1, label="Min genomes")
                window_input = gr.Slider(100, 600, value=300, step=50, label="Window size (bp)")
                design_input = gr.Slider(300, 1200, value=600, step=100, label="Design region (bp)")
                top_n_input = gr.Slider(1, 10, value=5, step=1, label="Primer candidates")
                top_regions_input = gr.Slider(1, 10, value=3, step=1, label="Candidate regions")

        run_btn = gr.Button("Run Pipeline", variant="primary", size="lg")

        gr.Markdown("### Results")
        results_text = gr.Textbox(label="Pipeline Summary", lines=20)

        gr.Markdown("### Downloads")
        with gr.Row():
            report_file = gr.File(label="HTML Report")
            primers_file = gr.File(label="Primers TSV")
            msa_file = gr.File(label="MSA (FASTA)")
        with gr.Row():
            conserved_file = gr.File(label="Conserved Region")
            regions_file = gr.File(label="Candidate Regions TSV")

        run_btn.click(
            fn=run_web_pipeline,
            inputs=[
                fasta_input, gff_input,
                k_input, min_genomes_input, window_input,
                design_input, top_n_input, top_regions_input,
            ],
            outputs=[
                results_text, report_file, primers_file,
                msa_file, conserved_file, regions_file,
            ],
        )

    return app


def main():
    app = create_app()
    app.launch(server_name="0.0.0.0", server_port=7860, share=False)


if __name__ == "__main__":
    main()
