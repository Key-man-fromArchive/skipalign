"""Gradio web interface for skipalign."""

from __future__ import annotations

import shutil
import tempfile
import threading
import time
import zipfile
from pathlib import Path

import gradio as gr
import psutil

from skipalign import __version__

# Persistent storage for jobs (auto-cleaned after MAX_AGE_HOURS)
JOBS_DIR = Path(tempfile.gettempdir()) / "skipalign_jobs"
JOBS_DIR.mkdir(exist_ok=True)
MAX_AGE_HOURS = 48  # auto-delete after 2 days
MAX_UPLOAD_MB = 500


def _cleanup_old_jobs():
    """Background thread: delete job directories older than MAX_AGE_HOURS."""
    while True:
        try:
            now = time.time()
            cutoff = now - MAX_AGE_HOURS * 3600
            for job_dir in JOBS_DIR.iterdir():
                if job_dir.is_dir() and job_dir.stat().st_mtime < cutoff:
                    shutil.rmtree(job_dir, ignore_errors=True)
        except Exception:
            pass
        time.sleep(3600)  # check every hour


def _get_system_stats() -> str:
    """Get current CPU and RAM usage."""
    cpu = psutil.cpu_percent(interval=0.5)
    mem = psutil.virtual_memory()
    used_gb = mem.used / (1024 ** 3)
    total_gb = mem.total / (1024 ** 3)
    disk = psutil.disk_usage("/")
    disk_used_gb = disk.used / (1024 ** 3)
    disk_total_gb = disk.total / (1024 ** 3)

    jobs_size = sum(f.stat().st_size for f in JOBS_DIR.rglob("*") if f.is_file()) / (1024 ** 2)

    return (
        f"CPU: {cpu:.1f}% | "
        f"RAM: {used_gb:.1f}/{total_gb:.1f} GB ({mem.percent}%) | "
        f"Disk: {disk_used_gb:.0f}/{disk_total_gb:.0f} GB | "
        f"Jobs cache: {jobs_size:.1f} MB (auto-clean: {MAX_AGE_HOURS}h)"
    )


def _create_zip(output_dir: Path) -> str:
    """Zip all output files for bulk download."""
    zip_path = output_dir / "skipalign_results.zip"
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in sorted(output_dir.iterdir()):
            if f.name != "skipalign_results.zip" and f.name != "work" and f.is_file():
                zf.write(f, f.name)
    return str(zip_path)


def run_web_pipeline(
    files: list[str],
    gff_files: list[str] | None,
    k: int,
    min_genomes: int,
    window_size: int,
    design_region: int,
    top_n: int,
    top_regions: int,
) -> tuple[str, str | None, str | None, str | None, str | None, str | None, str | None]:
    """Run the pipeline from uploaded files and return results."""
    from skipalign.pipeline import run_pipeline

    if not files:
        return "No FASTA files uploaded.", None, None, None, None, None, None

    # Create job directory (persistent, cleaned by background thread)
    job_id = f"job_{int(time.time())}_{id(files)}"
    job_dir = JOBS_DIR / job_id
    input_dir = job_dir / "genomes"
    output_dir = job_dir / "results"
    input_dir.mkdir(parents=True)

    # Copy uploaded FASTA files
    for f in files:
        src = Path(f)
        shutil.copy(src, input_dir / src.name)

    # Copy GFF3 files if provided
    ann_dir = None
    if gff_files:
        ann_dir_path = job_dir / "annotations"
        ann_dir_path.mkdir()
        for f in gff_files:
            shutil.copy(Path(f), ann_dir_path / Path(f).name)
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
        return f"Pipeline failed: {e}", None, None, None, None, None, None

    if "error" in summary:
        return f"Pipeline stopped: {summary['error']}", None, None, None, None, None, None

    # Collect output files
    results_text = _format_summary(summary)

    report = output_dir / "report.html"
    primers = output_dir / "primers.tsv"
    regions = output_dir / "candidate_regions.tsv"
    msa_files = sorted(output_dir.glob("msa_region_*.fasta"))
    conserved_files = sorted(output_dir.glob("conserved_region_*.fasta"))

    # Create zip for bulk download
    zip_path = _create_zip(output_dir)

    return (
        results_text,
        str(report) if report.exists() else None,
        str(primers) if primers.exists() else None,
        str(msa_files[0]) if msa_files else None,
        str(conserved_files[0]) if conserved_files else None,
        str(regions) if regions.exists() else None,
        zip_path,
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

        # System monitor
        system_stats = gr.Textbox(label="System", interactive=False, lines=1)

        with gr.Row():
            with gr.Column(scale=2):
                fasta_input = gr.File(
                    label="Genome FASTA files (drag & drop, up to 500MB)",
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
        refresh_btn = gr.Button("Refresh System Stats", size="sm")

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
            zip_file = gr.File(label="All Results (ZIP)")

        run_btn.click(
            fn=run_web_pipeline,
            inputs=[
                fasta_input, gff_input,
                k_input, min_genomes_input, window_input,
                design_input, top_n_input, top_regions_input,
            ],
            outputs=[
                results_text, report_file, primers_file,
                msa_file, conserved_file, regions_file, zip_file,
            ],
        )

        refresh_btn.click(fn=_get_system_stats, outputs=system_stats)
        app.load(fn=_get_system_stats, outputs=system_stats)

    return app


def main():
    # Start cleanup thread
    cleaner = threading.Thread(target=_cleanup_old_jobs, daemon=True)
    cleaner.start()

    app = create_app()
    app.launch(
        server_name="0.0.0.0",
        server_port=7860,
        share=False,
        max_file_size=f"{MAX_UPLOAD_MB}mb",
    )


if __name__ == "__main__":
    main()
