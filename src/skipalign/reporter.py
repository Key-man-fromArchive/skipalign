"""HTML + TSV report generation for skipalign results."""

from __future__ import annotations

import base64
import io
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

from skipalign.mapper import Hit
from skipalign.primer import PrimerProbeSet
from skipalign.scorer import ScoredWindow


def _fig_to_base64(fig: plt.Figure) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("utf-8")


def plot_conservation_landscape(
    windows: list[ScoredWindow],
    region_start: int,
    region_end: int,
    genome_count_total: int,
) -> str:
    """Generate conservation landscape plot as base64 PNG."""
    if not windows:
        return ""

    fig, ax = plt.subplots(figsize=(12, 4))
    starts = [w.start for w in windows]
    scores = [w.genome_count for w in windows]

    ax.bar(starts, scores, width=max(10, (max(starts) - min(starts)) / len(starts)),
           color="#4a90d9", alpha=0.7, label="Genome count")

    # Highlight selected region
    ax.axvspan(region_start, region_end, alpha=0.15, color="red", label="Design region")

    ax.set_xlabel("Genomic position (bp)")
    ax.set_ylabel("Genomes with unitig hits")
    ax.set_title("Conservation Landscape — Window Scoring")
    ax.legend(loc="upper right")
    ax.set_ylim(0, genome_count_total + 1)

    return _fig_to_base64(fig)


def plot_conservation_heatmap(
    conservation_scores: list[float],
    region_start: int,
) -> str:
    """Generate MSA conservation heatmap as base64 PNG."""
    if not conservation_scores:
        return ""

    fig, ax = plt.subplots(figsize=(12, 1.5))
    scores = np.array(conservation_scores).reshape(1, -1)
    im = ax.imshow(scores, aspect="auto", cmap="RdYlGn", vmin=0, vmax=1,
                   extent=[region_start, region_start + len(conservation_scores), 0, 1])
    ax.set_xlabel("Position (bp)")
    ax.set_yticks([])
    ax.set_title("MSA Conservation (green = conserved, red = variable)")
    fig.colorbar(im, ax=ax, label="Conservation", shrink=0.6)

    return _fig_to_base64(fig)


def generate_html_report(
    output_path: Path,
    summary: dict,
    windows: list[ScoredWindow],
    primers: list[PrimerProbeSet],
    conservation_scores: list[float],
    region_start: int,
    region_end: int,
    validation_results: list | None = None,
) -> None:
    """Generate a self-contained HTML report."""
    genome_count = summary["input"]["genomes"]

    landscape_img = plot_conservation_landscape(windows, region_start, region_end, genome_count)
    heatmap_img = plot_conservation_heatmap(conservation_scores, region_start)

    primers_html = _render_primers_table(primers)
    summary_html = _render_summary_table(summary)
    validation_html = _render_validation_table(validation_results) if validation_results else ""
    windows_html = _render_windows_table(windows[:20])

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>skipalign Report</title>
<style>
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
         background: #f5f7fa; color: #1a1a2e; line-height: 1.6; padding: 2rem; }}
  .container {{ max-width: 1100px; margin: 0 auto; }}
  h1 {{ font-size: 1.8rem; margin-bottom: 0.5rem; color: #16213e; }}
  h2 {{ font-size: 1.3rem; margin: 2rem 0 1rem; color: #0f3460;
        border-bottom: 2px solid #e2e8f0; padding-bottom: 0.5rem; }}
  .subtitle {{ color: #64748b; margin-bottom: 2rem; }}
  .card {{ background: #fff; border-radius: 8px; padding: 1.5rem;
           margin-bottom: 1.5rem; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }}
  table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
  th {{ background: #f1f5f9; text-align: left; padding: 0.6rem 0.8rem;
       border-bottom: 2px solid #e2e8f0; font-weight: 600; }}
  td {{ padding: 0.5rem 0.8rem; border-bottom: 1px solid #f1f5f9; }}
  tr:hover {{ background: #f8fafc; }}
  .issue {{ color: #dc2626; font-size: 0.8rem; }}
  .ok {{ color: #16a34a; }}
  img {{ max-width: 100%; height: auto; }}
  .stats {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 1rem; margin-bottom: 1.5rem; }}
  .stat {{ background: #fff; border-radius: 8px; padding: 1rem;
           box-shadow: 0 1px 3px rgba(0,0,0,0.08); text-align: center; }}
  .stat-value {{ font-size: 1.8rem; font-weight: 700; color: #0f3460; }}
  .stat-label {{ font-size: 0.8rem; color: #64748b; }}
  .footer {{ text-align: center; color: #94a3b8; font-size: 0.8rem; margin-top: 3rem; }}
</style>
</head>
<body>
<div class="container">
  <h1>skipalign Report</h1>
  <p class="subtitle">Alignment-free conserved region discovery &amp; primer-probe design &mdash; v{summary['version']}</p>

  <div class="stats">
    <div class="stat"><div class="stat-value">{summary['input']['genomes']}</div><div class="stat-label">Genomes</div></div>
    <div class="stat"><div class="stat-value">{summary['results']['unique_kmers']:,}</div><div class="stat-label">{summary['params']['k']}-mers</div></div>
    <div class="stat"><div class="stat-value">{summary['results']['unitigs']}</div><div class="stat-label">Unitigs</div></div>
    <div class="stat"><div class="stat-value">{summary['results']['primer_candidates']}</div><div class="stat-label">Primer Sets</div></div>
    <div class="stat"><div class="stat-value">{summary['runtime_seconds']}s</div><div class="stat-label">Runtime</div></div>
  </div>

  <h2>1. Pipeline Summary</h2>
  <div class="card">{summary_html}</div>

  <h2>2. Conservation Landscape</h2>
  <div class="card">
    {"<img src='data:image/png;base64," + landscape_img + "' alt='Conservation landscape'>" if landscape_img else "<p>No data</p>"}
  </div>

  <h2>3. MSA Conservation</h2>
  <div class="card">
    {"<img src='data:image/png;base64," + heatmap_img + "' alt='Conservation heatmap'>" if heatmap_img else "<p>MAFFT not available</p>"}
  </div>

  <h2>4. Primer-Probe Candidates</h2>
  <div class="card">{primers_html if primers_html else "<p>No candidates generated</p>"}</div>

  <h2>5. In-silico PCR Validation (MFEprimer)</h2>
  <div class="card">{validation_html if validation_html else "<p>MFEprimer not available — skipped</p>"}</div>

  <h2>6. Top Scored Windows</h2>
  <div class="card">{windows_html}</div>

  <div class="footer">
    Generated by <strong>skipalign</strong> v{summary['version']}
    &mdash; k={summary['params']['k']}, window={summary['params']['window']}bp, min_genomes={summary['params']['min_genomes']}
  </div>
</div>
</body>
</html>"""

    output_path.write_text(html)


def _render_summary_table(summary: dict) -> str:
    r = summary["results"]
    region = r.get("conserved_region", {})
    rows = [
        ("Input genomes", f"{summary['input']['genomes']}"),
        ("Total base pairs", f"{summary['input']['total_bp']:,}"),
        ("k-mer length", f"{summary['params']['k']}"),
        ("Unique k-mers", f"{r['unique_kmers']:,}"),
        ("Conserved k-mers", f"{r['conserved_kmers']:,}"),
        ("Unitigs", f"{r['unitigs']}"),
        ("Mapping hits", f"{r['mapping_hits']:,}"),
        ("Scored windows", f"{r['scored_windows']}"),
        ("Conserved region", f"{region.get('start', '?')}-{region.get('end', '?')} bp"),
        ("Feature", f"{region.get('feature') or 'N/A'}"),
        ("Genome coverage", f"{region.get('genome_count', '?')}/{summary['input']['genomes']}"),
        ("Primer candidates", f"{r['primer_candidates']}"),
    ]
    html = "<table><tr><th>Metric</th><th>Value</th></tr>"
    for label, value in rows:
        html += f"<tr><td>{label}</td><td>{value}</td></tr>"
    html += "</table>"
    return html


def _render_primers_table(primers: list[PrimerProbeSet]) -> str:
    if not primers:
        return ""
    html = """<table>
    <tr><th>Set</th><th>Type</th><th>Sequence</th><th>Length</th>
    <th>Tm (°C)</th><th>GC%</th><th>Degeneracy</th><th>Issues</th></tr>"""
    for i, pps in enumerate(primers, 1):
        for oligo_type, seq, tm, gc, deg, issues in [
            ("Forward", pps.forward, pps.forward_tm, pps.forward_gc, pps.forward_degeneracy, pps.forward_issues),
            ("Probe", pps.probe, pps.probe_tm, pps.probe_gc, pps.probe_degeneracy, pps.probe_issues),
            ("Reverse", pps.reverse, pps.reverse_tm, pps.reverse_gc, pps.reverse_degeneracy, pps.reverse_issues),
        ]:
            issue_str = "; ".join(issues) if issues else '<span class="ok">OK</span>'
            if issues:
                issue_str = f'<span class="issue">{issue_str}</span>'
            html += f"<tr><td>PS_{i:03d}</td><td>{oligo_type}</td><td><code>{seq}</code></td>"
            html += f"<td>{len(seq)}</td><td>{tm}</td><td>{gc}</td><td>{deg}</td><td>{issue_str}</td></tr>"
    html += "</table>"
    return html


def _render_validation_table(validation_results: list) -> str:
    if not validation_results:
        return ""
    html = """<table>
    <tr><th>Primer Set</th><th>Forward</th><th>Reverse</th>
    <th>Hits</th><th>Total</th><th>Coverage</th><th>Amplicons</th></tr>"""
    for v in validation_results:
        cov_class = "ok" if v.coverage_percent >= 80 else ("issue" if v.coverage_percent < 50 else "")
        cov_str = f'<span class="{cov_class}">{v.coverage_percent}%</span>' if cov_class else f"{v.coverage_percent}%"
        sizes = sorted(set(a.size for a in v.amplicons)) if v.amplicons else []
        size_str = f"{min(sizes)}-{max(sizes)}bp" if sizes else "—"
        html += f"<tr><td>{v.primer_set_id}</td>"
        html += f"<td><code>{v.forward}</code></td>"
        html += f"<td><code>{v.reverse}</code></td>"
        html += f"<td>{v.hit_sequences}</td><td>{v.total_db_sequences}</td>"
        html += f"<td>{cov_str}</td><td>{len(v.amplicons)} ({size_str})</td></tr>"
    html += "</table>"
    return html


def _render_windows_table(windows: list[ScoredWindow]) -> str:
    if not windows:
        return "<p>No windows</p>"
    html = """<table>
    <tr><th>Rank</th><th>Start</th><th>End</th><th>Coverage Score</th>
    <th>Genome Count</th><th>Feature</th></tr>"""
    for i, w in enumerate(windows, 1):
        html += f"<tr><td>{i}</td><td>{w.start}</td><td>{w.end}</td>"
        html += f"<td>{w.coverage_score:,}</td><td>{w.genome_count}</td>"
        html += f"<td>{w.feature or 'N/A'}</td></tr>"
    html += "</table>"
    return html
