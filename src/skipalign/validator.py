"""MFEprimer in-silico PCR validation."""

from __future__ import annotations

import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from skipalign.primer import PrimerProbeSet


@dataclass
class ValidationResult:
    primer_set_id: str
    forward: str
    reverse: str
    total_db_sequences: int
    hit_sequences: int
    coverage_percent: float
    amplicons: list[AmpliconHit] = field(default_factory=list)


@dataclass
class AmpliconHit:
    target: str
    size: int
    forward_tm: float
    reverse_tm: float
    product_tm: float
    position: str  # e.g. "TestVirus_A:1084-1177"


def check_mfeprimer() -> bool:
    """Check if MFEprimer is available in PATH."""
    return shutil.which("mfeprimer") is not None


def index_database(fasta_path: Path) -> None:
    """Index a FASTA database for MFEprimer."""
    index_file = Path(f"{fasta_path}.primerqc")
    if index_file.exists():
        return  # already indexed

    result = subprocess.run(
        ["mfeprimer", "index", "-i", str(fasta_path)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"MFEprimer index failed: {result.stderr}")


def _count_db_sequences(fasta_path: Path) -> int:
    """Count number of sequences in a FASTA file."""
    count = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def _write_primer_tsv(
    primers: list[PrimerProbeSet],
    output_path: Path,
) -> None:
    """Write primers in MFEprimer TSV format: name \\t forward \\t reverse."""
    with open(output_path, "w") as f:
        for i, pps in enumerate(primers, 1):
            f.write(f"PS_{i:03d}\t{pps.forward}\t{pps.reverse}\n")


def run_validation(
    primers: list[PrimerProbeSet],
    db_fasta: Path,
    output_dir: Path,
    max_amplicon_size: int = 500,
) -> list[ValidationResult]:
    """Run MFEprimer specificity check on primer-probe sets.

    Args:
        primers: list of PrimerProbeSet candidates
        db_fasta: path to genome FASTA (will be indexed if needed)
        output_dir: directory for output files
        max_amplicon_size: max amplicon size to search for
    """
    if not check_mfeprimer():
        raise RuntimeError(
            "MFEprimer not found in PATH. "
            "Download from: https://github.com/quwubin/MFEprimer-3.0/releases"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    # Index DB
    index_database(db_fasta)

    # Write primer TSV
    primer_tsv = output_dir / "primers_mfe.tsv"
    _write_primer_tsv(primers, primer_tsv)

    # Run MFEprimer spec
    result_txt = output_dir / "mfeprimer_result.txt"
    result_json = output_dir / "mfeprimer_result.txt.json"

    result = subprocess.run(
        [
            "mfeprimer", "spec",
            "-i", str(primer_tsv),
            "-d", str(db_fasta),
            "-o", str(result_txt),
            "-j",
            "-S", str(max_amplicon_size),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"MFEprimer spec failed: {result.stderr}")

    # Parse JSON results
    total_db = _count_db_sequences(db_fasta)
    return _parse_results(result_json, primers, total_db)


def _parse_results(
    json_path: Path,
    primers: list[PrimerProbeSet],
    total_db_sequences: int,
) -> list[ValidationResult]:
    """Parse MFEprimer JSON output into ValidationResult objects."""
    with open(json_path) as f:
        data = json.load(f)

    amp_list = data.get("AmpList") or []

    # Group amplicons by primer set ID extracted from F.Seq.ID (e.g. "PS_001_fp")
    amp_by_set: dict[str, list[dict]] = {}
    for amp in amp_list:
        fwd_id = amp.get("F", {}).get("Seq", {}).get("ID", "")
        # "PS_001_fp" → "PS_001"
        set_id = fwd_id.rsplit("_fp", 1)[0] if "_fp" in fwd_id else "unknown"
        amp_by_set.setdefault(set_id, []).append(amp)

    results = []
    for i, pps in enumerate(primers, 1):
        set_id = f"PS_{i:03d}"
        set_amps = amp_by_set.get(set_id, [])

        # Count unique target sequences hit
        hit_targets: set[str] = set()
        amplicon_hits: list[AmpliconHit] = []

        for amp in set_amps:
            hit_id = amp.get("Hid", "")
            hit_targets.add(hit_id)

            f_start = amp.get("F", {}).get("Start", 0)
            r_end = amp.get("R", {}).get("End", 0)
            f_tm = amp.get("F", {}).get("Tm", 0.0)
            r_tm = amp.get("R", {}).get("Tm", 0.0)

            amplicon_hits.append(AmpliconHit(
                target=hit_id,
                size=amp.get("Size", 0),
                forward_tm=round(f_tm, 1),
                reverse_tm=round(r_tm, 1),
                product_tm=round(amp.get("GC", 0.0), 1),  # product Tm not directly available, use GC
                position=f"{hit_id}:{f_start}-{r_end}",
            ))

        coverage = len(hit_targets) / total_db_sequences * 100 if total_db_sequences > 0 else 0

        results.append(ValidationResult(
            primer_set_id=set_id,
            forward=pps.forward,
            reverse=pps.reverse,
            total_db_sequences=total_db_sequences,
            hit_sequences=len(hit_targets),
            coverage_percent=round(coverage, 1),
            amplicons=amplicon_hits,
        ))

    return results


def format_validation_summary(results: list[ValidationResult]) -> str:
    """Format validation results as a human-readable summary."""
    lines = ["MFEprimer In-silico PCR Validation", "=" * 40, ""]
    for r in results:
        lines.append(f"{r.primer_set_id}: {r.hit_sequences}/{r.total_db_sequences} "
                      f"genomes ({r.coverage_percent}% coverage)")
        lines.append(f"  Forward:  {r.forward}")
        lines.append(f"  Reverse:  {r.reverse}")
        if r.amplicons:
            sizes = sorted(set(a.size for a in r.amplicons))
            lines.append(f"  Amplicons: {len(r.amplicons)} hits, sizes {min(sizes)}-{max(sizes)}bp")
        lines.append("")
    return "\n".join(lines)
