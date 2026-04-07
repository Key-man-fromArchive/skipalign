"""TaqMan primer-probe design with degeneracy-aware rules."""

from __future__ import annotations

import subprocess
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from tempfile import NamedTemporaryFile

from Bio import SeqIO

IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}


def degeneracy(seq: str) -> int:
    """Calculate total degeneracy of an IUPAC sequence."""
    result = 1
    for base in seq.upper():
        result *= len(IUPAC_EXPAND.get(base, base))
    return result


@dataclass
class PrimerProbeSet:
    forward: str
    reverse: str
    probe: str
    amplicon_length: int
    forward_tm: float = 0.0
    reverse_tm: float = 0.0
    probe_tm: float = 0.0
    forward_degeneracy: int = 0
    reverse_degeneracy: int = 0
    probe_degeneracy: int = 0


def check_mafft() -> bool:
    """Check if MAFFT is available in PATH."""
    return shutil.which("mafft") is not None


def run_mafft(input_fasta: Path, output_fasta: Path) -> None:
    """Run MAFFT alignment on input FASTA."""
    if not check_mafft():
        raise RuntimeError(
            "MAFFT not found in PATH. Install via: conda install -c bioconda mafft"
        )
    result = subprocess.run(
        ["mafft", "--auto", str(input_fasta)],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"MAFFT failed: {result.stderr}")
    output_fasta.write_text(result.stdout)


def validate_forward_primer(seq: str) -> list[str]:
    """Validate forward primer against TaqMan design rules."""
    issues = []
    if degeneracy(seq) > 100:
        issues.append(f"Degeneracy {degeneracy(seq)} exceeds cap of 100")
    if "GGGG" in seq.upper():
        issues.append("Contains >4 consecutive G residues")
    # Check 3' end conservation (last 5 bases should have minimal degeneracy)
    tail = seq[-5:]
    if degeneracy(tail) > 4:
        issues.append(f"3' end degeneracy {degeneracy(tail)} too high")
    return issues


def validate_probe(seq: str) -> list[str]:
    """Validate probe against TaqMan design rules."""
    issues = []
    if degeneracy(seq) > 100:
        issues.append(f"Degeneracy {degeneracy(seq)} exceeds cap of 100")
    if seq[0].upper() == "G":
        issues.append("Probe starts with 5' G")
    if len(seq) > 1 and seq[1].upper() == "G":
        issues.append("G at second position from 5' end")
    if "GGGG" in seq.upper():
        issues.append("Contains >4 consecutive G residues")
    if "AAAAAA" in seq.upper():
        issues.append("Contains >=6 consecutive A residues")
    return issues


def validate_reverse_primer(seq: str) -> list[str]:
    """Validate reverse primer against TaqMan design rules."""
    issues = []
    if degeneracy(seq) > 100:
        issues.append(f"Degeneracy {degeneracy(seq)} exceeds cap of 100")
    if "GGGG" in seq.upper():
        issues.append("Contains >4 consecutive G residues")
    return issues
