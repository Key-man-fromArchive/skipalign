"""TaqMan primer-probe design with degeneracy-aware rules."""

from __future__ import annotations

import subprocess
import shutil
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.Seq import Seq

import primer3

IUPAC_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "GC", "W": "AT",
    "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
    "H": "ACT", "V": "ACG", "N": "ACGT",
}

BASES_TO_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "AC": "M", "AG": "R", "AT": "W", "CG": "S", "CT": "Y", "GT": "K",
    "ACG": "V", "ACT": "H", "AGT": "D", "CGT": "B",
    "ACGT": "N",
}

# Design constants from the paper
MAX_DEGENERACY = 100
CONSENSUS_THRESHOLD = 0.50


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
    forward_gc: float = 0.0
    reverse_gc: float = 0.0
    probe_gc: float = 0.0
    forward_issues: list[str] = field(default_factory=list)
    reverse_issues: list[str] = field(default_factory=list)
    probe_issues: list[str] = field(default_factory=list)


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


def extract_conserved_region(
    genomes: dict[str, str],
    region_start: int,
    region_end: int,
    output_fasta: Path,
) -> None:
    """Extract a conserved region from all genomes and write as multi-FASTA."""
    with open(output_fasta, "w") as f:
        for name, seq in sorted(genomes.items()):
            start = max(0, region_start)
            end = min(len(seq), region_end)
            region = seq[start:end]
            f.write(f">{name}:{start}-{end}\n{region}\n")


def compute_consensus(alignment_path: Path, threshold: float = CONSENSUS_THRESHOLD) -> str:
    """Compute IUPAC consensus from a MAFFT alignment.

    At each position:
    - If one base has frequency >= threshold, use that base
    - Otherwise, encode using IUPAC ambiguity codes
    """
    alignment = AlignIO.read(str(alignment_path), "fasta")
    length = alignment.get_alignment_length()
    consensus = []

    for i in range(length):
        column = [str(rec.seq[i]).upper() for rec in alignment]
        column = [b for b in column if b in "ACGT"]  # skip gaps
        if not column:
            consensus.append("N")
            continue

        counts = Counter(column)
        total = len(column)

        # Check if any single base exceeds threshold
        top_base, top_count = counts.most_common(1)[0]
        if top_count / total >= threshold:
            consensus.append(top_base)
        else:
            # Use IUPAC ambiguity code for bases above a minimum presence
            present = sorted(b for b, c in counts.items() if c / total >= 0.1)
            key = "".join(present)
            consensus.append(BASES_TO_IUPAC.get(key, "N"))

    return "".join(consensus)


def compute_conservation_scores(alignment_path: Path) -> list[float]:
    """Compute per-position conservation score (fraction of most common base)."""
    alignment = AlignIO.read(str(alignment_path), "fasta")
    length = alignment.get_alignment_length()
    scores = []

    for i in range(length):
        column = [str(rec.seq[i]).upper() for rec in alignment]
        column = [b for b in column if b in "ACGT"]
        if not column:
            scores.append(0.0)
            continue
        counts = Counter(column)
        scores.append(counts.most_common(1)[0][1] / len(column))

    return scores


def _gc_content(seq: str) -> float:
    """Calculate GC content of a sequence (handling IUPAC codes)."""
    gc = sum(1 for b in seq.upper() if b in "GCS")
    return gc / len(seq) * 100 if seq else 0.0


def _reverse_complement(seq: str) -> str:
    """Reverse complement of a DNA sequence."""
    return str(Seq(seq).reverse_complement())


def design_primers_from_consensus(
    consensus: str,
    conservation_scores: list[float],
    top_n: int = 5,
) -> list[PrimerProbeSet]:
    """Design primer-probe sets using primer3 on the consensus sequence.

    Uses primer3 to generate candidates, then validates against TaqMan rules.
    """
    # Strip gaps and replace IUPAC ambiguity codes with N for primer3
    # (primer3 only accepts ACGTN)
    clean_consensus = ""
    for base in consensus:
        if base in "ACGTacgt":
            clean_consensus += base.upper()
        elif base != "-":
            clean_consensus += "N"

    if len(clean_consensus) < 80:
        return []

    # Use primer3 to design primer pairs with internal oligo (probe)
    design_result = primer3.design_primers(
        seq_args={
            "SEQUENCE_ID": "conserved_region",
            "SEQUENCE_TEMPLATE": clean_consensus,
        },
        global_args={
            "PRIMER_NUM_RETURN": top_n * 3,  # request more, filter later
            "PRIMER_PRODUCT_SIZE_RANGE": [[70, 200]],
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_OPT_SIZE": 22,
            "PRIMER_MAX_SIZE": 28,
            "PRIMER_MIN_TM": 50.0,
            "PRIMER_OPT_TM": 57.0,
            "PRIMER_MAX_TM": 63.0,
            "PRIMER_INTERNAL_MIN_SIZE": 18,
            "PRIMER_INTERNAL_OPT_SIZE": 24,
            "PRIMER_INTERNAL_MAX_SIZE": 30,
            "PRIMER_INTERNAL_MIN_TM": 55.0,
            "PRIMER_INTERNAL_OPT_TM": 60.0,
            "PRIMER_INTERNAL_MAX_TM": 68.0,
            "PRIMER_PICK_INTERNAL_OLIGO": 1,
        },
    )

    num_returned = design_result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    candidates: list[PrimerProbeSet] = []

    for i in range(num_returned):
        fwd = design_result.get(f"PRIMER_LEFT_{i}_SEQUENCE", "")
        rev = design_result.get(f"PRIMER_RIGHT_{i}_SEQUENCE", "")
        probe = design_result.get(f"PRIMER_INTERNAL_{i}_SEQUENCE", "")
        fwd_tm = design_result.get(f"PRIMER_LEFT_{i}_TM", 0.0)
        rev_tm = design_result.get(f"PRIMER_RIGHT_{i}_TM", 0.0)
        probe_tm = design_result.get(f"PRIMER_INTERNAL_{i}_TM", 0.0)
        product_size = design_result.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)

        if not (fwd and rev and probe):
            continue

        pps = PrimerProbeSet(
            forward=fwd,
            reverse=rev,
            probe=probe,
            amplicon_length=product_size,
            forward_tm=round(fwd_tm, 1),
            reverse_tm=round(rev_tm, 1),
            probe_tm=round(probe_tm, 1),
            forward_degeneracy=degeneracy(fwd),
            reverse_degeneracy=degeneracy(rev),
            probe_degeneracy=degeneracy(probe),
            forward_gc=round(_gc_content(fwd), 1),
            reverse_gc=round(_gc_content(rev), 1),
            probe_gc=round(_gc_content(probe), 1),
            forward_issues=validate_forward_primer(fwd),
            reverse_issues=validate_reverse_primer(rev),
            probe_issues=validate_probe(probe),
        )
        candidates.append(pps)

    # Sort: fewest total issues first, then by amplicon length
    candidates.sort(key=lambda p: (
        len(p.forward_issues) + len(p.reverse_issues) + len(p.probe_issues),
        p.amplicon_length,
    ))

    return candidates[:top_n]


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
