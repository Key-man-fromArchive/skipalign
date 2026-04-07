"""FASTA and GFF3 I/O utilities."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO


def load_genomes(input_dir: str | Path) -> dict[str, str]:
    """Load all FASTA files from a directory.

    Returns:
        dict mapping genome name (filename stem) to uppercase sequence string.
    """
    input_path = Path(input_dir)
    genomes: dict[str, str] = {}

    for fasta_file in sorted(input_path.glob("*.fasta")) + sorted(input_path.glob("*.fa")):
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = fasta_file.stem
            genomes[name] = str(record.seq).upper()
            break  # Take the first record per file (complete genome)

    if not genomes:
        raise FileNotFoundError(f"No FASTA files found in {input_path}")

    return genomes


def load_all_gff3(annotations_dir: str | Path) -> dict[str, str]:
    """Load all GFF3 file paths from a directory.

    Returns:
        dict mapping genome name (filename stem) to GFF3 file path.
    """
    ann_path = Path(annotations_dir)
    gff_files: dict[str, str] = {}

    for gff_file in sorted(ann_path.glob("*.gff3")) + sorted(ann_path.glob("*.gff")):
        gff_files[gff_file.stem] = str(gff_file)

    return gff_files
