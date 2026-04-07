"""Generate synthetic test genomes with a planted conserved region.

Creates 8 synthetic genomes (~2000bp) where:
- Positions 1000-1300 contain a 300bp conserved region (>95% identity)
- The rest is divergent random sequence
- A GFF3 annotation marks a "NS5-like" gene spanning positions 800-1600

This allows the pipeline to be tested end-to-end: the conserved region
should be discoverable by k-mer analysis and window scoring.
"""

import random
from pathlib import Path

random.seed(42)

GENOME_LENGTH = 2000
CONSERVED_START = 1000
CONSERVED_END = 1300
CONSERVED_LENGTH = CONSERVED_END - CONSERVED_START
GENE_START = 800
GENE_END = 1600
NUM_GENOMES = 8
MUTATION_RATE_CONSERVED = 0.03  # 3% mutation in conserved region
GENOME_NAMES = [
    "TestVirus_A", "TestVirus_B", "TestVirus_C", "TestVirus_D",
    "TestVirus_E", "TestVirus_F", "TestVirus_G", "TestVirus_H",
]


def random_seq(length: int) -> str:
    return "".join(random.choice("ACGT") for _ in range(length))


def mutate(seq: str, rate: float) -> str:
    bases = list(seq)
    for i in range(len(bases)):
        if random.random() < rate:
            bases[i] = random.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def main():
    out_dir = Path(__file__).parent / "data"
    genomes_dir = out_dir / "genomes"
    annotations_dir = out_dir / "annotations"
    genomes_dir.mkdir(parents=True, exist_ok=True)
    annotations_dir.mkdir(parents=True, exist_ok=True)

    # Generate the "ancestral" conserved region (shared across all genomes)
    conserved_ancestor = random_seq(CONSERVED_LENGTH)

    for name in GENOME_NAMES:
        # Random flanking regions (unique per genome)
        left_flank = random_seq(CONSERVED_START)
        right_flank = random_seq(GENOME_LENGTH - CONSERVED_END)

        # Conserved region with slight mutations
        conserved = mutate(conserved_ancestor, MUTATION_RATE_CONSERVED)

        genome = left_flank + conserved + right_flank

        # Write FASTA
        fasta_path = genomes_dir / f"{name}.fasta"
        with open(fasta_path, "w") as f:
            f.write(f">{name} synthetic test genome\n")
            for i in range(0, len(genome), 80):
                f.write(genome[i:i+80] + "\n")

    # Write GFF3 for first genome (to test annotation feature)
    gff_path = annotations_dir / f"{GENOME_NAMES[0]}.gff3"
    with open(gff_path, "w") as f:
        f.write("##gff-version 3\n")
        f.write(f"{GENOME_NAMES[0]}\t.\tgene\t{GENE_START+1}\t{GENE_END}\t.\t+\t.\tID=gene_NS5;Name=NS5\n")
        f.write(f"{GENOME_NAMES[0]}\t.\tCDS\t{GENE_START+1}\t{GENE_END}\t.\t+\t0\tID=cds_NS5;Parent=gene_NS5;Name=NS5;product=RNA-dependent RNA polymerase\n")

    # Write expected results for validation
    expected_path = out_dir / "expected_results.txt"
    with open(expected_path, "w") as f:
        f.write(f"conserved_region_start={CONSERVED_START}\n")
        f.write(f"conserved_region_end={CONSERVED_END}\n")
        f.write(f"conserved_region_length={CONSERVED_LENGTH}\n")
        f.write(f"gene_name=NS5\n")
        f.write(f"num_genomes={NUM_GENOMES}\n")
        f.write(f"conserved_ancestor={conserved_ancestor}\n")

    print(f"Generated {NUM_GENOMES} synthetic genomes in {genomes_dir}")
    print(f"Conserved region: {CONSERVED_START}-{CONSERVED_END} ({CONSERVED_LENGTH}bp)")
    print(f"Gene annotation: NS5 ({GENE_START}-{GENE_END})")
    print(f"GFF3: {gff_path}")


if __name__ == "__main__":
    main()
