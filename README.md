# skipalign

Alignment-free conserved region discovery and RT-qPCR primer-probe design for highly divergent viral genera.

Based on the alignment-free k-mer guided approach from [Sayasit et al. (2026)](https://doi.org/10.64898/2026.03.17.712358).

## What it does

Traditional primer design requires multiple sequence alignment (MSA) across all target genomes — this breaks down when nucleotide divergence exceeds 25-30%. **skipalign** uses k-mer analysis and compacted De Bruijn graphs to discover conserved regions *without* alignment, then applies MSA only on the small conserved region for final primer-probe design.

```
Input: FASTA genomes → k-mer counting → presence-absence matrix
→ De Bruijn graph → unitigs → genome mapping → window scoring
→ conserved region → MAFFT (local MSA) → primer-probe candidates
```

## Install

### pip (requires MAFFT in PATH)

```bash
pip install -e .
skipalign --help
```

### Conda (includes MAFFT)

```bash
conda env create -f environment.yml
conda activate skipalign
skipalign --help
```

### Docker

```bash
docker build -t skipalign .
docker run -v ./data:/data skipalign run -i /data/genomes -o /data/results
```

## Usage

```bash
# Full pipeline
skipalign run \
  --input genomes/ \
  --annotations gff3/ \
  --k 19 \
  --min-genomes 3 \
  --window 300 \
  --output results/

# Discover optimal k (when targeting a new viral genus)
skipalign find-k --input genomes/ --k-min 9 --k-max 51 -o acf_result.tsv
```

## Default parameters

| Parameter | Default | Source |
|-----------|---------|--------|
| k-mer length | 19 | ACF analysis of 51 Orthoflavivirus genomes |
| Window size | 300 bp | ~2-4x typical TaqMan amplicon |
| Min genomes | 3 | High-confidence unitig threshold |
| Degeneracy cap | 100 | Per-oligonucleotide variant limit |
| Consensus threshold | 50% | Per-position base frequency |

## External dependency

**MAFFT** is the only external binary required (used in the final MSA step). Install via `conda install -c bioconda mafft` or `apt-get install mafft`.

## License

MIT
