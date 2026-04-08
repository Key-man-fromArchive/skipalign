# skipalign

[![CI](https://github.com/Key-man-fromArchive/skipalign/actions/workflows/ci.yml/badge.svg)](https://github.com/Key-man-fromArchive/skipalign/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Alignment-free conserved region discovery and RT-qPCR primer-probe design for highly divergent viral genera.

Based on the alignment-free k-mer guided approach from [Sayasit et al. (2026)](https://doi.org/10.64898/2026.03.17.712358).

## Why skipalign?

Traditional primer design requires multiple sequence alignment (MSA) across all target genomes — this breaks down when nucleotide divergence exceeds 25-30%. **skipalign** discovers conserved regions *without* alignment using k-mer analysis and compacted De Bruijn graphs, then applies MSA only on the small conserved region for final primer-probe design.

This makes it possible to design pan-genus diagnostic assays for highly variable viral groups like Orthoflavivirus (dengue, Zika, JEV).

## Pipeline

![Pipeline Overview](docs/pipeline_overview.svg)

The pipeline discovers conserved regions without global alignment, then designs TaqMan primer-probe sets with in-silico validation:

1. **k-mer Counting** — Decompose genomes into canonical 19-mers
2. **De Bruijn Graph** — Merge overlapping conserved k-mers into unitigs
3. **Genome Mapping** — Map unitigs back to genomic coordinates
4. **Window Scoring** — Score 300bp windows by conservation across genomes
5. **Primer Design** — MAFFT local MSA + primer3 with TaqMan rules
6. **Validate & Report** — MFEprimer in-silico PCR + HTML report

## Install

### pip

```bash
pip install -e .
```

Requires [MAFFT](https://mafft.cbrc.jp/alignment/software/) in PATH for primer design.
Optionally install [MFEprimer](https://github.com/quwubin/MFEprimer-3.0/releases) for in-silico PCR validation.

### Conda (includes MAFFT)

```bash
conda env create -f environment.yml
conda activate skipalign
```

### Docker (batteries included)

```bash
docker build -t skipalign .
docker run -v ./data:/data skipalign run -i /data/genomes -o /data/results
```

## Quick Start

```bash
# Place your genome FASTA files in a folder
ls genomes/
# DENV1.fasta  DENV2.fasta  ZIKV.fasta  JEV.fasta ...

# Run the full pipeline
skipalign run --input genomes/ --output results/

# Check results
open results/report.html        # Interactive HTML report
cat results/primers.tsv          # Primer-probe candidates
cat results/pipeline_summary.json # Pipeline statistics
```

### Terminal Output

```
skipalign v0.1.0 — Alignment-free primer design

  [1/6] Loading genomes ... ✓  8 genomes loaded
  [2/6] Counting k-mers & building PA matrix ... ✓  14,875 unique 19-mers, 264 conserved
  [3/6] Extracting unitigs ... ✓  103 unitigs (19-29bp)
  [4/6] Mapping to genomes ... ✓  447 hits across 8 genomes
  [5/6] Scoring windows ... ✓  60 windows above threshold
  [6/6] Designing primers ... ✓  3 candidate sets
  [7/7] Validating with MFEprimer ... ✓  coverage: 62.5%, 62.5%, 62.5%

  ✓ Pipeline complete in 0.9s
```

## Example Results

See [`tests/example_output/`](tests/example_output/) for a complete set of output files from 8 synthetic test genomes.

### Conservation Landscape

Sliding window scores across the genome — the peak at ~1000-1300bp is the discovered conserved region:

![Conservation Landscape](docs/conservation_landscape.png)

### MSA Conservation Heatmap

Per-position conservation in the extracted 600bp design region (green = conserved, red = variable):

![MSA Heatmap](docs/msa_heatmap.png)

## Commands

### `skipalign run` — Full pipeline

```bash
skipalign run \
  --input genomes/ \          # Directory with genome FASTA files (required)
  --annotations gff3/ \       # GFF3 annotations (optional)
  --k 19 \                    # k-mer length (default: 19)
  --min-genomes 3 \           # Min genomes for unitig filtering (default: 3)
  --window 300 \              # Scoring window size in bp (default: 300)
  --design-region 600 \       # Conserved region extraction size (default: 600)
  --top 5 \                   # Number of primer candidates (default: 5)
  --output results/           # Output directory (default: results/)
```

### `skipalign find-k` — Discover optimal k-mer length

```bash
skipalign find-k \
  --input genomes/ \
  --k-min 9 --k-max 51 \
  --output acf_result.tsv
```

Use this when targeting a new viral genus where the optimal k is unknown.

## Output Files

| File | Description |
|------|-------------|
| `report.html` | Self-contained HTML report with conservation plots, MSA heatmap, primer candidates, and MFEprimer validation results |
| `primers.tsv` | Primer-probe candidates with Tm, GC%, degeneracy, and TaqMan rule validation |
| `conserved_region.fasta` | Extracted conserved region sequences from all genomes |
| `msa_alignment.fasta` | MAFFT alignment of the conserved region |
| `pipeline_summary.json` | Pipeline statistics and parameters |
| `validation_summary.txt` | MFEprimer in-silico PCR coverage results |

## Default Parameters

| Parameter | Default | Source |
|-----------|---------|--------|
| k-mer length | 19 | ACF analysis of 51 Orthoflavivirus genomes |
| Window size | 300 bp | ~2-4x typical TaqMan amplicon (70-150bp) |
| Design region | 600 bp | Flanking context for primer placement |
| Min genomes | 3 | High-confidence unitig threshold |
| Degeneracy cap | 100 | Per-oligonucleotide variant limit |
| Consensus threshold | 50% | Per-position base frequency |

## External Dependencies

| Tool | Required | Purpose | Install |
|------|----------|---------|---------|
| **MAFFT** | Yes (for primer design) | Local MSA of conserved region | `conda install -c bioconda mafft` |
| **MFEprimer** | No (optional) | In-silico PCR validation | [GitHub releases](https://github.com/quwubin/MFEprimer-3.0/releases) |

If MAFFT is not installed, the pipeline will still discover conserved regions but skip primer design. If MFEprimer is not installed, validation is skipped.

## Citation

If you use skipalign, please cite **both the software and the original method paper**:

**Software:**
> skipalign: Alignment-free conserved region discovery and RT-qPCR primer-probe design (2026). https://github.com/Key-man-fromArchive/skipalign

**Method:**
> Sayasit K, Chaimayo C, Nuwong W, et al. Alignment-Free Guided Design of a Pan-Orthoflavivirus RT-qPCR Assay. *bioRxiv* (2026). https://doi.org/10.64898/2026.03.17.712358

GitHub's "Cite this repository" button (sidebar) provides BibTeX and APA formats via [`CITATION.cff`](CITATION.cff).

## License

MIT
