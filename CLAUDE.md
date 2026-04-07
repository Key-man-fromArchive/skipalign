# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Implementation of the **Alignment-Free Guided Design Pipeline** for pan-Orthoflavivirus RT-qPCR assay development, based on [Sayasit et al. (2026) bioRxiv](https://doi.org/10.64898/2026.03.17.712358). The pipeline replaces traditional multiple sequence alignment (MSA) with k-mer analysis and compacted De Bruijn graphs to discover conserved diagnostic targets across highly divergent viral genomes.

Reference paper: `docs/paper/Aligment-free-guided.pdf`

## Pipeline Architecture (7 stages)

The pipeline executes sequentially; each stage produces artifacts consumed by downstream stages:

```
Stage 1: Data Acquisition & Preprocessing
  NCBI RefSeq download (18,678 accessions) → filter by keywords ("complete genome", etc.)
  → 11,846 genomes (10,472 non-segmented + 1,374 segmented)
  → filename convention: <taxid>.fasta or <taxid_accession>.fasta
  Taxonomic annotation via TaxonKit v0.20.0

Stage 2: AF Phylogeny (Mash-style Jaccard)
  KITSUNE dmatrix module, k=11 fixed, Jaccard → Mash transform D = (-1/k)ln(2J/(1+J))
  → pairwise distance matrix → neighbor-joining tree (Pestivirus outgroup)
  Purpose: clade visualization/validation only, NOT inferential phylogenetics

Stage 3: ACF-guided k Selection (KITSUNE)
  KITSUNE ACF with --fast option (Jellyfish, no frequency filter)
  Each Flaviviridae genome queried against full RefSeq reference set
  Subset to Orthoflavivirus → plot ACF vs k (k=9..51)
  Key finding: k=19 is shortest k where Orthoflavivirus has specific common k-mers
  Output: 526,787 unique 19-mers across 51 orthoflaviviral genomes

Stage 4: Compacted De Bruijn Graph → Unitigs
  Input: PA (presence-absence) matrix of 19-mers × 51 genomes
  Bifrost for cDBG construction → 4,814 unitigs total
  High-confidence filter: unitigs present in ≥3 genomes → 339 unitigs (19–41 bp)

Stage 5: Exact Mapping & Feature Annotation
  Bowtie v1.0.0 (exhaustive exact-match, all zero-mismatch hits)
  SAMtools v1.16.1 (SAM→BAM) → BEDtools v2.30.0 (BAM→BED)
  BEDtools intersect with GFF3 annotations → functional context per unitig

Stage 6: Conserved Window Selection & Scoring
  Sliding 300bp windows across features, scored by:
    (i)  total unitig bp coverage
    (ii) number of genomes with unitig hits
  Top windows at cutoff ≥15 genomes → combine adjacent windows → 600bp design region
  Key result: NS5 gene, positions ~8.7–9.3 kbp

Stage 7: Sequence Extraction & Primer-Probe Design
  Extract 600bp from each genome → MAFFT local MSA
  Manual primer-probe design rules (TaqMan):
    - Forward primer: conserved 3' end, ≤100 degenerate variants, no >4 consecutive G
    - Probe: no 5' G, no >4 G-run, no ≥6 A-run, Tm > primer Tm
    - IUPAC ambiguity at ≥50% frequency threshold
  In silico validation: BLASTN vs NCBI nr (≤2 mismatches)
  Final: 23nt F-primer, 24nt R-primer, 25nt probe (Cy5/BHQ2), 135bp amplicon
```

## External Tools & Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| KITSUNE | - | ACF calculation, k-mer phylogeny, dmatrix |
| Jellyfish | 2.3.1 | k-mer counting |
| Bifrost | - | Compacted De Bruijn graph, unitig extraction |
| Bowtie | 1.0.0 | Exact-match short sequence mapping |
| SAMtools | 1.16.1 | SAM/BAM manipulation |
| BEDtools | 2.30.0 | Genomic interval operations, GFF3 intersection |
| MAFFT | - | Local MSA of extracted conserved regions |
| TaxonKit | 0.20.0 | NCBI taxonomy annotation |
| Python 3 | ≥3.10 | Pipeline orchestration, analysis, visualization |

## Key Biological Parameters

- **k-mer length for genus-level specificity**: k=19 (ACF-derived)
- **k-mer length for clade visualization**: k=11 (KITSUNE benchmark)
- **High-confidence unitig threshold**: present in ≥3 of 51 genomes
- **Conserved window size**: 300bp (scoring), 600bp (design extraction)
- **Genome cutoff for conservation**: ≥15 of 51 genomes
- **Target region**: NS5 gene (~8.7–9.3 kbp on Orthoflavivirus genome)
- **Degeneracy cap**: ≤100 variants per oligonucleotide
- **Consensus base threshold**: ≥50% per-position frequency
- **Max mismatches for BLAST specificity check**: ≤2

## Design Decisions & Rationale

- **NS5 over UTRs**: UTRs show high conservation but contain sfRNA artifacts that inflate copy numbers; NS5 provides stoichiometrically accurate viral genome quantification.
- **k=19 selection**: Shortest k-mer length where Orthoflavivirus genus retains specific common k-mers across all 51 genomes — balances specificity against sensitivity.
- **300bp window (not amplicon size)**: 2-4x typical TaqMan amplicon (70–150bp); gives room for primer/probe placement while maintaining fine genomic resolution.
- **Compacted De Bruijn graph**: Merges overlapping k-mers into unitigs to eliminate trivial sliding-window duplicates; reduces 526,787 k-mers → 4,814 unitigs.
- **Exact matching only**: Bowtie in zero-mismatch mode ensures unitigs faithfully represent sequences actually present in genomes.

## NotebookLM

Project notebook alias: (to be configured)
