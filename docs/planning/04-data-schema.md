# skipalign 데이터 스키마

> CLI 도구이므로 데이터베이스 대신 **중간 파일 형식**을 정의합니다.

## 1. 입력 데이터

### FASTA (필수)
```
genomes/
├── DENV1.fasta      # 단일 레코드 (complete genome)
├── DENV2.fasta
├── ZIKV.fasta
└── JEV.fasta
```
- 파일명 stem = 게놈 식별자
- 파일당 첫 번째 레코드만 사용

### GFF3 (선택적)
```
annotations/
├── DENV1.gff3       # FASTA와 동일한 stem
├── DENV2.gff3
└── ...
```

## 2. 중간 산출물 (work directory)

### k-mer counts (`work/kmers/`)
```json
{
  "genome": "DENV1",
  "k": 19,
  "total_kmers": 10542,
  "unique_kmers": 9876,
  "kmers": ["ACGTACGTACGTACGTACG", "..."]
}
```

### Presence-absence matrix (`work/pa_matrix.npz`)
- scipy sparse CSR format
- rows: k-mers (index in `work/kmer_index.tsv`)
- columns: genomes (index in `work/genome_names.tsv`)

### Unitigs (`work/unitigs.tsv`)
| unitig_id | sequence | length | genome_count | genomes |
|-----------|----------|--------|-------------|---------|
| UTG_001 | ACGT... | 23 | 15 | DENV1,DENV2,... |

### Mapping hits (`work/hits.tsv`)
| unitig_id | genome | start | end | strand | feature |
|-----------|--------|-------|-----|--------|---------|
| UTG_001 | DENV1 | 8742 | 8765 | + | NS5 |

### Scored windows (`work/windows.tsv`)
| rank | start | end | coverage_score | genome_count | feature |
|------|-------|-----|---------------|-------------|---------|
| 1 | 8700 | 9000 | 4523 | 18 | NS5 |

## 3. 최종 출력 (results directory)

### primers.tsv
| primer_set_id | type | sequence | length | tm | gc_percent | degeneracy | issues |
|--------------|------|----------|--------|-----|-----------|------------|--------|
| PS_001 | forward | ACGT... | 23 | 56.1 | 48.2 | 8 | |
| PS_001 | reverse | TGCA... | 24 | 56.3 | 45.0 | 12 | |
| PS_001 | probe | GCTA... | 25 | 58.4 | 50.7 | 4 | |

### conserved_region.fasta
- 각 게놈에서 추출된 보존 영역 서열 (multi-FASTA)

### msa_alignment.fasta
- MAFFT 정렬 결과

### pipeline_summary.json
```json
{
  "version": "0.1.0",
  "params": { "k": 19, "window": 300, "min_genomes": 3 },
  "input": { "genomes": 51, "total_bp": 561000 },
  "results": {
    "unique_kmers": 526787,
    "unitigs": 4814,
    "high_confidence_unitigs": 339,
    "conserved_windows": 2,
    "conserved_region": { "gene": "NS5", "start": 8700, "end": 9300 },
    "primer_candidates": 5
  },
  "runtime_seconds": 120
}
```

### report.html
- 보존 영역 위치 시각화 (genome coordinate plot)
- MSA conservation heatmap
- Primer-probe 후보 테이블 (Tm, GC%, degeneracy)
- 파이프라인 실행 통계
