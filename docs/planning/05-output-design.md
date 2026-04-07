# skipalign 출력 설계

> CLI 도구이므로 UI 디자인 대신 **터미널 출력 + HTML 리포트 설계**를 정의합니다.

## 1. 터미널 출력 (Rich)

### Progress 표시
```
skipalign v0.1.0 — Alignment-free primer design pipeline

  Input: 51 genomes from genomes/
  Parameters: k=19, window=300bp, min_genomes=3

  [1/6] Counting k-mers ...................... ✓  526,787 unique 19-mers
  [2/6] Building PA matrix ................... ✓  339 conserved k-mers
  [3/6] Extracting unitigs ................... ✓  339 unitigs (19-41bp)
  [4/6] Mapping to genomes ................... ✓  1,247 hits across 51 genomes
  [5/6] Scoring windows ...................... ✓  2 windows above threshold
  [6/6] Designing primers .................... ✓  5 candidate sets

  ✓ Pipeline complete in 1m 42s
  Results: results/report.html
```

### 경고/에러 스타일
```
  ⚠ MAFFT not found — install via: conda install -c bioconda mafft
  ⚠ No conserved windows found at threshold ≥15 — try lowering --min-genomes
  ✗ No FASTA files found in genomes/ — provide .fasta or .fa files
```

## 2. HTML 리포트 구조

### 헤더
- 프로젝트 제목, 실행 시간, 파라미터 요약

### 섹션 1: Pipeline Summary
- 입력 게놈 수, 총 bp, k 값
- 각 단계별 요약 통계 (테이블)

### 섹션 2: Conservation Landscape
- X축: pseudo-genomic coordinate (0-11kb)
- Y축: genome count (0-51)
- 각 300bp window의 score를 bar chart로 표시
- 선택된 보존 영역을 하이라이트

### 섹션 3: MSA Visualization
- 보존 영역 600bp MSA heatmap
- 색상: conservation level (green=conserved, red=variable)
- Primer/probe binding sites를 박스로 표시

### 섹션 4: Primer-Probe Candidates
- 후보별 테이블: sequence, length, Tm, GC%, degeneracy, amplicon size
- TaqMan 규칙 위반 시 경고 표시
- MFEprimer 결과 (있을 경우)

### 섹션 5: Validation (선택적)
- MFEprimer in-silico PCR 결과
- 특이성 히트 수, cross-reactivity 경고

## 3. TSV 출력 설계

모든 TSV는 header row 포함, tab-delimited, UTF-8 인코딩.

| 파일 | 용도 |
|------|------|
| `primers.tsv` | 기계 파싱용 primer 후보 전체 데이터 |
| `windows.tsv` | 모든 scored window (threshold 이하 포함) |
| `unitigs.tsv` | 고신뢰 unitig 목록 |
