# skipalign TRD (Technical Requirements Document)

## 1. 기술 스택

### Core
- **Language**: Python ≥3.10
- **CLI Framework**: Typer + Rich (터미널 출력)
- **Build System**: Hatchling (pyproject.toml)

### Bioinformatics
- **FASTA I/O**: BioPython (SeqIO)
- **k-mer counting**: Pure Python (canonical k-mer sets) / Jellyfish (optional, large datasets)
- **De Bruijn Graph**: Pure Python (networkx 또는 custom)
- **MSA**: MAFFT (외부 바이너리, 유일한 필수 의존성)
- **Primer Design**: primer3-py
- **In-silico PCR**: MFEprimer (외부 바이너리, 선택적)

### Data Processing
- **Matrix**: scipy.sparse (CSR) + numpy
- **Data manipulation**: 표준 라이브러리 (dataclasses, pathlib)

### Report
- **HTML 생성**: Jinja2 템플릿
- **시각화**: matplotlib 또는 plotly (conservation plot, MSA heatmap)
- **TSV 출력**: csv 표준 라이브러리

## 2. 아키텍처

```
CLI Layer (Typer)
    │
    ▼
Pipeline Orchestrator (pipeline.py)
    │
    ├── kmer.py        ── k-mer counting, canonical form
    ├── matrix.py      ── PA matrix construction, conservation filter
    ├── unitig.py      ── cDBG → unitig extraction
    ├── mapper.py      ── exact matching + GFF3 annotation
    ├── scorer.py      ── 300bp window scoring
    ├── primer.py      ── MAFFT + TaqMan rules + primer3
    ├── validator.py   ── MFEprimer in-silico PCR (optional)
    ├── reporter.py    ── HTML + TSV report generation
    └── io.py          ── FASTA/GFF3 parsing
```

- **패턴**: Pipeline (sequential stages with intermediate artifacts)
- **데이터 흐름**: 각 단계의 출력이 다음 단계의 입력 (functional pipeline)
- **상태 관리**: Stateless — 모든 중간 결과를 work directory에 파일로 저장

## 3. 성능 요구사항

| 시나리오 | 게놈 수 | k | 예상 런타임 | 메모리 |
|----------|---------|---|-------------|--------|
| Small (논문 재현) | 51 | 19 | < 5분 | < 1GB |
| Medium | 200 | 19 | < 30분 | < 8GB |
| Large | 1000 | 19 | < 2시간 | < 32GB |

## 4. 배포

### Tier 1: pip
```bash
pip install skipalign
# MAFFT를 별도 설치 필요
```

### Tier 2: Conda
```bash
conda env create -f environment.yml
# MAFFT 포함 자동 설치
```

### Tier 3: Docker
```bash
docker run -v ./data:/data skipalign run -i /data/genomes -o /data/results
# 모든 의존성 번들
```

## 5. 개발 환경

- **Linter**: ruff
- **Test**: pytest + pytest-cov
- **CI**: GitHub Actions (lint + test + Docker build)
- **Python**: 3.10, 3.11, 3.12 매트릭스 테스트
