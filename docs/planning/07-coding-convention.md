# skipalign 코딩 컨벤션

## 1. 프로젝트 구조

```
skipalignment-ppd/
├── src/skipalign/       # 소스 코드
│   ├── cli.py           # Typer CLI entry point
│   ├── pipeline.py      # Full pipeline orchestration
│   ├── kmer.py          # k-mer counting, canonical form
│   ├── matrix.py        # PA matrix construction
│   ├── unitig.py        # cDBG → unitig extraction
│   ├── mapper.py        # Exact matching + GFF3
│   ├── scorer.py        # Window scoring
│   ├── primer.py        # TaqMan design + primer3
│   ├── validator.py     # MFEprimer in-silico PCR
│   ├── reporter.py      # HTML + TSV report
│   └── io.py            # FASTA/GFF3 I/O
├── tests/               # pytest 테스트
│   ├── data/            # 테스트 데이터 (소형 게놈 세트)
│   └── test_*.py
├── docs/planning/       # 기획 문서
├── pyproject.toml       # 빌드 설정, 의존성
├── Dockerfile           # Docker 이미지
├── environment.yml      # Conda 환경
└── README.md
```

## 2. 네이밍 규칙

| 대상 | 규칙 | 예시 |
|------|------|------|
| 모듈 | snake_case | `kmer.py`, `pa_matrix.py` |
| 함수 | snake_case | `count_kmers()`, `build_pa_matrix()` |
| 클래스 | PascalCase | `ScoredWindow`, `PrimerProbeSet` |
| 상수 | UPPER_SNAKE_CASE | `IUPAC_EXPAND`, `MAX_DEGENERACY` |
| 변수 | snake_case | `genome_count`, `kmer_index` |
| dataclass 필드 | snake_case | `coverage_score`, `amplicon_length` |

## 3. 타입 힌트

- 모든 public 함수에 타입 힌트 필수
- `from __future__ import annotations` 사용 (PEP 604 `X | Y` 문법)
- 복잡한 타입은 `TypeAlias` 사용

```python
from __future__ import annotations
from scipy.sparse import csr_matrix

def build_pa_matrix(kmer_sets: dict[str, set[str]]) -> tuple[csr_matrix, list[str], list[str]]:
    ...
```

## 4. Docstring

- Google style docstring
- public 함수에만 작성 (private 함수는 코드가 self-documenting이면 생략)

```python
def count_kmers(sequence: str, k: int) -> set[str]:
    """Extract all canonical k-mers from a sequence.

    Args:
        sequence: Input nucleotide sequence (ACGT).
        k: k-mer length.

    Returns:
        Set of canonical k-mer strings.
    """
```

## 5. Linter / Formatter

- **ruff**: lint + format (pyproject.toml에 설정)
- `ruff check .` (lint)
- `ruff format .` (format)
- Line length: 100

```toml
[tool.ruff]
target-version = "py310"
line-length = 100

[tool.ruff.lint]
select = ["E", "F", "I", "N", "W"]
```

## 6. 테스트

- **pytest** 사용
- 파일 구조: `tests/test_{module}.py`
- 테스트 함수: `test_{기능}_{시나리오}()`
- 테스트 데이터: `tests/data/` (소형 FASTA 세트)

```python
def test_count_kmers_skips_n():
    """k-mers spanning N characters should be excluded."""
    kmers = count_kmers("ACGTNACGT", 4)
    assert all("N" not in k for k in kmers)
```

## 7. Git 커밋 메시지

Conventional Commits:

```
feat: add MFEprimer in-silico validation
fix: handle empty FASTA files gracefully
docs: update README with Docker usage
test: add integration test for full pipeline
refactor: extract window scoring into separate module
```

## 8. 에러 처리

- 사용자 입력 오류: `typer.echo()` + `raise typer.Exit(1)`
- 외부 도구 실패: `RuntimeError` with 설치 안내 메시지
- 내부 로직: 표준 Python exceptions (`ValueError`, `FileNotFoundError`)
- 생물학적 경고 (primer 규칙 위반 등): 경고 출력 후 계속 진행

## 9. 생물학적 규칙 상수

`primer.py`에 논문 기반 상수를 모아 정의:

```python
MAX_DEGENERACY = 100           # per-oligonucleotide variant cap
CONSENSUS_THRESHOLD = 0.50     # minimum per-position base frequency
MAX_CONSECUTIVE_G = 4          # forward primer / probe
MAX_CONSECUTIVE_A_PROBE = 5    # probe (>=6 prohibited)
DEFAULT_K = 19                 # ACF-derived for Orthoflavivirus
DEFAULT_WINDOW = 300           # scoring window (bp)
DEFAULT_DESIGN_REGION = 600    # extraction region (bp)
```
