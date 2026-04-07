# skipalign PRD (Product Requirements Document)

## 1. 제품 개요

- **서비스명**: skipalign
- **한줄 설명**: Alignment-free conserved region discovery and RT-qPCR primer-probe design
- **핵심 가치**: MSA가 실패하는 고도 분기 바이러스 속(>25% divergence)에서도 보존 진단 타겟을 발굴
- **타겟 사용자**: 임상 진단 개발자 (RT-qPCR 키트 개발)
- **라이선스**: MIT

## 2. 핵심 기능

### 기능 1: 보존 영역 발굴 (Alignment-Free Discovery)
- **설명**: 입력 게놈을 k-mer로 분해 → presence-absence matrix → compacted De Bruijn graph → unitig 추출 → sliding window scoring으로 보존 영역 발굴
- **Why**: 전통적 MSA는 25-30% 이상 divergence에서 false conservation을 생성. k-mer 기반 접근은 정렬 편향 없이 진짜 보존 영역을 찾음
- **우선순위**: MVP 필수

### 기능 2: Primer-Probe 자동 설계 (TaqMan Design)
- **설명**: 발굴된 보존 영역(600bp)에서 MAFFT local MSA → TaqMan 설계 규칙 적용 → primer3로 후보 생성
- **Why**: 보존 영역을 찾은 뒤에도 degeneracy, Tm, 구조 등 복잡한 규칙 적용이 필요. 자동화하여 반복 실험을 줄임
- **우선순위**: MVP 필수

### 기능 3: In-silico 검증 및 리포트 (Validation & Report)
- **설명**: MFEprimer로 in-silico PCR 검증 + HTML 시각 리포트(보존 영역 맵, MSA heatmap, primer 후보 테이블) + TSV 데이터
- **Why**: 임상 진단 개발자는 후보 primer의 신뢰성을 빠르게 판단해야 함. 시각적 리포트와 기계 파싱용 데이터를 동시에 제공
- **우선순위**: MVP 필수

## 3. 사용자 스토리

- As a **진단 키트 개발자**, I want to **입력 폴더에 바이러스 게놈 FASTA를 넣고 한 커맨드로 실행**, so that **보존 영역과 primer 후보를 자동으로 얻을 수 있다**
- As a **바이러스학 연구자**, I want to **새로운 바이러스 속에도 동일한 파이프라인을 적용**, so that **특정 바이러스에 종속되지 않는 범용 진단 설계가 가능하다**
- As a **감시 기관 담당자**, I want to **HTML 리포트로 결과를 빠르게 검토**, so that **primer 후보의 coverage와 특이성을 한눈에 판단할 수 있다**

## 4. 제약 조건

- **외부 바이너리**: MAFFT 1개 + MFEprimer 1개 (선택적)
- **하드웨어**: 서버/HPC (64GB+ RAM) 타겟, 수백 게놈까지 지원
- **배포**: pip / Conda / Docker 3-tier
- **입력**: FASTA 필수, GFF3 선택적

## 5. MVP 범위

Full pipeline + 리포트:
1. k-mer counting → PA matrix → unitig extraction
2. Genome mapping → window scoring → conserved region extraction
3. MAFFT MSA → primer-probe design (primer3 + TaqMan rules)
4. MFEprimer in-silico validation (선택적)
5. HTML + TSV 결과 리포트

## 6. 성공 지표

- 논문 재현: Orthoflavivirus 51개 게놈에서 NS5 영역 동일 발굴
- 범용성: 다른 바이러스 속(예: Alphavirus)에도 적용 가능
- 사용성: `skipalign run -i genomes/ -o results/` 한 줄로 end-to-end 실행
