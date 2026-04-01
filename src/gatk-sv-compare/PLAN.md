# `gatk-sv-compare` — Implementation Plan

## 1. Overview

A standalone Python CLI package for comparing two GATK-SV VCFs. It replaces ~3,500 lines of
R/bash spread across 8+ scripts and a multi-task WDL orchestration (`MainVcfQc.wdl`,
`CollectSiteLevelBenchmarking`, `PlotQcVcfWide`, `plot_callset_comparison.R`,
`plot_sv_vcf_distribs.R`, `plot_sv_perSample_distribs.R`, etc.) with a single pip-installable
Python package.

The tool has two subcommands:

1. **`preprocess`** — Runs GATK `SVConcordance` (twice, swapping eval/truth) and
   `SVRegionOverlap` to produce two fully annotated VCFs.
2. **`analyze`** — Streams the annotated VCFs and produces all tables, plots, and summary
   statistics.

> **This document is the build plan.** Each section is designed to be implemented and tested
> independently, in order. Later sections depend on earlier ones, but within a section the
> modules are independent of each other.


---
## 2. Guiding Principles

### 2.1 Modularity
Every analysis module is a self-contained unit with a single public entry function
`run(data: AggregatedData, config: AnalysisConfig, output_dir: Path) -> None`. Modules share
no mutable state. All shared logic lives in utility layers (`dimensions`, `vcf_reader`,
`plot_utils`).

### 2.2 DRY / Single Source of Truth
- Bucketing logic (AF bins, size bins, SV type classification, genomic context mapping) is
  defined **once** in `dimensions.py`.
- Plot styling (colors, fonts, figure sizes, SV type → color map) is defined **once** in
  `plot_utils.py`.
- VCF field extraction (SVTYPE, SVLEN, AF, STATUS, GQ, genotype counts) is defined **once** in
  `vcf_reader.py`.

### 2.3 Streaming & Memory Discipline
- Never load an entire VCF into memory. Stream with pysam `VariantFile`, one chromosome at a
  time.
- Build DataFrames from pre-allocated lists of dicts, then construct the DataFrame in one call.
  Never use `DataFrame.append()` or `pd.concat()` in a loop.
- For genotype-level operations on large cohorts (100K samples), extract only the fields needed
  (GQ, GT as allele count) per variant using numpy vectorization over the pysam genotype
  arrays, not Python loops over samples.

### 2.4 Parallelism
- Chromosome-level parallelism via `multiprocessing.Pool` (or `concurrent.futures.ProcessPoolExecutor`)
  for the aggregation pass.
- Module-level parallelism: analysis modules are independent and can run concurrently.
- Plotting: each figure is independent; parallelize figure generation with a thread pool (matplotlib
  is GIL-bound for rendering but I/O-bound for writing PNGs).

### 2.5 Testing
- Every module gets a companion `test_<module>.py`.
- Use small synthetic VCFs (5 variants, 3 samples) for unit tests — generated
  programmatically with pysam, not checked-in files.
- Integration tests use a small real-ish VCF fixture (50 variants, 10 samples).
- Pytest fixtures for `AggregatedData` objects with known expected outputs.

### 2.6 Type Safety
- All public functions have full type annotations.
- Use `dataclasses` or `NamedTuple` for structured data (not raw dicts).
- Run `mypy --strict` in CI.


---
## 3. VCF Format Specification & `validate` Subcommand

### 3.1 The Problem

GATK-SV produces VCFs in two internal styles that are structurally incompatible in subtle ways
(see `format_gatk_vcf_for_svtk.py` and `format_svtk_vcf_for_gatk.py`):

The table below is derived from examining real VCFs at each pipeline stage (see
`/Users/markw/Work/talkowski/sv-pipe-testing/mw_compare_vcfs/vcf_examples/`).

| Property | Final annotated (`all_batches.annotated`) | Cleaned (`all_batches.cleaned`) | Regenotyped (`chr1.regenotyped`) | Genotyped PESR / Depth | Cluster batch |
|---|---|---|---|---|---|
| REF | `N` | `N` | `N` | `N` | `N` |
| ALT (biallelic) | `<DEL>`, `<DUP>`, `<INS>`, `<INS:ME:ALU>`, `<INS:ME:LINE1>`, `<INS:ME:SVA>`, `<INS:ME>`, `<INV>`, `<BND>`, `<CPX>`, `<CTX>` | same | `<DEL>`, `<INS:UNK>`, `<INS:ME:*>`, `<INV>`, `<BND>`, `<CPX>`, `<DUP>` | same as cleaned | same |
| ALT (multiallelic CNV) | `<CNV>` (single symbolic, GT=`.`) | `<CNV>` | `<CN0>,<CN1>,<CN2>,<CN3>` | `<CNV>` | — |
| SVTYPE values observed | DEL(24432), INS(20938), DUP(11093), BND(9063), CPX(769), CNV(375), INV(100), CTX(1) | same | same minus CNV (uses DUP) | same minus CNV | no CNV |
| FILTER values | `PASS`, `UNRESOLVED`, `HIGH_NCR`, `HIGH_ALGORITHM_FDR`, `MULTIALLELIC`, combinations | `PASS`, `UNRESOLVED`, `MULTIALLELIC` | `PASS` only | `.` (empty) | `.` (empty) |
| ECN FORMAT | ✅ Present | ✅ Present | ❌ Absent | ❌ Absent | ❌ Absent |
| GQ range | 0–99 | 0–99 | 0–99 | 0–999 | — |
| OGQ FORMAT | ✅ (original pre-recalibration GQ) | ❌ | ❌ | ❌ | ❌ |
| SL FORMAT | ✅ (log-odds genotype score) | ❌ | ❌ | ❌ | ❌ |
| END for INS | `END == POS` (confirmed) | `END == POS` | `END` may differ from POS | `END` may differ from POS | `END == POS` |
| END for BND | `END == POS`, partner in `END2`+`CHR2` | same | same | END can be < POS (pysam warns) | END may be partner coord |
| END for CPX (dDUP) | `END == POS` (dispersed dup) | same | same | — | — |
| END for CPX (other) | `END > POS` (span of complex event) | same | same | — | — |
| SVLEN | Always present, positive, `Number=1` | same | same | Present, positive | Present |
| CHR2 | Present on BND/CTX only | Present on ALL records | Present on ALL | Present on ALL | Present on ALL |
| STRANDS | Optional (not always present) | Present on all | Present on all | Present on all | Present on all |
| Pre-computed geno counts | `N_BI_GENOS`, `N_HOMREF`, `N_HET`, `N_HOMALT`, `FREQ_*` | ❌ | ❌ | ❌ | ❌ |
| Pre-computed AF stats | `AC`, `AF`, `AN` (biallelic only) | ❌ | ❌ | ❌ | ❌ |
| Sex-stratified stats | `*_MALE`, `*_FEMALE` variants of all above | ❌ | ❌ | ❌ | ❌ |
| CNV copy-state stats | `CN_NUMBER`, `CN_COUNT`, `CN_FREQ`, `CN_NONREF_*` | ❌ | ❌ | ❌ | ❌ |
| Concordance truth fields | `TRUTH_DISTANCE_START/END`, `TRUTH_RECIPROCAL_OVERLAP`, `TRUTH_SIZE_SIMILARITY` | ❌ | ❌ | ❌ | ❌ |
| gnomAD AF annotations | `gnomad_v4.1_sv_AF`, `gnomad_v4.1_sv_*_AF` | ❌ | ❌ | ❌ | ❌ |
| MEMBERS INFO | ❌ | ❌ | ✅ | ✅ | ✅ |
| varGQ INFO | ❌ | ❌ | ✅ | ✅ | ❌ |
| MULTIALLELIC INFO flag | ❌ (in FILTER) | ❌ (in FILTER) | ✅ (INFO flag) | ✅ (INFO flag) | ❌ |
| NCN/NCR INFO | ✅ | ❌ | ❌ | ❌ | ❌ |

### 3.2 Expected Input Format

`gatk-sv-compare` expects **GATK-style VCFs**, consistent with the final output of the
GATK-SV pipeline (post-`FilterGenotypes`, post-annotation). The canonical reference is
`all_batches.annotated.vcf.gz`. Specifically:

1. **Symbolic ALT alleles.** Biallelic sites use a single symbolic allele (`<DEL>`, `<DUP>`,
   `<INS>`, `<INS:ME:ALU>`, `<INS:ME:LINE1>`, `<INS:ME:SVA>`, `<INS:ME>`, `<INV>`, `<BND>`,
   `<CPX>`, `<CTX>`). **CNV sites use `<CNV>` with FILTER=MULTIALLELIC** — these have `GT=.`
   (no standard genotype) and use `CN`/`CNQ` FORMAT fields instead. The regenotyped
   intermediate uses `<CN0>,<CN1>,<CN2>,<CN3>` multiallelic encoding which is NOT expected.
2. **`SVTYPE`** INFO field present on every record. Observed values: DEL, DUP, INS, INV,
   BND, CPX, CTX, CNV.
3. **`SVLEN`** INFO field present and **positive** for all types except BND/CTX (where it is
   absent). For CPX with `dDUP` type, SVLEN is the insert length. For INS, SVLEN is the
   insertion length.
4. **`END` semantics** are type-dependent:
   - **DEL/DUP/INV**: `END > POS`, defines the reference span.
   - **INS**: `END == POS` (point insertion; length from SVLEN).
   - **BND/CTX**: `END == POS`; partner position in `END2`, partner chromosome in `CHR2`.
   - **CPX (dDUP, dDUP_iDEL)**: `END == POS` (dispersed duplication to a distant locus).
   - **CPX (others: delINV, INVdup, dupINVdup, etc.)**: `END > POS` (span of complex event).
5. **`ECN`** FORMAT field present (expected copy number per sample per contig).
6. **`GQ`** FORMAT field present, range **0–99**. Additionally, `OGQ` (original GQ before
   recalibration) and `SL` (log-odds score) are present in the final annotated VCF.
7. **No breakend notation** in ALT alleles — BND uses `<BND>` symbolic allele.
8. **FILTER** is populated. Observed values: `PASS`, `UNRESOLVED`, `HIGH_NCR`,
   `HIGH_ALGORITHM_FDR`, `MULTIALLELIC`, and combinations thereof.
9. **Pre-computed genotype count INFO fields** are present in the final annotated VCF:
   `N_BI_GENOS`, `N_HOMREF`, `N_HET`, `N_HOMALT`, `FREQ_HOMREF`, `FREQ_HET`,
   `FREQ_HOMALT`, plus sex-stratified versions (`*_MALE`, `*_FEMALE`).
   **These can be used directly** instead of re-computing genotype counts from the GT
   arrays, which is a major performance optimization for large cohorts.
10. **Pre-computed allele frequency fields**: `AC`, `AF`, `AN` for biallelic sites. For CNV
    sites, `AC=0`, `AF=.`, `AN=0` — use `CN_NONREF_FREQ` instead.
11. **CNV-specific fields**: `CN_NUMBER`, `CN_COUNT`, `CN_STATUS`, `CN_FREQ`,
    `CN_NONREF_COUNT`, `CN_NONREF_FREQ` provide copy-state distributions.
12. **No-call genotypes** (`./. ` or `.`) are common, especially on sites with `HIGH_NCR`
    filter. The `NCN` (no-call count) and `NCR` (no-call rate) INFO fields quantify this.

**Key optimization:** When pre-computed genotype count fields are present (`N_BI_GENOS`,
`N_HOMREF`, `N_HET`, `N_HOMALT`, `AC`, `AF`, `AN`), use them directly rather than
iterating over all sample genotypes. This avoids the O(n_samples) per-variant cost and
is critical for 100K-sample VCFs. The `vcf_reader` should detect these fields and use
them preferentially, falling back to GT-based counting only when absent.

### 3.3 `validate` Subcommand

```
gatk-sv-compare validate --vcf input.vcf.gz [--fix --out fixed.vcf.gz]
```

**Checks (report all issues, do not stop at first):**

| Check ID | Severity | Description |
|---|---|---|
| `MULTI_ALLELIC_NON_CNV` | ERROR | Non-CNV record has >2 alleles (e.g. `<CN0>,<CN1>,<CN2>,<CN3>` encoding from regenotyped stage) |
| `MISSING_SVTYPE` | ERROR | `SVTYPE` INFO field absent |
| `UNKNOWN_SVTYPE` | WARN | `SVTYPE` not in {DEL, DUP, INS, INV, BND, CTX, CPX, CNV} |
| `MISSING_SVLEN` | WARN | `SVLEN` absent for DEL/DUP/INS/INV/CPX |
| `SVLEN_SIGN` | WARN | `SVLEN` is negative (GATK-style expects positive for all types) |
| `SVLEN_NUMBER_DOT` | WARN | `SVLEN` header has `Number=.` instead of `Number=1` |
| `INS_END_MISMATCH` | WARN | INS record has `END != POS` (intermediate svtk-style has END at SR2POS) |
| `BND_END_MISMATCH` | WARN | BND/CTX record has `END != POS` and no `END2` |
| `CPX_DDDUP_END` | WARN | CPX dDUP/dDUP_iDEL record has `END != POS` |
| `BREAKEND_NOTATION` | ERROR | ALT allele uses breakend notation (`]` or `[`) |
| `MISSING_ECN` | WARN | `ECN` FORMAT field absent (present in cleaned+annotated, absent in genotyped/regenotyped) |
| `GQ_RANGE` | WARN | GQ values > 99 detected (genotyped-stage VCFs use 0–999 scale) |
| `EMPTY_FILTER` | WARN | FILTER is `.` (empty) — common in genotyped/cluster-batch stages |
| `MISSING_GT` | ERROR | `GT` FORMAT field absent |
| `CHR2_ON_NON_BND` | INFO | `CHR2` present on non-BND/CTX record (normal in intermediates, stripped in final) |
| `MEMBERS_PRESENT` | INFO | `MEMBERS` INFO field present (indicates pre-CleanVcf intermediate) |
| `VARGQ_PRESENT` | INFO | `varGQ` INFO field present (indicates pre-FilterGenotypes intermediate) |
| `MULTIALLELIC_INFO_FLAG` | INFO | `MULTIALLELIC` is an INFO flag rather than FILTER value (regenotyped-stage pattern) |
| `MISSING_PRECOMPUTED_COUNTS` | INFO | Pre-computed genotype count fields (`N_BI_GENOS`, `N_HOMREF`, etc.) absent — will fall back to GT-based counting (slower) |
| `CNV_NO_GT` | INFO | CNV record has `GT=.` — expected for final annotated CNVs, will use `CN_NONREF_FREQ` for frequency |

**`--fix` mode** (with `--out`): apply automatic corrections where possible:
- `SVLEN_SIGN`: take absolute value.
- `SVLEN_NUMBER_DOT`: rewrite header line.
- `GQ_RANGE`: scale down by dividing by 10 (floor), matching `format_svtk_vcf_for_gatk.py` logic.
- `EMPTY_FILTER`: set to `PASS`.
- `BREAKEND_NOTATION`: rewrite ALT to `<BND>`.
- `BND_END_MISMATCH`: move END to `END2`, set `END = POS`.
- `INS_END_MISMATCH`: set `END = POS` (SVLEN already carries the insert length).
- `CPX_DDDUP_END`: set `END = POS` for dDUP/dDUP_iDEL CPX types.
- `MULTI_ALLELIC_NON_CNV`: if ALT is `<CN0>,<CN1>,<CN2>,<CN3>` and SVTYPE is DUP,
  convert to biallelic `<DUP>` and rewrite GT from CN-based to diploid (matching
  `format_gatk_vcf_for_svtk.py` CN→GT logic, but in reverse).
- `CHR2_ON_NON_BND`: remove `CHR2` from non-BND/CTX records.

Unfixable errors (`MISSING_SVTYPE`, `MISSING_GT`) cause the tool to exit non-zero with a
diagnostic report.

**Pipeline stage detection:** The validate subcommand should auto-detect which pipeline
stage the VCF comes from by checking for sentinel fields:

| Stage | Sentinel Fields |
|---|---|
| Final annotated | `N_BI_GENOS`, `OGQ`, `SL`, gnomAD annotations |
| Cleaned (post-CleanVcf) | `ECN` present, no `N_BI_GENOS`, no `MEMBERS` |
| Regenotyped | `<CN0>,<CN1>,<CN2>,<CN3>` alleles, `varGQ`, no `ECN` |
| Genotyped (PESR/depth) | `varGQ`, `MEMBERS`, no `ECN`, GQ 0–999 |
| Cluster batch | `MEMBERS`, no `GQ`, no genotype annotations |

This detection is reported in the validation summary and used to tailor fix recommendations.

### 3.4 Implementation

```
gatk_sv_compare/
  validate.py          # ~200 lines
  vcf_format.py        # shared checks used by validate + vcf_reader
```

**Build order:** Implement `validate` first. It is the foundation — every downstream module
depends on the VCF being well-formed.


---
## 4. Package Structure

```
src/gatk-sv-compare/
├── PLAN.md                          # this file
├── README.md                        # existing design notes
├── pyproject.toml                   # build config, entry points, dependencies
├── src/
│   └── gatk_sv_compare/
│       ├── __init__.py              # package version
│       ├── cli.py                   # argparse CLI: preprocess, validate, analyze
│       │
│       ├── # ── Layer 0: Constants & Config ──
│       ├── dimensions.py            # bucketing constants, categorize_variant()
│       ├── config.py                # AnalysisConfig dataclass, I/O paths
│       │
│       ├── # ── Layer 1: VCF I/O ──
│       ├── vcf_format.py            # format checks, GATK vs svtk detection
│       ├── validate.py              # validate subcommand
│       ├── vcf_reader.py            # streaming pysam extraction → per-chrom DataFrames
│       │
│       ├── # ── Layer 2: Preprocessing ──
│       ├── preprocess.py            # SVConcordance + SVRegionOverlap orchestration
│       │
│       ├── # ── Layer 3: Aggregation ──
│       ├── aggregate.py             # AggregatedData dataclass, parallel chrom aggregation
│       │
│       ├── # ── Layer 4: Analysis Modules ──
│       ├── modules/
│       │   ├── __init__.py
│       │   ├── base.py              # AnalysisModule ABC
│       │   ├── binned_counts.py
│       │   ├── overall_counts.py
│       │   ├── genotype_dist.py
│       │   ├── counts_per_genome.py
│       │   ├── site_overlap_plots.py
│       │   ├── site_overlap_tables.py
│       │   ├── overlap_heatmaps.py
│       │   ├── allele_freq.py
│       │   ├── carrier_freq.py
│       │   ├── genotype_quality.py
│       │   ├── genotype_exact_match.py
│       │   └── genotype_concordance.py
│       │
│       ├── # ── Layer 5: Plotting Utilities ──
│       └── plot_utils.py            # colors, styles, stacked bars, scatter, ternary, heatmap
│
└── tests/
    ├── conftest.py                  # shared fixtures (synthetic VCFs, AggregatedData)
    ├── test_dimensions.py
    ├── test_vcf_format.py
    ├── test_validate.py
    ├── test_vcf_reader.py
    ├── test_aggregate.py
    ├── test_preprocess.py
    └── test_modules/
        ├── test_binned_counts.py
        ├── test_overall_counts.py
        ├── ...
        └── test_genotype_concordance.py
```


---
## 5. Layer 0: Constants & Configuration

### 5.1 `dimensions.py`

Single source of truth for all bucketing. Every other module imports from here.

```python
# --- SV Types ---
SVTYPES_CORE = ["DEL", "DUP", "INS", "INV", "BND", "CPX", "CTX", "CNV"]
SVTYPE_INS_MEI = "INS:MEI"
# Observed ALT alleles in real pipeline output:
#   <INS:ME:ALU> (6609), <INS:ME:LINE1> (2126), <INS:ME:SVA> (434), <INS:ME> (4)
# The <INS:UNK> allele in regenotyped intermediates is NOT an MEI.
MEI_ALT_PATTERNS = frozenset({"<INS:ME:", "<INS:ME>"})  # substring match on ALT allele
INS_NON_MEI_ALTS = frozenset({"<INS>", "<INS:UNK>"})    # explicitly not MEI

# --- AF Buckets ---
# Label, (min_exclusive, max_inclusive), special_rule
AF_BUCKETS = [
    ("AC=1",     None),           # special: AC == 1
    ("AC>1,AF<1%", (0.0, 0.01)), # AC > 1 and AF < 0.01
    ("1-10%",    (0.01, 0.10)),
    ("10-50%",   (0.10, 0.50)),
    (">50%",     (0.50, 1.00)),
]

# --- Size Buckets ---
# Applied to all SV types with a defined SVLEN; BND/CTX get "N/A"
SIZE_BUCKETS = [
    ("<100bp",       (0, 100)),
    ("100-500bp",    (100, 500)),
    ("500bp-2.5kb",  (500, 2500)),
    ("2.5-10kb",     (2500, 10000)),
    ("10-50kb",      (10000, 50000)),
    (">50kb",        (50000, float("inf"))),
]

# --- Genomic Context ---
# Populated from SVRegionOverlap INFO annotations
GENOMIC_CONTEXTS = ["segdup", "simple_repeat", "repeatmasker", "none"]

# --- FILTER handling ---
# Observed FILTER values in final annotated VCF:
#   PASS (53690), UNRESOLVED (9092), HIGH_NCR (2416), HIGH_ALGORITHM_FDR (879),
#   MULTIALLELIC (375), HIGH_ALGORITHM_FDR;HIGH_NCR (315), HIGH_NCR;UNRESOLVED (4)
# Default: analyze ALL variants. --pass-only flag restricts to PASS-only.
# CNV sites always have FILTER=MULTIALLELIC.
FILTER_PASS = "PASS"
FILTER_MULTIALLELIC = "MULTIALLELIC"

# --- SVConcordance STATUS values ---
STATUS_MATCHED = "MATCHED"
STATUS_UNMATCHED = "UNMATCHED"
```

**Key function:**

```python
@dataclass(frozen=True)
class VariantCategory:
    svtype: str          # e.g. "DEL", "INS:MEI", "BND"
    size_bucket: str     # e.g. "100-500bp", "N/A"
    af_bucket: str       # e.g. "AC=1"
    genomic_context: str # e.g. "segdup"

def categorize_variant(record_info: dict) -> VariantCategory:
    """Pure function. Takes a dict of extracted VCF fields, returns a VariantCategory.

    For CNV sites (svtype='CNV'), uses cn_nonref_freq as the AF for bucketing,
    since AC/AF/AN are uninformative (AC=0, AF='.', AN=0 in annotated VCF).
    For BND/CTX, size_bucket is always 'N/A' (no SVLEN).
    """
```

This function is called once per variant during the aggregation pass. It must be fast.

### 5.2 CNV Handling Strategy

CNV records (`SVTYPE=CNV`, `FILTER=MULTIALLELIC`) require special treatment throughout
the tool because they differ fundamentally from biallelic SV types:

| Aspect | Biallelic SV | CNV |
|---|---|---|
| ALT allele | `<DEL>`, `<DUP>`, etc. | `<CNV>` |
| GT field | `0/0`, `0/1`, `1/1`, `./.` | `.` (always missing) |
| Frequency metric | `AF` (from AC/AN) | `CN_NONREF_FREQ` |
| Sample count | `N_BI_GENOS` | `CN_NUMBER` |
| Genotype counts | `N_HOMREF`, `N_HET`, `N_HOMALT` | `CN_COUNT` (per copy state) |
| HWE analysis | Yes | **No** (not biallelic, HWE undefined) |
| GQ field | `GQ` | `GQ` (present but derived from `CNQ`) |
| Per-genome counting | Count from GT | Count from CN != ECN |

**Module-level handling:**
- `genotype_dist`: **Exclude** CNV from ternary plots and HWE tests.
- `counts_per_genome`: Count CNV carriers via `CN != ECN` rather than GT.
- `allele_freq`: Use `CN_NONREF_FREQ` as the frequency for AF correlation.
- `binned_counts` / `site_overlap_*`: Include CNV as its own SVTYPE category.
- `genotype_quality`: Use `GQ` (which is present on CNV records in the annotated VCF).
- `genotype_exact_match`: Compare `CN` values, not GT, for CNV sites.

### 5.3 `config.py`

```python
@dataclass
class AnalysisConfig:
    vcf_a_path: Path
    vcf_b_path: Path
    vcf_a_label: str            # human-readable label for VCF A
    vcf_b_label: str
    output_dir: Path
    reference_dict: Path
    contigs: list[str]          # from reference dict or contig list
    n_workers: int              # for multiprocessing
    modules: list[str] | None   # None = all; or subset of module names
    pass_only: bool             # if True, restrict to FILTER=PASS variants
    per_chrom: bool             # if True, produce per-contig stratified plots
    # Optional inputs
    seg_dup_track: Path | None
    simple_repeat_track: Path | None
    repeatmasker_track: Path | None
    gatk_path: str              # path to gatk executable
```


---
## 6. Layer 1: VCF I/O

### 6.1 `vcf_format.py`

Contains the low-level check functions used by both `validate.py` and `vcf_reader.py`.

```python
class FormatIssue:
    check_id: str      # e.g. "MISSING_SVTYPE"
    severity: str      # "ERROR", "WARN", or "INFO"
    record_id: str     # variant ID, or "" for header-level issues
    message: str

@unique
class PipelineStage(Enum):
    """Auto-detected pipeline stage based on header sentinel fields."""
    FINAL_ANNOTATED = "final_annotated"    # N_BI_GENOS, OGQ, SL, gnomAD fields
    CLEANED = "cleaned"                     # ECN, no N_BI_GENOS, no MEMBERS
    REGENOTYPED = "regenotyped"             # <CN0..CN3> alleles, varGQ, no ECN
    GENOTYPED = "genotyped"                 # varGQ, MEMBERS, no ECN, GQ 0-999
    CLUSTER_BATCH = "cluster_batch"         # MEMBERS, no GQ annotations
    UNKNOWN = "unknown"

def check_record(record: pysam.VariantRecord) -> list[FormatIssue]: ...
def detect_pipeline_stage(vcf: pysam.VariantFile) -> PipelineStage: ...
def is_mei(alt_allele: str) -> bool: ...
def is_cnv(record: pysam.VariantRecord) -> bool: ...
def has_precomputed_counts(vcf: pysam.VariantFile) -> bool: ...
```

### 6.2 `validate.py`

CLI-facing module. Streams the VCF, collects all `FormatIssue`s, prints a summary report,
and optionally writes a corrected VCF.

### 6.3 `vcf_reader.py`

The workhorse for the `analyze` subcommand. Provides a single streaming extractor:

```python
@dataclass
class SiteRecord:
    """Flat representation of one VCF record — only the fields we need."""
    variant_id: str
    contig: str
    pos: int
    end: int
    svtype: str             # already classified (incl. INS:MEI)
    svlen: int | None
    alt_allele: str         # raw ALT string, e.g. "<INS:ME:ALU>"
    filter: str             # FILTER value, e.g. "PASS", "HIGH_NCR"
    af: float
    ac: int
    an: int
    # Pre-computed genotype counts (from INFO fields if available, else from GT)
    n_bi_genos: int         # N_BI_GENOS or computed from GT
    n_hom_ref: int          # N_HOMREF
    n_het: int              # N_HET
    n_hom_alt: int          # N_HOMALT
    n_no_call: int          # NCN or computed
    # CNV-specific (only for SVTYPE=CNV)
    cn_nonref_freq: float | None    # CN_NONREF_FREQ
    # SVConcordance fields (only after preprocessing)
    status: str | None      # SVConcordance STATUS, or None if not annotated
    truth_vid: str | None   # matched truth variant ID
    # SVRegionOverlap fields (only after preprocessing)
    genomic_context: str    # from SVRegionOverlap, or "none"
    # Evidence and algorithm metadata
    algorithms: list[str] | None   # ALGORITHMS field
    evidence: list[str] | None     # EVIDENCE field
    # Optional per-sample data (expensive — only extracted when needed)
    gq_array: np.ndarray | None       # only if GQ modules requested
    concordance_metrics: dict | None   # VAR_PPV, VAR_SENSITIVITY, etc.

def iter_contig(
    vcf_path: Path,
    contig: str,
    extract_gq: bool = False,
    extract_concordance: bool = False,
    sample_indices: np.ndarray | None = None,  # subset of samples
) -> Iterator[SiteRecord]:
    """Yields SiteRecords for one contig. Zero-copy where possible."""
    ...
```

**Performance notes — critical for 100K-sample VCFs:**

- **Use pre-computed INFO fields first.** The final annotated VCF contains `N_BI_GENOS`,
  `N_HOMREF`, `N_HET`, `N_HOMALT`, `AC`, `AF`, `AN` as INFO fields. Reading these is
  O(1) per variant vs. O(n_samples) for GT iteration. The reader must detect these fields
  in the header and use them preferentially:
  ```python
  has_precomputed = "N_BI_GENOS" in vcf.header.info
  if has_precomputed:
      n_hom_ref = record.info["N_HOMREF"]
      ...
  else:
      # Fall back to GT-based counting
  ```
- **For CNV sites** (SVTYPE=CNV, FILTER=MULTIALLELIC), GT is `.` (missing). Use
  `CN_NONREF_FREQ` for frequency and `CN_NUMBER` for sample count. Do NOT try to
  parse GT for these records.
- **For GQ extraction**, allocate one numpy array per variant and fill it in one pass.
  Only extract GQ for non-CNV biallelic sites (CNV uses CNQ instead).
- **`sample_indices`** allows extracting genotypes for only the overlapping samples (for
  the genotype concordance modules), avoiding O(n_all_samples) when only a subset is needed.


---
## 7. Layer 2: Preprocessing (`preprocess.py`)

### 7.1 Workflow

```
Input: vcf_a, vcf_b, reference_dict, contig_list, tracks (optional)

Step 1: SVConcordance (A=eval, B=truth) → concordance_a.vcf.gz
Step 2: SVConcordance (B=eval, A=truth) → concordance_b.vcf.gz
Step 3: SVRegionOverlap(concordance_a) → annotated_a.vcf.gz   [if tracks provided]
Step 4: SVRegionOverlap(concordance_b) → annotated_b.vcf.gz   [if tracks provided]
```

Each step is scattered per-contig and run in parallel, then concatenated.

### 7.2 Implementation

```python
def run_preprocess(config: AnalysisConfig) -> tuple[Path, Path]:
    """
    Returns paths to the two fully annotated VCFs.
    All intermediate files go into config.output_dir / "preprocess/".
    """
    ...
```

Subprocess calls via `subprocess.run()` with proper error handling, logging, and
timeout. Each GATK invocation is wrapped in a helper that:
1. Constructs the command list.
2. Logs the full command.
3. Runs it and captures stdout/stderr.
4. Checks return code and raises a clear exception on failure.

### 7.3 Contig-Level Scatter

```python
def _run_sv_concordance_contig(
    eval_vcf: Path, truth_vcf: Path, contig: str, reference_dict: Path,
    output_path: Path, gatk_path: str,
    clustering_config: Path | None, stratification_config: Path | None,
    track_names: list[str] | None, track_intervals: list[Path] | None,
) -> Path:
    ...

def _run_sv_region_overlap_contig(
    vcf: Path, contig: str, reference_dict: Path, output_path: Path,
    gatk_path: str, tracks: dict[str, Path],
) -> Path:
    ...
```

These are `pickle`-safe top-level functions (not closures) so they work with
`multiprocessing.Pool.starmap()`.

### 7.4 Concatenation

Use pysam to concatenate per-contig VCFs (header from first, records streamed in order).
Index with tabix.


---
## 8. Layer 3: Aggregation (`aggregate.py`)

### 8.1 Overview

A single parallel pass over the two annotated VCFs that extracts all per-variant features
into chromosome-level DataFrames, then concatenates them.

### 8.2 `AggregatedData`

```python
@dataclass
class AggregatedData:
    """The output of the aggregation pass. Input to all analysis modules."""

    # Per-variant site-level data for each VCF
    # Columns: variant_id, contig, pos, svtype, svlen, af, ac, an, filter,
    #          status, truth_vid, genomic_context, size_bucket, af_bucket,
    #          n_hom_ref, n_het, n_hom_alt, n_no_call, is_cnv,
    #          cn_nonref_freq (NaN for non-CNV)
    # NOTE: For CNV sites (is_cnv=True), af is populated from cn_nonref_freq
    # rather than AC/AN, and n_hom_ref/n_het/n_hom_alt are NaN.
    # The --pass-only filter is applied here during aggregation (rows with
    # FILTER != PASS are dropped before DataFrame construction).
    sites_a: pd.DataFrame
    sites_b: pd.DataFrame

    # Per-variant GQ summaries (optional, only if GQ modules requested)
    # Stored as a dict of contig → numpy array (n_variants × n_samples)
    # to avoid materializing a giant DataFrame
    gq_a: dict[str, np.ndarray] | None
    gq_b: dict[str, np.ndarray] | None

    # Per-variant concordance metrics (optional, only if concordance modules requested)
    concordance_a: pd.DataFrame | None
    concordance_b: pd.DataFrame | None

    # Overlapping sample info
    shared_samples: list[str]
    sample_indices_a: np.ndarray | None  # indices of shared samples in VCF A
    sample_indices_b: np.ndarray | None

    # Labels
    label_a: str
    label_b: str
```

### 8.3 Aggregation Logic

```python
def aggregate(config: AnalysisConfig) -> AggregatedData:
    """
    Parallel aggregation over both VCFs.
    Returns AggregatedData ready for analysis modules.
    """
    # 1. Open both VCFs to read headers (sample lists, contigs)
    # 2. Determine shared samples
    # 3. Determine which optional extractions are needed based on config.modules
    # 4. Parallel map: for each contig, extract SiteRecords → list[dict]
    # 5. Concatenate per-contig lists → single DataFrame per VCF
    # 6. Add derived columns (size_bucket, af_bucket) via dimensions.categorize_variant
    # 7. Return AggregatedData
    ...
```

**Per-contig worker:**

```python
def _extract_contig(
    vcf_path: Path, contig: str, extract_gq: bool, extract_concordance: bool,
    sample_indices: np.ndarray | None,
) -> tuple[list[dict], np.ndarray | None]:
    """Returns (list of row dicts, optional GQ matrix) for one contig."""
    rows = []
    gq_arrays = []
    for site in iter_contig(vcf_path, contig, extract_gq, extract_concordance, sample_indices):
        rows.append(asdict(site))  # or a to_dict() method
        if extract_gq:
            gq_arrays.append(site.gq_array)
    gq_matrix = np.stack(gq_arrays) if gq_arrays else None
    return rows, gq_matrix
```

Building the DataFrame:

```python
# GOOD: one DataFrame construction from list of dicts
df = pd.DataFrame(rows)

# BAD: never do this
# df = pd.DataFrame()
# for row in rows:
#     df = df.append(row)  # O(n²)
```


---
## 9. Layer 4: Analysis Modules

All modules implement:

```python
class AnalysisModule(ABC):
    """Base class for all analysis modules."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Directory name for outputs, e.g. 'binned_counts'."""

    @property
    def requires_shared_samples(self) -> bool:
        """Override to True if the module only runs when samples overlap."""
        return False

    @abstractmethod
    def run(self, data: AggregatedData, config: AnalysisConfig) -> None:
        """Execute the analysis. Write outputs to config.output_dir / self.name /."""

    def output_dir(self, config: AnalysisConfig) -> Path:
        d = config.output_dir / self.name
        d.mkdir(parents=True, exist_ok=True)
        return d
```

### Module Registry

```python
ALL_MODULES: list[type[AnalysisModule]] = [
    BinnedCountsModule,
    OverallCountsModule,
    GenotypeDistModule,
    CountsPerGenomeModule,
    SiteOverlapPlotsModule,
    SiteOverlapTablesModule,
    OverlapHeatmapsModule,
    AlleleFreqModule,
    CarrierFreqModule,
    GenotypeQualityModule,
    GenotypeExactMatchModule,
    GenotypeConcordanceModule,
]
```

The CLI iterates the registry, checks `requires_shared_samples`, and calls `run()`.

---

### 9.1 `binned_counts` — Cross-Tabulated Count Tables

**Purpose:** Compressed tables of variant counts across all dimension combinations.

**Outputs:**
- `counts_a.tsv`, `counts_b.tsv` — one row per (svtype × size_bucket × af_bucket × genomic_context) cell,
  columns: `svtype`, `size_bucket`, `af_bucket`, `genomic_context`, `n_variants`, `n_matched`, `n_unmatched`.
- `counts_a.parquet`, `counts_b.parquet` — same, for programmatic downstream use.

**Implementation:**
```python
df.groupby(["svtype", "size_bucket", "af_bucket", "genomic_context"]).agg(
    n_variants=("variant_id", "count"),
    n_matched=("status", lambda s: (s == STATUS_MATCHED).sum()),
    n_unmatched=("status", lambda s: (s == STATUS_UNMATCHED).sum()),
).reset_index()
```

**Test:** Assert row counts sum to total variants. Assert known cell values on synthetic data.

---

### 9.2 `overall_counts` — Distribution Visualizations

**Purpose:** Bar charts of variant counts per VCF, sliced by single dimensions.

**Plots (per VCF, side by side):**
1. SV count by type (bar chart, one bar per SVTYPE, colored by type).
2. Size distribution by type (log-scale density/histogram, overlaid per type).
3. AF distribution by type (histogram, overlaid per type).
4. Genomic context by type (stacked bar, % in each context per type).

**Layout:** 2 × 4 panel figure (VCF A top row, VCF B bottom row), plus supporting individual PNGs.

**Implementation:** Thin wrappers around `plot_utils` functions. Data prep is pure pandas
`groupby` + `value_counts`.

---

### 9.3 `genotype_dist` — Hardy-Weinberg Ternary Plots

**Purpose:** Ternary plots showing (frac_hom_ref, frac_het, frac_hom_alt) per variant,
colored by HWE test result.

**Methodology** (matching `plot_sv_vcf_distribs.R` lines 1018–1090):
1. Restrict to **biallelic, autosomal** sites. **Exclude CNV** (SVTYPE=CNV) — they lack
   standard GT and HWE is undefined for multiallelic copy-number states.
2. Compute per-variant genotype fractions using **pre-computed INFO fields** (if available):
   `AA = N_HOMREF / N_BI_GENOS`, `AB = N_HET / N_BI_GENOS`, `BB = N_HOMALT / N_BI_GENOS`.
   Fall back to GT-based counting if pre-computed fields are absent.
3. Chi-squared HWE test per variant:
   - **Pass:** p ≥ 0.05 → dark green
   - **Nominal fail:** p < 0.05 but p ≥ 0.05/n_variants → light green
   - **Bonferroni fail:** p < 0.05/n_variants → magenta
4. Plot on equilateral triangle (Ref → Het → Hom vertices).
5. Draw HWE equilibrium curve (theoretical parabola) and Bonferroni boundary.

**Ternary projection** (no external dependency):
```python
def ternary_to_cartesian(aa: float, ab: float, bb: float) -> tuple[float, float]:
    """Convert (frac_ref, frac_het, frac_alt) to (x, y) in equilateral triangle."""
    x = 0.5 * (2 * bb + ab)
    y = (sqrt(3) / 2) * ab
    return x, y
```

**Panels:** All SV (both VCFs), then per size bucket: <100bp, 100–500bp, 500bp–2.5kb, 2.5–10kb,
10–50kb, >50kb.

**Additional plot:** Carrier frequency vs. allele frequency scatter (matching R line 1211).

---

### 9.4 `counts_per_genome` — Per-Sample SV Counts

**Purpose:** Box/violin plots of how many SV sites and alleles each sample carries.

**Computation:**
- For each sample, count sites where GT has at least one alt allele.
- For each sample, sum alt allele counts across all sites.
- **For CNV sites:** count carriers via `CN != ECN` (not GT, which is `.` for CNVs).
  The per-sample CN values must be read from FORMAT/CN, and ECN from FORMAT/ECN.
- Stratify by SV type, size bucket, AF bucket.

**Performance concern:** With 100K samples and 10M variants, we cannot iterate
sample-by-sample. Instead, during the aggregation pass, maintain per-sample running counts
using numpy vectorization over the GT arrays:

```python
# Per variant: gt_array shape = (n_samples,) with values 0, 1, 2 (alt allele count)
per_sample_sites += (gt_array > 0).astype(int)
per_sample_alleles += gt_array
```

These accumulators are stored in the `AggregatedData` or computed in a dedicated sub-pass.

**Plots:**
- Box/swarm plot: sites per genome by SV type.
- Box/swarm plot: alleles per genome by SV type.
- Heatmap grid: sites per genome by (size_bucket × SV class) and (af_bucket × SV class),
  showing median of N samples.
- Sites/alleles vs. GQ threshold curves (median across samples).

---

### 9.5 `site_overlap_plots` — Overlap Bar Charts

**Purpose:** Stacked bar charts showing overlap status between the two callsets.

**Overlap categorization** uses the `STATUS` INFO field from SVConcordance:
- `MATCHED` → overlapping site
- `UNMATCHED` → exclusive to this callset

For each VCF, produce stacked bars showing fraction matched vs. unmatched, sliced by:
1. SV type
2. Size bucket
3. AF bucket
4. Genomic context

**Color scheme:** Match existing pipeline (greens for matched, gray for unmatched).

**Layout:** 1 × 4 panel + legend (class, size, AF, context), for each VCF direction.

---

### 9.6 `site_overlap_tables` — Full Cross-Tabulated Overlap Metrics

**Purpose:** TSV/parquet tables with counts and percentages for every cell in the
svtype × size × AF × context grid.

**Columns:** `svtype`, `size_bucket`, `af_bucket`, `genomic_context`, `n_total_a`, `n_matched_a`,
`pct_matched_a`, `n_total_b`, `n_matched_b`, `pct_matched_b`.

---

### 9.7 `overlap_heatmaps` — Cross-Dimension Heatmap Grids

**Purpose:** The annotated heatmap grids shown in the example plots (Size × Freq, Size × Class,
Freq × Class) with count and percentage annotations per cell.

**Plots (per VCF):**
1. Size × Freq overlap heatmap
2. Size × Class overlap heatmap
3. Freq × Class overlap heatmap

**Each cell shows:** `n_matched / n_total (pct%)` with color intensity proportional to overlap
rate.

**Implementation:** `matplotlib.pyplot.imshow()` or `seaborn.heatmap()` with custom text
annotations. Cell color on a diverging colormap (e.g. YlGn).

---

### 9.8 `allele_freq` — AF Correlation Analysis

**Purpose:** For matched sites, scatter plot of AF in VCF A vs. AF in VCF B.

**Methodology** (matching `plot_callset_comparison.R` lines 305–400):
1. Join sites_a and sites_b on `(variant_id, truth_vid)` — matched pairs.
2. Scatter plot with:
   - Points colored by density.
   - Linear regression line (red).
   - Rolling mean (dark red).
   - Pearson R, Spearman ρ, p-value in title.
3. Identity line (gray diagonal) for reference.

**Plots:**
- Overall scatter.
- Per SV type.
- Per size bucket.
- Per AF bucket.
- Per genomic context.

**Tables:** For each (svtype × size × AF × context) cell: n_matched, R, ρ, mean absolute
AF difference.

---

### 9.9 `carrier_freq` — Carrier Frequency vs. Allele Frequency

**Purpose:** QC diagnostic scatter of carrier frequency (fraction of samples with ≥1 alt
allele) vs. allele frequency, per VCF. Deviation from the expected curve implies genotyping
bias (e.g. systematic het↔hom miscalls).

**Computation:**
- Carrier freq = `(n_het + n_hom_alt) / (n_het + n_hom_alt + n_hom_ref)`
- AF = `ac / an`

**Plot:** Scatter with identity line, rolling mean, and theoretical HWE relationship curve.
Restrict to biallelic autosomal sites.

---

### 9.10 `genotype_quality` — GQ Distributions

**Purpose:** Histograms and summary statistics of genotype quality.

**Plots:**
- Overall GQ histogram (alt genotypes only) per VCF, overlaid.
- Per SV type.
- Per size bucket.
- Per AF bucket.

**Tables:** Mean, median, Q25, Q75 GQ for each (svtype × size × AF × context) cell.

**Performance:** GQ arrays are large (n_variants × n_samples). Process per-contig, accumulate
histogram bin counts (not raw values). Use `np.histogram()` with fixed bins [0, 5, 10, ..., 95, 99].

---

### 9.11 `genotype_exact_match` — Genotype Concordance (Shared Samples Only)

**Requires:** `shared_samples` is non-empty.

**Purpose:** For matched sites, compare genotypes of shared samples between the two VCFs.

**Metrics per site:**
- **Exact match rate:** fraction of shared samples with identical GT.
- **Hom-ref ↔ Het mismatch rate.**
- **Hom-alt ↔ Hom-ref mismatch rate.**
- **Het ↔ Hom-ref mismatch rate.**

**Plots:** Histograms of each metric, stratified by svtype, size, AF, context (each
stratification as a separate panel row).

**Tables:** Mean/median metric for each (svtype × size × AF × context) cell.

---

### 9.12 `genotype_concordance` — SVConcordance INFO Metrics (Shared Samples Only)

**Requires:** `shared_samples` is non-empty and SVConcordance INFO fields present.

**Purpose:** Distribution plots of SVConcordance's per-variant concordance metrics.

**Fields:**
- `VAR_SENSITIVITY`
- `VAR_PPV`
- `VAR_SPECIFICITY`
- `HET_SENSITIVITY`
- `HET_PPV`
- `HOMVAR_SENSITIVITY`
- `HOMVAR_PPV`
- `GENOTYPE_CONCORDANCE`
- `NON_REF_GENOTYPE_CONCORDANCE`
- `CNV_CONCORDANCE`

**Plots:** Histogram of each field, stratified by svtype, size, AF, context.

**Tables:** Mean/median per cell.


---
## 10. Layer 5: Plotting Utilities (`plot_utils.py`)

### 10.1 Style Constants

```python
# SV type colors (consistent with existing pipeline)
SVTYPE_COLORS = {
    "DEL": "#D43925",
    "DUP": "#2376B2",
    "INS": "#7B2D8E",
    "INV": "#F57E20",
    "BND": "#4DAF4A",
    "CTX": "#4DAF4A",
    "CPX": "#E31A8B",
    "CNV": "#984EA3",
    "INS:MEI": "#AF5FA0",
    "OTH": "#999999",
}

# Overlap status colors
OVERLAP_COLORS = {
    "matched": "#4DAC26",
    "unmatched": "#696969",
}

# HWE colors
HWE_COLORS = {
    "pass": "#4DAC26",
    "nominal": "#81F850",
    "bonferroni": "#AC26A1",
}

FIGURE_DPI = 300
DEFAULT_FIGSIZE = (12, 8)  # inches
```

### 10.2 Reusable Plot Functions

```python
def plot_stacked_bars(ax, matrix, colors, labels, scaled=True, annotate_counts=True): ...
def plot_scatter_af(ax, x, y, label_x, label_y, show_lm=True, show_rolling=True): ...
def plot_ternary(ax, aa, ab, bb, colors, draw_hwe_curve=True, alpha=0.05): ...
def plot_heatmap_annotated(ax, matrix, row_labels, col_labels, fmt="{n}/{total}\n({pct}%)"): ...
def plot_histogram(ax, values, bins, color, label): ...
def plot_boxplot_grouped(ax, data_dict, colors): ...
def save_figure(fig, path: Path, dpi=FIGURE_DPI): ...
```

### 10.3 Multi-Panel Composite Figures

```python
def composite_overlap_panel(data: AggregatedData, config: AnalysisConfig, output_path: Path):
    """6-panel overlap figure: by-class + by-size + by-AF + AF-correlation + 2 heatmaps."""
    fig = plt.figure(figsize=(24, 6))
    gs = fig.add_gridspec(1, 6, width_ratios=[3, 4, 3, 3, 4, 4])
    ...

def composite_genotype_panel(data: AggregatedData, config: AnalysisConfig, output_path: Path):
    """8-panel figure: carrier-freq-vs-AF + ternary (all + 6 size buckets)."""
    fig = plt.figure(figsize=(24, 4))
    gs = fig.add_gridspec(1, 8, width_ratios=[2, 2, 1, 1, 1, 1, 1, 1])
    ...

def composite_per_genome_panel(data: AggregatedData, config: AnalysisConfig, output_path: Path):
    """Sites/alleles per genome: boxplots + GQ curves + heatmaps."""
    ...
```


---
## 11. CLI Design (`cli.py`)

```
gatk-sv-compare <subcommand> [options]

Subcommands:
  validate     Check VCF format and optionally fix issues
  preprocess   Run SVConcordance + SVRegionOverlap
  analyze      Run all analysis modules on preprocessed VCFs
  run          Run preprocess + analyze end-to-end
```

### `validate`
```
gatk-sv-compare validate \
  --vcf input.vcf.gz \
  [--fix --out fixed.vcf.gz]
```

### `preprocess`
```
gatk-sv-compare preprocess \
  --vcf-a callset_a.vcf.gz \
  --vcf-b callset_b.vcf.gz \
  --reference-dict ref.dict \
  --contig-list contigs.list \
  --output-dir /path/to/output \
  [--seg-dup-track seg_dups.bed] \
  [--simple-repeat-track simple_repeats.bed] \
  [--repeatmasker-track repeatmasker.bed] \
  [--gatk-path /path/to/gatk] \
  [--num-workers 4]
```

### `analyze`
```
gatk-sv-compare analyze \
  --vcf-a preprocessed_a.vcf.gz \
  --vcf-b preprocessed_b.vcf.gz \
  --label-a "Callset A" \
  --label-b "Callset B" \
  --output-dir /path/to/output \
  [--modules binned_counts,site_overlap_plots,...] \
  [--pass-only]              # restrict to FILTER=PASS variants only \
  [--per-chrom]              # additionally stratify all plots by contig \
  [--num-workers 4]
```

### `run`
```
gatk-sv-compare run \
  --vcf-a callset_a.vcf.gz \
  --vcf-b callset_b.vcf.gz \
  --label-a "Callset A" \
  --label-b "Callset B" \
  --reference-dict ref.dict \
  --contig-list contigs.list \
  --output-dir /path/to/output \
  [--seg-dup-track seg_dups.bed] \
  [--simple-repeat-track simple_repeats.bed] \
  [--repeatmasker-track repeatmasker.bed] \
  [--gatk-path /path/to/gatk] \
  [--modules binned_counts,site_overlap_plots,...] \
  [--pass-only]              # restrict to FILTER=PASS variants only \
  [--per-chrom]              # additionally stratify all plots by contig \
  [--num-workers 4]
```


---
## 12. Output Directory Structure

```
output_dir/
├── preprocess/
│   ├── concordance_a_eval.per_contig/    # intermediate per-contig VCFs
│   ├── concordance_b_eval.per_contig/
│   ├── concordance_a.vcf.gz              # concatenated concordance VCF (A as eval)
│   ├── concordance_a.vcf.gz.tbi
│   ├── concordance_b.vcf.gz              # concatenated concordance VCF (B as eval)
│   ├── concordance_b.vcf.gz.tbi
│   ├── annotated_a.vcf.gz               # final annotated VCF A
│   ├── annotated_a.vcf.gz.tbi
│   ├── annotated_b.vcf.gz
│   └── annotated_b.vcf.gz.tbi
│
├── main_plots/                           # composite multi-panel PNGs
│   ├── overlap_summary.{label_a}.png
│   ├── overlap_summary.{label_b}.png
│   ├── genotype_distributions.{label_a}.png
│   ├── genotype_distributions.{label_b}.png
│   ├── per_genome_summary.{label_a}.png
│   ├── per_genome_summary.{label_b}.png
│   ├── af_correlation.png
│   └── overlap_heatmaps.png
│
├── binned_counts/
│   ├── counts_a.tsv
│   ├── counts_b.tsv
│   └── ...
├── overall_counts/
│   ├── sv_count_by_type.{label_a}.png
│   ├── sv_count_by_type.{label_b}.png
│   ├── size_distribution.{label_a}.png
│   ├── ...
│   └── tables/
│       └── ...
├── genotype_dist/
│   ├── ternary.all.{label_a}.png
│   ├── ternary.all.{label_b}.png
│   ├── ternary.by_size.{label_a}.png
│   ├── carrier_freq_vs_af.{label_a}.png
│   ├── ...
│   └── tables/
│       └── hwe_stats.tsv
├── counts_per_genome/
│   ├── sites_per_genome.by_type.png
│   ├── alleles_per_genome.by_type.png
│   ├── sites_per_genome.by_size_x_class.png
│   ├── ...
│   └── tables/
│       └── ...
├── site_overlap_plots/
│   ├── overlap.by_class.{label_a}.png
│   ├── overlap.by_size.{label_a}.png
│   ├── ...
├── site_overlap_tables/
│   └── overlap_metrics.tsv
├── overlap_heatmaps/
│   ├── heatmap.size_x_freq.{label_a}.png
│   ├── heatmap.size_x_class.{label_a}.png
│   ├── heatmap.freq_x_class.{label_a}.png
│   ├── ...
├── allele_freq/
│   ├── af_correlation.overall.png
│   ├── af_correlation.by_type.png
│   ├── ...
│   └── tables/
│       └── af_correlation_stats.tsv
├── carrier_freq/
│   ├── carrier_vs_allele.{label_a}.png
│   └── carrier_vs_allele.{label_b}.png
├── genotype_quality/
│   ├── gq_histogram.overall.png
│   ├── ...
│   └── tables/
│       └── gq_summary.tsv
├── genotype_exact_match/            # only if shared samples
│   ├── exact_match.by_type.png
│   ├── ...
│   └── tables/
│       └── genotype_match_rates.tsv
└── genotype_concordance/            # only if shared samples
    ├── var_sensitivity.by_type.png
    ├── ...
    └── tables/
        └── concordance_metrics.tsv
```


---
## 13. Build Order & Implementation Phases

### Phase 1: Foundation (no GATK dependency needed for testing)

| Step | Module | Est. Lines | Depends On | Test With |
|------|--------|------------|------------|-----------|
| 1.1 | `pyproject.toml` | 40 | — | `pip install -e .` |
| 1.2 | `dimensions.py` | 120 | — | unit tests with dicts |
| 1.3 | `config.py` | 60 | 1.2 | unit tests |
| 1.4 | `vcf_format.py` | 150 | 1.2 | synthetic pysam records |
| 1.5 | `validate.py` | 200 | 1.4 | synthetic VCFs |
| 1.6 | `plot_utils.py` | 350 | 1.2 | visual spot checks |
| 1.7 | `cli.py` (skeleton) | 100 | 1.2, 1.3 | `--help` smoke test |

### Phase 2: Data Pipeline

| Step | Module | Est. Lines | Depends On | Test With |
|------|--------|------------|------------|-----------|
| 2.1 | `vcf_reader.py` | 250 | 1.2, 1.4 | synthetic VCFs |
| 2.2 | `aggregate.py` | 200 | 2.1 | synthetic VCFs |
| 2.3 | `modules/base.py` | 40 | 1.3 | — |

### Phase 3: Site-Level Modules (independent of each other)

| Step | Module | Est. Lines | Depends On | Test With |
|------|--------|------------|------------|-----------|
| 3.1 | `binned_counts.py` | 80 | 2.2, 2.3 | synthetic AggregatedData |
| 3.2 | `overall_counts.py` | 150 | 2.2, 1.6 | synthetic AggregatedData |
| 3.3 | `site_overlap_plots.py` | 150 | 2.2, 1.6 | synthetic AggregatedData |
| 3.4 | `site_overlap_tables.py` | 80 | 2.2 | synthetic AggregatedData |
| 3.5 | `overlap_heatmaps.py` | 150 | 2.2, 1.6 | synthetic AggregatedData |
| 3.6 | `allele_freq.py` | 180 | 2.2, 1.6 | synthetic AggregatedData |
| 3.7 | `carrier_freq.py` | 100 | 2.2, 1.6 | synthetic AggregatedData |
| 3.8 | `genotype_dist.py` | 200 | 2.2, 1.6 | synthetic AggregatedData |
| 3.9 | `genotype_quality.py` | 120 | 2.2, 1.6 | synthetic AggregatedData |
| 3.10 | `counts_per_genome.py` | 180 | 2.2, 1.6 | synthetic AggregatedData |

### Phase 4: Genotype-Level Modules (require shared samples)

| Step | Module | Est. Lines | Depends On | Test With |
|------|--------|------------|------------|-----------|
| 4.1 | `genotype_exact_match.py` | 200 | 2.2, 1.6 | synthetic matched VCFs |
| 4.2 | `genotype_concordance.py` | 150 | 2.2, 1.6 | synthetic concordance data |

### Phase 5: Preprocessing & Integration

| Step | Module | Est. Lines | Depends On | Test With |
|------|--------|------------|------------|-----------|
| 5.1 | `preprocess.py` | 250 | 1.3 | mock subprocess + real GATK integration test |
| 5.2 | `cli.py` (full wiring) | 200 | all | end-to-end on test VCFs |
| 5.3 | Composite figure functions | 150 | 3.x, 4.x, 1.6 | visual review |

### Phase 6: Polish

| Step | Task | Notes |
|------|------|-------|
| 6.1 | Logging | Structured logging with `logging` module, progress bars with tqdm |
| 6.2 | Error messages | Clear diagnostics for common failures (missing fields, wrong format) |
| 6.3 | Performance profiling | cProfile on a medium VCF (~10K variants, 1K samples) |
| 6.4 | README / docstrings | User-facing documentation |
| 6.5 | CI | pytest + mypy + ruff in GitHub Actions |


---
## 14. Dependencies

### `pyproject.toml` (key sections)

```toml
[project]
name = "gatk-sv-compare"
version = "0.1.0"
requires-python = ">=3.10"
dependencies = [
    "pysam>=0.21",
    "pandas>=2.0",
    "numpy>=1.24",
    "matplotlib>=3.7",
    "scipy>=1.10",          # for chi-squared HWE test
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "pytest-cov",
    "mypy>=1.0",
    "ruff>=0.1",
    "tqdm",
]

[project.scripts]
gatk-sv-compare = "gatk_sv_compare.cli:main"
```

**No dependency on:** seaborn (matplotlib only), python-ternary (hand-rolled projection),
scikit-learn, or any GATK Python bindings (subprocess calls only).


---
## 15. Key Risks & Mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| VCF format heterogeneity | Crashes or wrong results | `validate` subcommand as gate; `vcf_format.py` defensive checks; test with known-bad VCFs |
| Memory blowup on 100K-sample VCFs | OOM | Stream genotypes as numpy arrays; accumulate histogram bins, not raw values; GQ stored per-contig not globally |
| SVConcordance STATUS field missing or unexpected values | Wrong overlap counts | Validate presence in `vcf_reader.py`; fall back gracefully with warning |
| Matplotlib ternary rendering correctness | Wrong HWE plots | Unit test ternary projection against known (AA=1,AB=0,BB=0) → vertex coordinates; compare against R `HWTernaryPlot` reference outputs |
| GATK not on PATH | `preprocess` fails | Check at startup; clear error message with install instructions |
| Multiprocessing pickle failures | Crash | All pool-mapped functions are top-level (not lambdas/closures); test with `n_workers=1` and `n_workers=4` |
| BND/CTX with no SVLEN | Division by zero or wrong bucketing | BND/CTX always assigned `size_bucket = "N/A"` in `dimensions.py` |
| INS SVLEN = -1 (unknown) | Wrong size bucket | Treat SVLEN ≤ 0 as unknown → `size_bucket = "unknown"` |
| CNV sites have no standard GT | Wrong genotype counts, crash on GT access | Detect SVTYPE=CNV, skip GT-based counting, use `CN_NONREF_FREQ` for AF and `CN_NUMBER` for sample count. Exclude from HWE analysis. |
| CNV AF=`.` in annotated VCF | Parse error or NaN propagation | For CNV, use `CN_NONREF_FREQ` as the frequency metric; assign to appropriate AF bucket |
| Pre-computed count fields absent (intermediate VCFs) | Slow GT-based fallback | Detect `N_BI_GENOS` in header; if absent, warn and fall back to GT iteration with progress logging |
| Genotyped-stage VCFs with GQ 0-999 | Wrong GQ distributions | `validate` detects and `--fix` rescales; `vcf_reader` also auto-detects via header `varGQ` sentinel |
| Regenotyped VCFs with multiallelic `<CN0>,<CN1>,<CN2>,<CN3>` alleles | pysam allele parsing confusion | `validate` detects; `--fix` converts to biallelic; `vcf_reader` raises clear error if encountered |
| dDUP CPX with END==POS | Incorrect size calculation from END-POS | For CPX, use SVLEN for size (not END-POS), since dDUP has END==POS but SVLEN carries the real length |
| FILTER combinations (`HIGH_ALGORITHM_FDR;HIGH_NCR`) | Unexpected FILTER string | Parse FILTER as a set (pysam returns tuple); check membership, don't string-compare |
| Empty VCF or VCF with no variants on a contig | Index errors | Guard all aggregation with `if len(rows) == 0: return empty_dataframe()` |


---
## 16. Glossary

| Term | Meaning |
|------|---------|
| **GATK-style VCF** | Final pipeline output format: biallelic symbolic alleles, ECN present, GQ in 0–99, END semantics per spec §3.2 |
| **svtk-style VCF** | Internal intermediate format: may have breakend notation, multiallelic CNV, GQ in 0–999, no ECN |
| **SVConcordance** | GATK tool that annotates an eval VCF with match status against a truth VCF. Adds `STATUS` (MATCHED/UNMATCHED), `TRUTH_VID`, and concordance metrics to INFO. |
| **SVRegionOverlap** | GATK tool that annotates variants with genomic context (seg dup, simple repeat, repeatmasker overlap). |
| **STATUS** | INFO field added by SVConcordance. Values: `MATCHED` (eval variant has a concordant truth variant), `UNMATCHED` (no match found). |
| **AggregatedData** | The central data structure: per-variant DataFrames for both VCFs plus optional GQ/concordance data. All analysis modules read from this. |
| **VariantCategory** | A frozen dataclass of (svtype, size_bucket, af_bucket, genomic_context) — the 4D bucketing key. |
