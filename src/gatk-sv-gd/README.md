# gatk-sv-gd

**Genomic disorder (GD) CNV detection from normalized read depth, optional BAF summaries, and a curated GD table.**

Part of the [GATK-SV](https://github.com/broadinstitute/gatk-sv) pipeline.

---

## Overview

`gatk-sv-gd` is a targeted CNV caller for recurrent NAHR-mediated genomic disorder loci.
It does not discover novel loci. Instead, it scores predefined GD entries from a curated GD
table against normalized read-depth data, optional B-allele-fraction (BAF) summaries, and
known breakpoint structure. The package emits per-bin posterior tables, per-entry carrier
calls, evaluation reports, and publication-ready plots.

GD loci are hard for generic SV callers because segmental duplications distort mapping,
multiple recurrent breakpoint pairs can coexist at the same locus, and depth signal varies
substantially across bins and samples. `gatk-sv-gd` addresses those failure modes by:

- preprocessing loci with exclusion masks, ploidy-aware bin QC, flank construction, and
    optional high-resolution fallback;
- fitting a cohort-wide Bayesian model jointly across all retained GD bins;
- calling carriers either from posterior marginals or with an optional Viterbi segmentation
    layer; and
- validating calls against curated or synthetic truth tables.

---

## Command overview

The package ships seven subcommands. Five are part of the main analysis workflow and two are
supporting utilities.

```
preprocess  →  infer  →  call  →  eval  →  plot
```

| Category | Subcommand | What it does |
|----------|------------|--------------|
| Core | `preprocess` | Filter depth data, estimate ploidy, build GD bin mappings, and optionally summarize BAF |
| Core | `infer` | Fit the Bayesian model and write per-bin CN and pair-state posteriors |
| Core | `call` | Convert posteriors into per-entry GD carrier calls |
| Core | `eval` | Score calls against curated or synthetic truth tables |
| Core | `plot` | Generate carrier, truth-stratified, and summary plots |
| Utility | `extract` | Pull putative GD DEL/DUP events out of indexed VCF/BCF files |
| Utility | `synthesize` | Spike synthetic GD events into count and optional BAF tables |

`preprocess` and `infer` can still be merged by calling `infer` directly on the raw depth
matrix, but the recommended workflow is to cache preprocessing outputs and reuse them across
inference sweeps.

For end-to-end runs, the repository also includes `run_gd.sh`, a convenience wrapper that
orchestrates `preprocess`, `infer`, `call`, optional `eval`, and `plot`, and wires common
annotation inputs including pseudoautosomal-region (PAR) intervals into preprocessing.

---

## Modeling overview

### Objective and decision target

The operational target is: for each sample and each curated GD entry, estimate whether the
sample is a carrier of the expected DEL or DUP event, quantify the confidence of that claim,
and identify the best matching breakpoint pair within the locus. The outputs are designed
for thresholded calling, cohort review, downstream plots, and truth-set evaluation.

### Proposed data-generating model

For bin $b$ and sample $s$, the model uses a latent unordered diploid pair state

$$
z_{b,s} \in \{(0,0), (0,1), (0,2), (1,1), (1,2), (2,2)\},
$$

where each tuple represents per-haplotype copy number and implies a total copy number

$$
c(z_{b,s}) = h_1 + h_2.
$$

Each bin has its own pair-state prior vector

$$
\boldsymbol{\pi}_b \sim \mathrm{Dirichlet}(\boldsymbol{\alpha}), \qquad
z_{b,s} \sim \mathrm{Categorical}(\boldsymbol{\pi}_b).
$$

Normalized depth is modeled as

$$
x_{b,s} \sim \mathcal{N}\!\left(c(z_{b,s})\,\mu_b,\; (v_s + w_b)\,\frac{\kappa}{L_b}\right),
$$

where:

- $\mu_b$ is a per-bin multiplicative bias term;
- $v_s$ is a per-sample variance component;
- $w_b$ is a per-bin variance component;
- $L_b$ is the genomic width of the bin; and
- $\kappa$ is the reference bin-size factor used to scale variance for smaller bins.

When BAF summaries are available from `preprocess`, the model adds a second observation on
the minor-allele fraction:

$$
m_{b,s} \sim \mathcal{N}\!\left(m(z_{b,s}),\; \lambda_{\mathrm{BAF}}\,\sigma^2_{b,s}\right),
$$

where $m(z_{b,s})$ is the expected minor BAF for the pair state and $\sigma^2_{b,s}$ is the
upstream variance estimate derived from the number of SNP sites in the bin. CN posteriors in
the output tables are produced by collapsing pair-state posteriors onto total copy number.

### Assumptions and justifications

- `domain-supported`: the GD table captures the recurrent breakpoint pairs that matter for the
    analysis. `gatk-sv-gd` is a targeted caller and will miss novel loci or novel breakpoint
    architectures that are absent from that table.
- `domain-supported`: the depth matrix has already been normalized upstream so that diploid
    autosomal depth is approximately 2.0. The package applies further QC and masking, but it is
    not a replacement for upstream normalization.
- `domain-supported`: sample ploidy can vary by contig, especially on sex chromosomes. The
    workflow estimates ploidy during preprocessing and uses PAR intervals to avoid misleading
    chrX filtering at pseudoautosomal bins.
- `convenience-driven`: conditional on the latent pair state, depth and BAF summaries are
    modeled with Gaussian likelihoods. This is a practical approximation for normalized depth
    and summarized minor-BAF values, not a claim that the raw assay is Gaussian.
- `convenience-driven`: bins are conditionally independent given the latent state and shared
    bias/variance terms. Correlation between adjacent bins is not modeled directly unless the
    optional Viterbi caller is used downstream.
- `implementation-driven`: supported pair states only cover per-haplotype copy number 0, 1,
    or 2. Extreme amplifications and more complex allelic structures fall outside the model's
    explicit state space.
- `unverified`: a shared per-bin pair-state prior is informative enough to pool signal across
    the cohort without over-concentrating on batch-specific artifacts. That assumption should be
    stress-tested on new cohorts.

### Likelihood, priors, and inference strategy

- Continuous latent variables are inferred for `bin_bias`, `sample_var`, `bin_var`, and the
    per-bin `pair_state_probs` vectors.
- The discrete `pair_state` variable is enumerated exactly with Pyro's
    `@config_enumerate(default="parallel")` machinery.
- SVI uses `TraceEnum_ELBO` by default, or `JitTraceEnum_ELBO` unless `--disable-jit` is
    supplied.
- Two guide families are exposed: `--guide-type diagonal` (default, `AutoDiagonalNormal`) and
    `--guide-type delta` (`AutoDelta` for MAP-style sweeps).
- Current CLI defaults are intentionally conservative for fitting stability and tuning:

| Flag | Current default | Meaning |
|------|-----------------|---------|
| `--alpha-ref` | `1.0` | Dirichlet concentration for the reference pair state `(1,1)` |
| `--alpha-non-ref` | `1.0` | Dirichlet concentration for all non-reference pair states |
| `--state-prior-weight` | `0.0` | Weight of the learned pair-state prior in analytical posterior reconstruction |
| `--baf-variance-scale` | `32.0` | Downweights the BAF likelihood by inflating its variance |
| `--var-bias-bin` | `0.01` | Scale for per-bin multiplicative bias |
| `--var-sample` | `0.001` | Scale for per-sample variance |
| `--var-bin` | `0.001` | Scale for per-bin variance |
| `--max-iter` | `5000` | Maximum SVI iterations |
| `--lr-init` / `--lr-min` / `--lr-decay` | `0.02 / 0.01 / 500` | Learning-rate schedule |
| `--n-discrete-samples` | `1000` | Samples used to reconstruct discrete posteriors |

If a stronger diploid prior is needed for sparse cohorts, raise `--alpha-ref` above
`--alpha-non-ref`; the model code supports that, but the CLI no longer assumes it by default.

### Robustness strategy for sparse or poor-quality data

- `preprocess` masks bins that overlap exclusion intervals such as segmental duplications,
    centromeres, or acrocentric arms, with a bypass rule for body intervals that would
    otherwise be almost entirely masked.
- Optional `--flank-exclusion-intervals` let you mask flanks without erasing the GD body.
- Under-covered body intervals can be re-queried from `--high-res-counts` before the hard
    `--min-bins-per-interval` check is enforced.
- Rebinning and minimum-coverage checks (`--max-bins-per-interval`, `--min-rebin-coverage`,
    `--min-flank-coverage`) make the locus representation more stable when original bins are
    uneven.
- Quality filters use ploidy-aware median/MAD thresholds, and depth values are hard-clamped by
    `--clamp-threshold` before inference.
- When BAF summaries are noisy or sparse, `--baf-variance-scale` reduces their influence
    without disabling depth-based calling.

### Scaling strategy for small and large datasets

- For repeated tuning runs, split the workflow into `preprocess` and `infer`; `--preprocessed-dir`
    reuses cached bins, mappings, and optional BAF summaries.
- `--region` can restrict preprocessing or synthesis to one contig or interval when debugging.
- `--guide-type delta` is useful for quick MAP-style sweeps; `--guide-type diagonal` preserves
    posterior uncertainty.
- `--device cuda` can accelerate larger runs when a compatible GPU is available.
- `synthesize` provides a controlled way to stress-test threshold choices and scaling behavior
    before applying the model to new cohorts.

### Confidence outputs and thresholding guidance

- `infer` writes per-bin `prob_cn_*` columns and per-bin `prob_pair_*` columns, along with
    `cn_map`, `pair_state_map`, and optional BAF summaries when available.
- `call` writes `gd_cnv_calls.tsv.gz` with the preferred `confidence_score` column plus
    supporting diagnostics such as `min_interval_confidence`,
    `min_flank_non_event_confidence`, `interval_coverage`, and `is_best_match`.
- The default `posterior-marginal` calling mode selects at most one DEL and one DUP per locus,
    requiring every covered interval to exceed
    `--min-posterior-interval-confidence` (default `0.80`) and each available flank to exceed
    `--min-flank-non-event-confidence` (default `0.90`).
- `viterbi` mode is available when `--transition-matrix` is supplied; it can optionally use a
    less sticky `--breakpoint-transition-matrix` at known SD breakpoints.
- `event_marginals.tsv.gz` is always emitted for downstream plotting. In Viterbi mode,
    `viterbi_paths.tsv.gz` provides the segmented path.

### Validation and falsification plan

- `eval` accepts either a curated BED-style truth table or the two-column truth table emitted
    by `synthesize`.
- Pass `preprocess/gd_table_filtered.tsv` to `call`, `plot`, and `eval` so that evaluation only
    scores the loci that actually survived preprocessing.
- `plot --eval-report` renders separate true-positive, false-positive, and false-negative PDF
    bundles instead of the default carrier PDF.
- `extract` can provide orthogonal evidence from VCF/BCF inputs when a cohort already has SV
    calls from other methods.
- Before trusting thresholds on a new cohort, check calibration with `synthesize` and review the
    resulting `confidence_distribution.png`, carrier plots, and truth report together.

### Expected failure modes

- Novel or atypical events that do not match any GD entry will not be called cleanly.
- Extremely noisy normalization can make the Gaussian depth model and learned bias terms absorb
    real signal or create false positives.
- Loci with too few surviving bins or too little flank coverage can fail preprocessing or yield
    weak posterior-marginal confidence.
- If BAF sample identifiers do not match the depth matrix exactly, BAF support will be dropped.
- Copy states outside the modeled pair-state space, especially high-level amplifications, are
    not represented explicitly.

---

## Installation

```bash
pip install -e src/gatk-sv-gd
```

With development dependencies:

```bash
pip install -e "src/gatk-sv-gd[dev]"
```

Runtime dependencies are declared in `pyproject.toml` and currently include Python 3.9+,
NumPy, pandas, pysam, PyTorch, Pyro, matplotlib, tqdm, and `intervaltree`.

For workflows that use optional high-resolution depth or BAF tables, follow the CLI help on
which files must be bgzipped and tabix-indexed.

---

## Usage

### Wrapper script

`run_gd.sh` is the simplest way to run the supported preprocess → infer → call → eval → plot
workflow with the current defaults:

```bash
src/gatk-sv-gd/run_gd.sh \
        --work-dir work/gd_run \
        --input-depth counts.rd.txt.gz \
        --high-res-depth counts.high_res.rd.txt.gz \
        --baf-table all_samples.baf.txt.gz \
        --gd-table gd_table.tsv \
        --segdup-bed segdups.bed.gz \
        --centromere-bed centromeres.bed.gz \
        --acrocentric-arm-bed acrocentric_arms.bed.gz \
        --gaps-bed gaps.bed.gz \
        --gtf genes.gtf.gz \
        --par-bed hg38_par.bed \
        --truth-table truth.tsv
```

The wrapper passes `preprocess/gd_table_filtered.tsv` and `preprocess/ploidy_estimates.tsv`
to downstream steps automatically and does not require transition matrices unless you opt into
`call --calling-mode viterbi` through `--call-args`.

### Raw subcommands

```bash
# Step 1: preprocess
gatk-sv-gd preprocess \
        --input counts.rd.txt.gz \
        --gd-table gd_table.tsv \
        --exclusion-intervals segdups.bed.gz centromeres.bed.gz acrocentric_arms.bed.gz \
        --par-intervals hg38_par.bed \
        --high-res-counts counts.high_res.rd.txt.gz \
        --baf-table all_samples.baf.txt.gz \
        --output-dir preprocess/

# Step 2: infer from cached preprocess output
gatk-sv-gd infer \
        --preprocessed-dir preprocess/ \
        --output-dir infer/

# Step 3: call from posterior marginals (default)
gatk-sv-gd call \
        --cn-posteriors infer/cn_posteriors.tsv.gz \
        --bin-mappings preprocess/bin_mappings.tsv.gz \
        --gd-table preprocess/gd_table_filtered.tsv \
        --ploidy-table preprocess/ploidy_estimates.tsv \
        --output-dir call/

# Step 4: evaluate against curated or synthetic truth
gatk-sv-gd eval \
        --calls call/gd_cnv_calls.tsv.gz \
        --truth-table truth.tsv \
        --gd-table preprocess/gd_table_filtered.tsv \
        --ploidy-table preprocess/ploidy_estimates.tsv \
        --output-dir eval/

# Step 5: plot calls (optionally stratified by eval report)
gatk-sv-gd plot \
        --calls call/gd_cnv_calls.tsv.gz \
        --cn-posteriors infer/cn_posteriors.tsv.gz \
        --gd-table preprocess/gd_table_filtered.tsv \
        --ploidy-table preprocess/ploidy_estimates.tsv \
        --event-marginals call/event_marginals.tsv.gz \
        --eval-report eval/truth_evaluation_report.tsv \
        --output-dir plot/
```

`infer` also supports a one-step mode that reruns preprocessing internally:

```bash
gatk-sv-gd infer \
        --input counts.rd.txt.gz \
        --gd-table gd_table.tsv \
        --exclusion-intervals segdups.bed.gz \
        --output-dir infer/
```

Run `gatk-sv-gd <subcommand> --help` for the complete option surface.

---

## Subcommand reference

### `preprocess`

Loads depth data, applies ploidy-aware QC and exclusion masking, constructs locus flanks,
optionally replaces under-covered regions from a high-resolution count table, and writes the
cached inputs consumed by `infer`.

```bash
gatk-sv-gd preprocess -i counts.rd.txt.gz -g gd_table.tsv -o preprocess/ [options]
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--exclusion-intervals` | `[]` | BED files masked in both body and flanks |
| `--flank-exclusion-intervals` | `[]` | BED files masked only in flanks |
| `--par-intervals` | `[]` | PAR BED files ignored during ploidy-aware filtering |
| `--high-res-counts` | — | bgzipped, tabix-indexed high-resolution depth table for fallback |
| `--baf-table` | — | bgzipped, tabix-indexed BAF table to summarize over retained bins |
| `--region` | — | Restrict preprocessing to specific contigs or intervals |
| `--exclusion-threshold` | `0.5` | Minimum overlap fraction to mask a bin |
| `--exclusion-bypass-threshold` | `0.6` | Skip masking when a body interval is mostly excluded |
| `--min-bins-per-interval` | `10` | Hard minimum bins per body interval |
| `--max-bins-per-interval` | `20` | Rebin body intervals down to at most this many bins |
| `--min-rebin-coverage` | `0.5` | Minimum coverage fraction for a rebinned bin |
| `--min-flank-bases` / `--max-flank-bases` | `50000 / 1000000` | Bounds on flank growth |
| `--min-flank-bins` / `--min-flank-coverage` | `10 / 0.5` | Minimum flank support |
| `--median-min` / `--median-max` / `--mad-max` | `1.0 / 3.0 / 1.0` | Per-bin QC thresholds |

**Outputs:**

| File | Description |
|------|-------------|
| `preprocessed_bins.tsv.gz` | Combined retained bin matrix |
| `bin_mappings.tsv.gz` | Per-bin locus, interval, and genomic coordinates |
| `locus_intervals.tsv.gz` | Retained locus-level interval metadata |
| `gd_entry_intervals.tsv.gz` | Retained per-entry breakpoint-interval metadata |
| `ploidy_estimates.tsv` | Per-sample, per-contig ploidy estimates |
| `gd_table_filtered.tsv` | GD table restricted to the loci that survived preprocessing |
| `preprocessed_baf.tsv.gz` | ROI-filtered BAF rows retained from the input BAF table |
| `preprocessed_baf_summary.tsv.gz` | Per-bin, per-sample BAF summaries used by `infer` |

Pass `gd_table_filtered.tsv` and `ploidy_estimates.tsv` to downstream `call`, `plot`, and
`eval` steps.

---

### `infer`

Fits the Bayesian model and writes posterior tables. With `--preprocessed-dir`, it reuses the
cached outputs from `preprocess`, including optional BAF summaries.

```bash
gatk-sv-gd infer \
        ( --input counts.rd.txt.gz --gd-table gd_table.tsv | --preprocessed-dir preprocess/ ) \
        --output-dir infer/ [options]
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--preprocessed-dir` | — | Reuse cached preprocessing outputs |
| `--alpha-ref` / `--alpha-non-ref` | `1.0 / 1.0` | Pair-state prior concentrations |
| `--state-prior-weight` | `0.0` | Weight of the learned pair-state prior in analytical posterior reconstruction |
| `--baf-variance-scale` | `32.0` | Downweight BAF evidence by inflating variance |
| `--var-bias-bin` / `--var-sample` / `--var-bin` | `0.01 / 0.001 / 0.001` | Continuous prior scales |
| `--bin-size-factor` | `10000.0` | Variance scaling reference bin size |
| `--guide-type` | `diagonal` | `diagonal` or `delta` |
| `--max-iter` | `5000` | Maximum SVI iterations |
| `--lr-init` / `--lr-min` / `--lr-decay` | `0.02 / 0.01 / 500` | Learning-rate schedule |
| `--early-stopping` | on | Disable with `--no-early-stopping` |
| `--patience` / `--min-delta` | `50 / 1000.0` | Early-stopping controls |
| `--n-discrete-samples` | `1000` | Samples used to reconstruct discrete posteriors |
| `--device` | `cpu` | `cpu` or `cuda` |

**Outputs:**

| File | When written | Description |
|------|-------------|-------------|
| `cn_posteriors.tsv.gz` | always | Per-bin, per-sample CN and pair-state posterior probabilities |
| `sample_posteriors.tsv.gz` | always | Per-sample variance estimates |
| `bin_posteriors.tsv.gz` | always | Per-bin bias, variance, and learned prior estimates |
| `bin_mappings.tsv.gz` | raw-input mode | Bin mapping metadata when `infer` performs preprocessing internally |
| `locus_intervals.tsv.gz` | raw-input mode | Retained locus intervals when `infer` performs preprocessing internally |
| `gd_entry_intervals.tsv.gz` | raw-input mode | Retained GD-entry intervals when `infer` performs preprocessing internally |
| `ploidy_estimates.tsv` | raw-input mode | Ploidy estimates when `infer` performs preprocessing internally |

---

### `call`

Converts posterior tables into per-entry GD carrier calls.

```bash
gatk-sv-gd call \
        --cn-posteriors infer/cn_posteriors.tsv.gz \
        --bin-mappings preprocess/bin_mappings.tsv.gz \
        --gd-table preprocess/gd_table_filtered.tsv \
        --ploidy-table preprocess/ploidy_estimates.tsv \
        --output-dir call/ [options]
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--calling-mode` | `posterior-marginal` | `posterior-marginal` or `viterbi` |
| `--min-posterior-interval-confidence` | `0.80` | Minimum interval confidence in posterior-marginal mode |
| `--min-flank-non-event-confidence` | `0.90` | Minimum flank non-event confidence in posterior-marginal mode |
| `--transition-matrix` | — | Required for `--calling-mode viterbi` |
| `--breakpoint-transition-matrix` | — | Optional breakpoint-specific transition matrix for Viterbi |
| `--min-mean-coverage` | `0.50` | Minimum per-interval support for Viterbi selection |
| `--ploidy-table` | — | Sample/contig ploidy table; diploid is assumed if omitted |

**Outputs:**

| File | Description |
|------|-------------|
| `gd_cnv_calls.tsv.gz` | One row per sample × GD entry with confidence and match metadata |
| `event_marginals.tsv.gz` | Per-bin DEL/DUP event probabilities used by `plot` |
| `viterbi_paths.tsv.gz` | Segmented path output for Viterbi runs |

Prefer `confidence_score` for thresholding; `log_prob_score` is retained as a compatibility
column.

---

### `eval`

Scores calls against either a curated truth table or the synthetic truth table generated by
`synthesize`.

```bash
gatk-sv-gd eval \
        --calls call/gd_cnv_calls.tsv.gz \
        --truth-table truth.tsv \
        --ploidy-table preprocess/ploidy_estimates.tsv \
        --output-dir eval/ [options]
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--gd-table` | — | Filtered GD table from `preprocess`; recommended for fair evaluation |
| `--ploidy-table` | *(required)* | Defines the batch sample set and contig ploidy |
| `--min-confidence` | — | Optional confidence threshold for counting predictions |

**Output:** `truth_evaluation_report.tsv`

---

### `plot`

Generates per-sample, per-locus plots and summary figures.

```bash
gatk-sv-gd plot \
        --calls call/gd_cnv_calls.tsv.gz \
        --cn-posteriors infer/cn_posteriors.tsv.gz \
        --gd-table preprocess/gd_table_filtered.tsv \
        --output-dir plot/ [options]
```

**Notable options:**

| Flag | Purpose |
|------|---------|
| `--ploidy-table` | Use per-contig ploidy in annotations and scaling |
| `--raw-counts` / `--high-res-counts` | Overlay raw depth profiles |
| `--gtf`, `--segdup-bed`, `--gaps-bed` | Overlay genes, segmental duplications, and assembly gaps |
| `--event-marginals` | Overlay per-bin DEL/DUP event probabilities |
| `--viterbi-paths` | Overlay Viterbi segmentation |
| `--eval-report` | Write `true_positives.pdf`, `false_positives.pdf`, and `false_negatives.pdf` |
| `--plot-all-samples`, `--sample`, `--locus` | Restrict plotting scope |
| `--carrier-confidence-threshold` | Filter which calls appear in `carrier_plots.pdf` |

Common outputs include `carrier_plots.pdf`, `carrier_summary.png`,
`confidence_distribution.png`, and, when `--eval-report` is supplied, truth-stratified PDFs.

---

### `extract`

Queries indexed VCF/BCF files over GD loci and annotates overlapping DEL/DUP records.

```bash
gatk-sv-gd extract \
        --vcf cohort.vcf.gz \
        --gd-table gd_table.tsv \
        --output-dir extract/
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--vcf` / `--vcf-list` | *(one required)* | One or more indexed VCF/BCF inputs |
| `--reciprocal-overlap` | `0.5` | Canonical NAHR overlap threshold |
| `--atypical-coverage` | `0.7` | Minimum GD coverage before labeling as canonical rather than atypical |
| `--non-nahr-overlap` | `0.01` | Minimum overlap for non-NAHR annotation |

**Outputs:** `gd_variants.vcf.gz`, `gd_variants.bed`

---

### `synthesize`

Spikes synthetic GD events into low-resolution counts, high-resolution counts, and optional
BAF tables, and writes the truth tables used by `eval`.

```bash
gatk-sv-gd synthesize \
        --lo-res-counts counts.rd.txt.gz \
        --hi-res-counts counts.high_res.rd.txt.gz \
        --baf-table all_samples.baf.txt.gz \
        --ploidy-table preprocess/ploidy_estimates.tsv \
        --gd-table gd_table.tsv \
        --output-dir synth/
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--gd-probability` | `0.5` | Chance that a sample gets a GD event |
| `--salted-event-probability` | `0.20` | Chance of adding nuisance background DEL/DUP events |
| `--viable-trisomy-probability` | `0.10` | Chance of adding a whole-chromosome viable trisomy/YY event |
| `--truth-table` | — | Reuse an existing synthetic truth assignment |
| `--region` | — | Restrict eligible GD entries |
| `--del-multiplier` / `--dup-multiplier` | `0.5 / 1.5` | Multipliers applied to spiked counts |

**Outputs:**

| File | Description |
|------|-------------|
| `lo_res_counts.synthesized.rd.txt.gz` | Low-resolution depth table with spike-ins |
| `hi_res_counts.synthesized.rd.txt.gz` | High-resolution depth table with spike-ins |
| `all_samples.synthesized.baf.txt.gz` | BAF table with spike-ins |
| `truth_table.tsv` | Sample-to-GD truth table accepted by `eval` |
| `background_events.tsv` | Additional nuisance/background events |

---

## Input formats

### Depth TSV (`-i / --input`)

Tab-separated. First three columns must be `Chr`, `Start`, `End`; remaining columns are
sample IDs. Values should be normalised so that diploid autosomal coverage ≈ 2.0.

```
Chr     Start   End     SAMPLE1  SAMPLE2  ...
chr1    0       500     1.98     2.01     ...
chr1    500     1000    2.03     1.97     ...
```

### GD table (`-g / --gd-table`)

Tab-separated. Required columns:

| Column | Description |
|--------|-------------|
| `chr` | Chromosome |
| `start_GRCh38`, `end_GRCh38` | GD interval coordinates (hg38) |
| `GD_ID` | Unique identifier for this DEL/DUP entry |
| `svtype` | `DEL` or `DUP` |
| `NAHR` | `yes` / `no` |
| `terminal` | `yes` / `no` (terminal deletions extending to the chromosome arm boundary) |
| `cluster` | Locus name grouping all entries at the same locus |
| `BP1`, `BP2` | Named breakpoints bounding this entry |

Multiple rows per `cluster` encode different event sizes (e.g. BP1–BP2, BP1–BP3, BP2–BP3
at 22q11.2). Breakpoints are sorted by genomic coordinate, not by name, so that complex
naming schemes (e.g. numeric IDs mixed with gene names) are handled correctly.

When you run `preprocess`, use `gd_table_filtered.tsv` downstream instead of the original GD
table so that calling, plotting, and evaluation stay aligned with the loci that survived bin
filtering.

### BAF table (`--baf-table`)

The optional BAF input must be bgzipped and tabix-indexed because `preprocess` queries it by
genomic interval rather than loading it whole. The expected logical columns are:

| Column | Description |
|--------|-------------|
| `Chr` | Chromosome |
| `Pos` | 1-based SNP position |
| `BAF` | B-allele fraction in `[0, 1]` |
| `Sample` | Sample identifier matching the depth matrix |

Headerless SiteDepthtoBAF-style input is accepted as long as the first four columns follow
that layout.

### Truth table (`--truth-table`)

`eval` accepts either:

- a curated BED-like table with columns `#chrom`, `start`, `end`, `name`, `svtype`,
  `samples`, `NAHR_GD`, and `NAHR_GD_atypical`; or
- the two-column synthetic truth table from `synthesize` with columns `sample_id` and `GD_ID`.

---

## Module reference

| Module | Description |
|--------|-------------|
| `cli.py` | CLI entry point; dispatches to subcommand `parse_args` / `main` functions |
| `models.py` | `GDLocus`, `GDTable` — locus and breakpoint data structures |
| `depth.py` | `ExclusionMask`, `DepthData`, `CNVModel` — Pyro model and SVI training |
| `bins.py` | Bin extraction, quality filtering, rebinning, interval assignment, I/O |
| `highres.py` | Tabix-indexed high-resolution depth queries and normalisation |
| `preprocess.py` | Per-locus bin collection pipeline; `preprocess` subcommand |
| `infer.py` | Model training and CN posterior computation; `infer` subcommand |
| `viterbi.py` | Viterbi HMM segmentation algorithm |
| `call.py` | CNV calling (posterior-marginal and Viterbi strategies); `call` subcommand |
| `output.py` | Posterior table writers and ploidy estimation |
| `plot.py` | Depth and CN-state visualisation; `plot` subcommand |
| `annotations.py` | Gene and genomic feature annotation overlays for plots |
| `eval.py` | Sensitivity/precision evaluation; `eval` subcommand |
| `extract.py` | VCF/BCF extraction utility; `extract` subcommand |
| `synthesize.py` | Synthetic spike-in generator; `synthesize` subcommand |

---

## Development

```bash
pip install -e "src/gatk-sv-gd[dev]"
pytest src/gatk-sv-gd/tests
```

---

## References

1. **Pyro** — Bingham, E. *et al.* (2019). Pyro: Deep Universal Probabilistic Programming.
   *Journal of Machine Learning Research*, 20(28), 1–6.
   <https://jmlr.org/papers/v20/18-403.html>

2. **Stochastic Variational Inference** — Hoffman, M. D., Blei, D. M., Wang, C., &
   Paisley, J. (2013). Stochastic Variational Inference.
   *Journal of Machine Learning Research*, 14, 1303–1347.

3. **Variational Bayes / ELBO** — Jordan, M. I., Ghahramani, Z., Jaakkola, T. S., & Saul,
   L. K. (1999). An Introduction to Variational Methods for Graphical Models.
   *Machine Learning*, 37, 183–233.

4. **Viterbi algorithm** — Viterbi, A. J. (1967). Error Bounds for Convolutional Codes and
   an Asymptotically Optimum Decoding Algorithm.
   *IEEE Transactions on Information Theory*, 13(2), 260–269.

5. **Genomic disorders and NAHR** — Lupski, J. R. (1998). Genomic disorders: structural
   features of the genome can lead to DNA rearrangements and human disease traits.
   *Trends in Genetics*, 14(10), 417–422.

6. **Segmental duplications and recurrent CNVs** — Stankiewicz, P., & Lupski, J. R. (2002).
   Genome architecture, rearrangements and genomic disorders.
   *Trends in Genetics*, 18(2), 74–82.

7. **GATK-SV** — Collins, R. L. *et al.* (2020). A structural variation reference for
   medical and population genetics. *Nature*, 581, 444–451.
   <https://doi.org/10.1038/s41586-020-2287-8>
