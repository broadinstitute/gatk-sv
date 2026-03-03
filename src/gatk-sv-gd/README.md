# gatk-sv-gd

**Genomic Disorder (GD) CNV detection from binned read counts.**

Part of the [GATK-SV](https://github.com/broadinstitute/gatk-sv) pipeline.

---

## Overview

`gatk-sv-gd` detects copy-number variants (CNVs) at genomic disorder (GD) loci using a
hierarchical Bayesian model trained on binned read-depth data. GD loci are recurrent CNV
hotspots whose deletions and duplications are driven by non-allelic homologous recombination
(NAHR) between flanking segmental duplications (SDs). Classic examples include 22q11.2
(DiGeorge / velocardiofacial syndrome), 7q11.23 (Williams–Beuren syndrome), 17p11.2
(Smith–Magenis syndrome), and 15q11–13 (Prader–Willi / Angelman syndrome).

Standard SV callers are challenged at these loci because: (i) SD boundaries create
read-mapping ambiguity; (ii) multiple distinct recurrent breakpoint pairs produce events of
varying sizes at the same locus; and (iii) large batch effects in read-depth data require
careful normalisation. This package addresses these challenges through locus-aware bin
preprocessing, exclusion masking, and a unified Bayesian model trained jointly across all
loci in a cohort.

---

## Pipeline overview

The five subcommands form a linear pipeline:

```
preprocess  →  infer  →  call  →  plot  →  eval
```

| Step | Subcommand | What it does |
|------|------------|--------------|
| 1 | `preprocess` | Load and filter depth data; collect per-locus bin matrices |
| 2 | `infer` | Train hierarchical Bayesian model; compute per-bin CN posteriors |
| 3 | `call` | Call carrier samples from posteriors (MLP or Viterbi strategy) |
| 4 | `plot` | Generate per-locus depth and CN-state visualisations |
| 5 | `eval` | Benchmark calls against a manually curated truth table |

`preprocess` and `infer` can be run as separate steps (recommended for large cohorts, where
preprocessing is slow but inference may be re-run with different hyperparameters) or merged
into a single `infer` invocation by passing the raw depth file directly.

---

## Background

### Genomic disorders and NAHR

Genomic disorders arise from NAHR between highly similar (≥ 97% identity) low-copy repeats
(LCRs). Because NAHR can occur between different pairs of LCR copies at the same locus, a
single disease region often harbours multiple recurrent event sizes corresponding to
different breakpoint-pair combinations — e.g. BP1–BP2, BP1–BP3, and BP2–BP3 events at
22q11.2. The GD table (see [Input formats](#input-formats)) encodes these combinations
explicitly as named breakpoints (BP1, BP2, BP3, …), and this package models each combination
as a separate GD entry within the same cluster.

### Read-depth normalisation

The input depth file is expected to contain *pre-normalised* coverage in units where a
diploid (CN = 2) autosomal sample has a median value of approximately 2.0. Per-sample
autosomal-median scaling and GC-bias correction are assumed to have been applied upstream
(e.g., by GATK-SV's `CollectCoverage` → `NormalizeCoverage` steps). Within this package,
additional quality filtering is applied (per-bin median depth and MAD thresholds), and
extreme values are hard-clamped before model training.

### Hierarchical Bayesian model

The core model (`CNVModel` in `depth.py`) is a hierarchical Bayesian model implemented in
[Pyro](https://pyro.ai/), a probabilistic programming framework built on PyTorch. All GD
loci across the entire cohort are pooled into a single joint model, enabling shared
hyperparameter estimation.

**Latent variables:**

| Variable | Prior | Interpretation |
|----------|-------|----------------|
| $\boldsymbol{\pi}_b \in \Delta^{K-1}$ | $\text{Dirichlet}(\boldsymbol{\alpha})$ | Per-bin CN state probability vector |
| $\text{cn}_{b,s} \in \{0,\ldots,K-1\}$ | $\text{Categorical}(\boldsymbol{\pi}_b)$ | Copy number of sample $s$ at bin $b$ |
| $\mu^{\text{bias}}_b > 0$ | $\text{LogNormal}(0,\, \sigma_{\text{bias}})$ | Per-bin multiplicative depth bias (GC, mappability) |
| $v_s > 0$ | $\text{Exponential}(1/\sigma_{\text{sample}})$ | Per-sample variance component |
| $w_b > 0$ | $\text{Exponential}(1/\sigma_{\text{bin}})$ | Per-bin variance component |

**Observation model:**

$$x_{b,s} \sim \mathcal{N}\!\bigl(\,\text{cn}_{b,s}\cdot\mu^{\text{bias}}_b,\;\; v_s + w_b\,\bigr)$$

where $x_{b,s}$ is the normalised depth at bin $b$ for sample $s$, and $K = 6$
(states CN 0–5).

The Dirichlet concentration vector $\boldsymbol{\alpha}$ is parameterised so that CN = 2
carries a large weight ($\alpha_{\text{ref}}$, default 50) while all other states share a
small weight ($\alpha_{\text{non-ref}}$, default 1), creating a strong prior toward diploid.
The per-bin distribution $\boldsymbol{\pi}_b$ is shared across all samples at that bin,
so cohort-level evidence is pooled.

**Inference:**

The discrete latent variable $\text{cn}_{b,s}$ is enumerated out exactly in parallel using
Pyro's `@config_enumerate` decorator and `TraceEnum_ELBO`. The remaining continuous
latents $(\mu^{\text{bias}}_b, v_s, w_b, \boldsymbol{\pi}_b)$ are approximated via
stochastic variational inference (SVI), maximising the evidence lower bound (ELBO):

$$\mathcal{L}(\phi) = \mathbb{E}_{q_\phi}\!\left[\log p(\mathbf{x}, \mathbf{z}) - \log q_\phi(\mathbf{z})\right]$$

Two variational guide families are supported:

- **`diagonal`** (default): `AutoDiagonalNormal` — mean-field Gaussian with learned mean
  and diagonal covariance.
- **`delta`**: `AutoDelta` — point-mass (MAP) guide; faster but no posterior uncertainty.

Optimisation uses Adam with an exponential learning-rate schedule
($\eta(k) = \eta_{\min} + (\eta_0 - \eta_{\min})\,e^{-k/\tau}$) and optional early
stopping based on ELBO stagnation.

### CNV calling

After inference two calling strategies are available:

#### 1. Mean log-probability (MLP, default)

For each GD entry the posterior probabilities $p(\text{cn}_{b,s} \neq \text{ploidy})$
are aggregated across the bins in the covered body interval(s) into a bin-count-weighted
mean log-probability score:

$$\ell = \frac{\sum_{i}\, n_i \log p_i(\text{CN} \neq \text{ploidy})}{\sum_{i} n_i}$$

where the sum runs over body intervals $i$ with $n_i$ bins. A call is made when
$\ell > \ell_{\text{threshold}}$ (default −0.3). Flanking regions (100 % of the locus
size on each side) are scored separately: a call whose flanks also exceed a second
threshold is reclassified as a large *spanning* CNV rather than a true GD carrier.

#### 2. Viterbi HMM segmentation

When a CN-state transition probability matrix is provided (`--transition-matrix`), the
Viterbi algorithm is run on the per-bin log-posteriors from the Pyro model, treating them
as per-position state priors. The path $\mathbf{s}^* = (s_1^*,\ldots,s_T^*)$ maximising

$$\log P(\mathbf{s}) = \sum_{t=1}^{T} \log p(s_t \mid \mathbf{x}_t)\;+\;\sum_{t=2}^{T} \log A(s_{t-1},\, s_t)$$

is found by standard Viterbi recursion, where $A$ is the user-supplied transition matrix and
$p(s_t \mid \mathbf{x}_t)$ is the Pyro-inferred CN posterior at bin $t$. The resulting
copy-number segmentation is compared against the expected breakpoint pattern for each GD
entry (body intervals at non-ploidy CN, flanks at ploidy CN). An optional
`--breakpoint-transition-matrix` applies a less-sticky transition matrix specifically at the
known SD-block boundaries to actively encourage state changes at breakpoints.

---

## Installation

```bash
pip install -e src/gatk-sv-gd
```

With conda:

```bash
conda activate <env>
pip install -e src/gatk-sv-gd
```

**Requirements:** Python ≥ 3.9, PyTorch, Pyro-PPL, NumPy, pandas, pysam, matplotlib, tqdm.

---

## Usage

### Full two-step pipeline (recommended for large cohorts)

```bash
# Step 1 — Preprocess: load, filter, and cache per-locus bin matrices
gatk-sv-gd preprocess \
    --input counts.tsv \
    --gd-table gd_table.tsv \
    --exclusion-intervals segdups.bed.gz centromeres.bed.gz \
    --output-dir preprocessed/

# Step 2 — Infer: train Bayesian model and compute CN posteriors
gatk-sv-gd infer \
    --preprocessed-dir preprocessed/ \
    --output-dir results/

# Step 3 — Call: call carrier samples from posteriors
gatk-sv-gd call \
    --cn-posteriors results/cn_posteriors.tsv.gz \
    --bin-mappings preprocessed/bin_mappings.tsv.gz \
    --gd-table gd_table.tsv \
    --output-dir results/

# Step 4 — Plot: visualise depth and CN-state posteriors
gatk-sv-gd plot \
    --calls results/gd_cnv_calls.tsv.gz \
    --cn-posteriors results/cn_posteriors.tsv.gz \
    --gd-table gd_table.tsv \
    --output-dir plots/

# Step 5 — Eval: benchmark against truth set (optional)
gatk-sv-gd eval \
    --calls results/gd_cnv_calls.tsv.gz \
    --truth-table truth.tsv \
    --output-dir results/
```

### Single-step inference (skips separate preprocessing)

```bash
gatk-sv-gd infer \
    --input counts.tsv \
    --gd-table gd_table.tsv \
    --exclusion-intervals segdups.bed.gz \
    --output-dir results/
```

Run `gatk-sv-gd <subcommand> --help` for full option details.

---

## Subcommand reference

### `preprocess`

Loads depth data, applies quality and exclusion filtering, estimates ploidy, and collects
per-locus bin matrices. Results are cached to disk for consumption by `infer`.

```
gatk-sv-gd preprocess -i counts.tsv -g gd_table.tsv -o preprocessed/ [options]
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `-i / --input` | *(required)* | Normalised depth TSV (bins × samples) |
| `-g / --gd-table` | *(required)* | GD locus definition table |
| `-e / --exclusion-intervals` | `[]` | BED files of regions to mask (segdups, centromeres, etc.) |
| `--high-res-counts` | — | bgzipped + tabix-indexed high-res count file; used as fallback for under-covered intervals |
| `--locus-padding` | 10 000 | Padding around locus outer boundaries (bp) |
| `--exclusion-threshold` | 1.0 | Bin overlap fraction required to apply exclusion mask |
| `--exclusion-bypass-threshold` | 0.6 | Skip exclusion masking for intervals where ≥ this fraction is masked (preserves coverage in heavily-duplicated regions) |
| `--min-bins-per-region` | 10 | Minimum bins per body interval; under-covered intervals trigger high-res fallback |
| `--max-bins-per-interval` | 20 | Maximum bins per interval after rebinning |
| `--min-rebin-coverage` | 0.5 | Minimum original-bin coverage fraction required for a rebinned bin |
| `--min-flank-bases` | 50 000 | Minimum cumulative bp each flank must cover |
| `--min-flank-bins` | 10 | Minimum number of bins per flank |
| `--median-min / --median-max` | 1.0 / 3.0 | Per-bin median depth quality bounds |
| `--mad-max` | 2.0 | Per-bin MAD upper bound |

**Outputs** (written to `--output-dir`):

| File | Description |
|------|-------------|
| `preprocessed_bins.tsv.gz` | Combined bin depth matrix (all loci concatenated) |
| `bin_mappings.tsv.gz` | Per-bin locus, interval, and coordinate metadata |
| `locus_metadata.tsv` | Per-locus inclusion/exclusion summary |
| `ploidy_table.tsv` | Per-sample ploidy estimates (used for sex-chromosome correction) |

---

### `infer`

Trains the hierarchical Bayesian model on the combined bin matrix and computes CN posteriors
for every bin × sample combination.

```
gatk-sv-gd infer  ( -i counts.tsv -g gd_table.tsv  |  --preprocessed-dir preprocessed/ )
                  -o results/ [options]
```

All preprocessing flags from `preprocess` are also accepted when running without
`--preprocessed-dir`.

**Model / inference options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--preprocessed-dir` | — | Load cached bins from a `preprocess` output directory; skips preprocessing |
| `--alpha-ref` | 50.0 | Dirichlet concentration for CN = 2 (strength of diploid prior) |
| `--alpha-non-ref` | 1.0 | Dirichlet concentration for non-diploid states |
| `--var-bias-bin` | 0.1 | Per-bin bias log-normal variance ($\sigma_{\text{bias}}$) |
| `--var-sample` | 0.2 | Per-sample variance scale ($\sigma_{\text{sample}}$) |
| `--var-bin` | 0.2 | Per-bin variance scale ($\sigma_{\text{bin}}$) |
| `--guide-type` | `diagonal` | Variational guide: `diagonal` (mean-field) or `delta` (MAP) |
| `--max-iter` | 1 000 | Maximum SVI iterations |
| `--lr-init / --lr-min / --lr-decay` | 0.01 / 0.001 / 500 | Learning rate schedule parameters |
| `--early-stopping` | on | Enabled by default; disable with `--no-early-stopping` |
| `--patience` | 50 | Early-stopping patience (epochs without improvement) |
| `--min-delta` | 1 000.0 | Minimum ELBO improvement to reset patience counter |
| `--clamp-threshold` | 5.0 | Clamp normalised depth values above this ceiling |
| `--n-discrete-samples` | — | Samples for discrete posterior integration |
| `--device` | `cpu` | `cpu` or `cuda` |
| `--disable-jit` | off | Disable Pyro JIT compilation (useful for debugging) |

**Outputs** (written to `--output-dir`):

| File | Description |
|------|-------------|
| `cn_posteriors.tsv.gz` | Per-bin × per-sample CN posterior probabilities (`prob_cn_0` … `prob_cn_5`) and MAP copy-number estimate |
| `sample_posteriors.tsv.gz` | Per-sample variance MAP estimates |
| `bin_posteriors.tsv.gz` | Per-bin bias and variance MAP estimates |

---

### `call`

Calls GD carrier samples from CN posteriors using either the mean log-probability (MLP) or
Viterbi strategy.

```
gatk-sv-gd call \
    --cn-posteriors results/cn_posteriors.tsv.gz \
    --bin-mappings preprocessed/bin_mappings.tsv.gz \
    --gd-table gd_table.tsv \
    --output-dir results/ [options]
```

**Key options:**

| Flag | Default | Description |
|------|---------|-------------|
| `--log-prob-threshold` | −0.3 | MLP score threshold for calling a carrier |
| `--flanking-log-prob-threshold` | −1.0 | MLP score in flanking regions to reclassify as spanning CNV |
| `--transition-matrix` | — | TSV of CN-state transition probabilities; enables Viterbi calling |
| `--breakpoint-transition-matrix` | — | More permissive transition matrix applied at SD-block boundaries |
| `--ploidy-table` | — | Per-sample ploidy overrides (sex chromosomes) |

**Output:** `gd_cnv_calls.tsv.gz` — one row per sample × GD entry with call status
(`is_carrier`, `is_spanning`), MLP/Viterbi score, and breakpoint assignment.

---

### `plot`

Generates per-locus visualisations showing normalised depth, CN-state posteriors, and call
annotations for each sample.

```
gatk-sv-gd plot \
    --calls results/gd_cnv_calls.tsv.gz \
    --cn-posteriors results/cn_posteriors.tsv.gz \
    --gd-table gd_table.tsv \
    --output-dir plots/ [options]
```

Optional overlays: `--gtf` (gene annotations), `--segdup-bed` (SD tracks), `--gaps-bed`
(assembly gaps), `--raw-counts` (raw depth beneath normalised signal).

---

### `eval`

Benchmarks calls against a manually curated truth table and reports per-site and aggregate
sensitivity and precision.

```
gatk-sv-gd eval \
    --calls results/gd_cnv_calls.tsv.gz \
    --truth-table truth.tsv \
    --output-dir results/
```

The truth table must be a BED-like TSV with columns `chrom`, `start`, `end`, `name`,
`svtype`, `samples`, `NAHR_GD`, `NAHR_GD_atypical`. Only canonical NAHR-mediated records
(`NAHR_GD = True`, `NAHR_GD_atypical = False`) are evaluated; atypical breakpoint records
are automatically dropped.

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
| `call.py` | CNV calling (MLP and Viterbi strategies); `call` subcommand |
| `output.py` | Posterior table writers and ploidy estimation |
| `plot.py` | Depth and CN-state visualisation; `plot` subcommand |
| `annotations.py` | Gene and genomic feature annotation overlays for plots |
| `eval.py` | Sensitivity/precision evaluation; `eval` subcommand |

---

## Development

```bash
pip install -e "src/gatk-sv-gd[dev]"
pytest src/gatk-sv-gd/
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
