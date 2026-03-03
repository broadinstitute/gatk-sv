# gatk-sv-gd

Genomic Disorder (GD) CNV detection from binned read counts.

Part of the [GATK-SV](https://github.com/broadinstitute/gatk-sv) pipeline.

## Installation

```bash
pip install -e src/gatk-sv-gd
```

Or with conda (after creating the environment):

```bash
conda activate <env>
pip install -e src/gatk-sv-gd
```

## Usage

The package provides a single CLI command `gatk-sv-gd` with four subcommands:

### 1. Infer — Bayesian CNV inference

```bash
gatk-sv-gd infer \
    --input counts.tsv \
    --gd-table gd_table.tsv \
    --output-dir results/
```

### 2. Call — Call CNVs from posteriors

```bash
gatk-sv-gd call \
    --cn-posteriors cn_posteriors.tsv.gz \
    --bin-mappings bin_mappings.tsv.gz \
    --gd-table gd_table.tsv \
    --output-dir results/
```

### 3. Plot — Generate visualisation plots

```bash
gatk-sv-gd plot \
    --calls gd_cnv_calls.tsv.gz \
    --cn-posteriors cn_posteriors.tsv.gz \
    --gd-table gd_table.tsv \
    --output-dir plots/
```

### 4. Eval — Evaluate against truth set

```bash
gatk-sv-gd eval \
    --calls gd_cnv_calls.tsv.gz \
    --truth-table truth.tsv \
    --output-dir results/
```

Run `gatk-sv-gd <subcommand> --help` for full option details.

## Development

```bash
pip install -e "src/gatk-sv-gd[dev]"
```
