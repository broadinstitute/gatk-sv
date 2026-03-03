# gatk-sv-ploidy

Whole-genome aneuploidy detection from binned read counts using a hierarchical
Bayesian model (Pyro).

## Installation

```bash
pip install -e .
```

## Usage

```bash
gatk-sv-ploidy <subcommand> [options]
```

### Subcommands

| Subcommand   | Description |
|-------------|-------------|
| `preprocess` | Read and normalise depth data, filter low-quality bins |
| `infer`      | Train Bayesian model and run discrete CN inference |
| `call`       | Assign sex karyotype and aneuploidy type per sample |
| `plot`       | Generate diagnostic and summary plots |
| `eval`       | Evaluate predictions against a truth set |

Run `gatk-sv-ploidy <subcommand> --help` for subcommand-specific options.

## Typical workflow

```bash
# 1. Preprocess raw depth data
gatk-sv-ploidy preprocess -i depth.tsv -o out/

# 2. Run Bayesian inference
gatk-sv-ploidy infer -i out/preprocessed_depth.tsv -o out/

# 3. Call sex and aneuploidy types
gatk-sv-ploidy call -c out/chromosome_stats.tsv -o out/

# 4. Generate plots
gatk-sv-ploidy plot -c out/chromosome_stats.tsv -o out/ \
    --bin-stats out/bin_stats.tsv.gz --training-loss out/training_loss.tsv

# 5. Evaluate against truth
gatk-sv-ploidy eval -p out/aneuploidy_type_predictions.tsv \
    -t truth.json -o out/
```
