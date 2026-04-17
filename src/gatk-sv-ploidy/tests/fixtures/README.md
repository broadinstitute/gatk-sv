## Real fixture tiers

This directory contains the data snapshots used by the rewritten `gatk-sv-ploidy` tests.

- `tiny/`: synthetic fixtures defined in `tests/conftest.py` for fast unit tests.
- `medium/`: downsampled real `preprocessed_depth.tsv` + `site_data.npz` used for the `<1 minute` integration path (`infer -> call -> eval`).
- `large/`: larger raw depth snapshot used for the `few minute` integration path (`preprocess -> infer`).

Each real fixture directory includes:

- `manifest.json`: provenance, selected samples, chromosome set, and retained size.
- `truth.json`: expected labels used by integration tests.
- `sex_truth.json`: expected sex labels used by integration tests.

## Regenerating fixtures

The real fixtures are derived from the example run at:

- `/Users/markw/Work/talkowski/sv-pipe-testing/mw_ploidy/gatk-sv-ploidy/site_depth`

Regenerate both tiers from `src/gatk-sv-ploidy` with:

```bash
PYTHONPATH=src /Users/markw/miniconda3/envs/gatk-sv/bin/python tests/fixtures/generate_real_fixtures.py
```

You can override the source run or output root:

```bash
PYTHONPATH=src /Users/markw/miniconda3/envs/gatk-sv/bin/python tests/fixtures/generate_real_fixtures.py \
  --source-run /path/to/site_depth \
  --output-root tests/fixtures
```

## Current downsampling targets

- `medium`: 6 samples, 69 bins across `chr13`, `chr18`, `chr21`, `chrX`, `chrY`, and at most 8 SNP sites per bin.
- `large`: 12 samples, raw depth snapshot retained for full-chromosome preprocessing, then tested with `--viable-only` to keep `chr13`, `chr18`, `chr21`, `chrX`, `chrY`.

If fixture contents change, rerun `pytest tests -q` from `src/gatk-sv-ploidy` to confirm the regenerated snapshots still satisfy the runtime and schema expectations.