# module21_StatCollection

This folder contains scripts developed to extract basic statistics of variant counts and distribution from VCF files.

## Goal of `CountVariantsByAFCutoff.wdl`

`CountVariantsByAFCutoff.wdl` summarizes per-sample variant burden from an input list of VCF files (for example, one VCF per contig).

Given an `af_cutoff`, it reports a table with:

1. `sample_id`
2. `total_variant_count`: total non-reference variant calls for that sample across all provided VCFs
3. `af_lt_cutoff_variant_count`: count of non-reference variant calls with `INFO/AF < af_cutoff`

This is intended as a lightweight stats-collection utility for downstream QC and cohort-level distribution checks.

## Inputs

- `Array[File] vcfs`
- `Float af_cutoff`
- Optional `output_prefix`
- Optional `bcftools_docker`

## Example input JSON

See `CountVariantsByAFCutoff.example.json`.
