# Long-Read SV Benchmarking Pipeline

This directory contains scripts for benchmarking structural variant (SV) callsets from long-read sequencing against short-read reference callsets.

## Overview

The benchmarking pipeline consists of two main user-facing scripts:

1. **`process_samples.sh`** - Processes individual samples through the SV evaluation pipeline
2. **`analyze_multi_caller_support_master.sh`** - Analyzes multi-caller support patterns across samples

## Prerequisites

Ensure the following tools are installed and available in your PATH:
- `bcftools` (VCF manipulation)
- `tabix` (VCF indexing)
- `svtk` (SV toolkit for VCF2BED conversion)
- `bgzip` (compression)
- `bedtools` (BED file operations)
- `python3` (Python scripts)
- `R` (for genomic context annotation)

### Required Reference Files

The reference annotation files are located in `benchmark_scripts/input/`:
- `hg38.RM.sorted.merged.bed.gz` - RepeatMasker annotations
- `hg38.SD.sorted.merged.bed.gz` - Segmental Duplications
- `hg38.SR.sorted.merged.bed.gz` - Simple Repeats

## Main Scripts

### 1. `process_samples.sh` - Main Processing Pipeline

This is the primary script for processing multiple samples through the benchmarking pipeline.

#### Usage

```bash
./process_samples.sh --mapping MAPPING_FILE --base-dir BASE_DIR [OPTIONS]
```

#### Required Arguments

- `--mapping MAPPING_FILE`: Path to tab-delimited mapping file with columns `SR_ID` and `LR_ID`
- `--base-dir BASE_DIR`: Base directory containing your input data

#### Optional Arguments

- `--num-samples NUM_SAMPLES`: Number of samples to process (default: 3)
- `--lr-caller LR_CALLER`: Long-read caller name (default: sniffles)
- `--comparison-directions DIRS`: Comparison directions: SR, LR, or SR,LR (default: SR,LR)
- `--sr-limit SR_LIMIT`: Maximum number of SR variants for testing (default: no limit)
- `--lr-limit LR_LIMIT`: Maximum number of LR variants for testing (default: no limit)
- `--distance DISTANCE`: Maximum distance between breakpoints for matching (default: 50bp)
- `--sr-vcf-path SR_VCF_PATH`: Path to SR VCF relative to base-dir (default: all_batches.filter_genotypes.sanitized.135.overlaps_tagged.vcf.gz)
- `--lr-vcf-path LR_VCF_PATH`: Directory for LR VCFs relative to base-dir (default: {lr-caller}_vcfs)
- `--lr-vcf-pattern PATTERN`: Filename pattern for LR VCFs with placeholders {LR_ID} and {LR_CALLER}
- `--output-dir OUTPUT_DIR`: Output directory (default: {lr-caller}_results)

#### Process Types

Control which processing steps to run with `--processes` (comma-separated):
- `PROCESS_SR`: Process short-read data (VCF filtering, annotation, BED conversion)
- `PROCESS_LR`: Process long-read data (VCF filtering, annotation, BED conversion)
- `COMPARE_CALLSETS`: Run callset comparisons
- `ANALYZE_RATES`: Run match rate analysis (requires comparison results)

Default: All processes are run

#### Examples

Basic usage with defaults:
```bash
./process_samples.sh --mapping Mapping.tsv --base-dir .
```

Custom SR VCF and distance parameter:
```bash
./process_samples.sh --mapping Mapping.tsv --base-dir . \
    --sr-vcf-path custom_sr_callset.vcf.gz --distance 100
```

Custom LR caller with pattern:
```bash
./process_samples.sh --mapping Mapping.tsv --base-dir . \
    --lr-caller pbsv --lr-vcf-pattern '{LR_ID}.pbsv.vcf.gz'
```

Process only 5 samples, run only analysis:
```bash
./process_samples.sh --mapping Mapping.tsv --base-dir . \
    --num-samples 5 --processes ANALYZE_RATES
```

#### Input File Structure

**Mapping file format** (tab-delimited):
```
SR_ID	LR_ID
Sample1_SR	Sample1_LR
Sample2_SR	Sample2_LR
```

**Expected directory structure**:
```
base-dir/
├── {sr-vcf-path}                           # Short-read VCF file
├── {lr-vcf-path}/                          # Long-read VCF directory
│   ├── {LR_ID}.hg38.{caller}.sv.vcf.gz    # Long-read VCF files
│   └── ...
└── {output-dir}/                           # Output directory (created)
    ├── {SR_ID}/                            # Per-sample results
    │   ├── *.sr.* files                    # Short-read processing results
    │   ├── *.lr.* files                    # Long-read processing results
    │   └── *_query.bed files               # Comparison results
    └── ...
```

#### Output Structure

For each sample, the following files are created in `{output-dir}/{SR_ID}/`:

**Short-read processing:**
- `{SR_ID}.sr.vcf.gz` - Extracted SR sample VCF
- `{SR_ID}.sr.filtered.vcf.gz` - Filtered SR variants
- `{SR_ID}.sr.filtered.revised.vcf.gz` - SVTYPE-annotated variants
- `{SR_ID}.sr.bed` - SR variants in BED format
- `{SR_ID}.sr.annotated.bed` - Genomic context annotations
- `{SR_ID}.sr.filter_info.tsv` - Filter information

**Long-read processing:**
- `{SR_ID}.lr.vcf.gz` - Filtered LR variants
- `{SR_ID}.lr.bed` - LR variants in BED format
- `{SR_ID}.lr.annotated.bed` - Genomic context annotations
- `{SR_ID}.lr.filter_info.tsv` - Filter information

**Comparison results:**
- `{SR_ID}.lr_query.bed` - LR query vs SR reference comparison
- `{SR_ID}.sr_query.bed` - SR query vs LR reference comparison

**Analysis results:**
- `lr_query_*` files - LR query analysis results
- `sr_query_*` files - SR query analysis results

### 2. `analyze_multi_caller_support_master.sh` - Multi-Caller Analysis

Analyzes multi-caller support rates across multiple samples and generates summary visualizations.

#### Usage

```bash
./analyze_multi_caller_support_master.sh --mapping MAPPING_FILE --input-dir INPUT_DIR --n-matches N_MATCHES --svtypes SVTYPES --output-dir OUTPUT_DIR [OPTIONS]
```

#### Required Arguments

- `--mapping MAPPING_FILE`: Sample mapping file
- `--input-dir INPUT_DIR`: Directory containing analysis results from multiple callers
- `--n-matches N_MATCHES`: Number of matches to analyze
- `--svtypes SVTYPES`: SV types to analyze (comma-separated, e.g., DEL,DUP,INS,INV)
- `--output-dir OUTPUT_DIR`: Output directory for results

#### Optional Arguments

- `--lr-callers LR_CALLERS`: Comma-separated list of LR callers (default: sniffles,cutesv)

#### Input Structure

The input directory should contain subdirectories for each caller:
```
input-dir/
├── sniffles_results/
│   ├── {SR_ID}/
│   │   └── lr_query_* files
│   └── ...
├── cutesv_results/
│   ├── {SR_ID}/
│   │   └── lr_query_* files
│   └── ...
└── ...
```

#### Output Files

- `multi_caller_support_rates_violin_plot.png` - Violin plot visualization
- `multi_caller_support_rates_heatmap.tsv` - Summary heatmap data
- `combined_analysis_results.tsv` - Combined analysis data

## Analysis Methods

### Overlap Detection

The pipeline uses reciprocal overlap criteria:
- **Deletions/Duplications**: 50% reciprocal overlap
- **Insertions**: Position-based matching within specified distance
- **Inversions**: Breakpoint-based matching within specified distance

### Filter Categories

Variants are classified into filter categories:
- **PASS**: High-quality variants that pass all filters
- **MULTIALLELIC**: Multi-allelic variants
- **HIGH_SR_BACKGROUND**: High short-read background
- **OTHER**: Other filter categories

## Troubleshooting

### Common Issues

1. **Missing reference files**: Ensure all three reference files are present in `benchmark_scripts/input/`
2. **File permissions**: Run `chmod +x benchmark_scripts/*.sh` to make scripts executable
3. **Tool availability**: Verify all required tools are installed and in PATH
4. **File format**: Ensure VCF files are properly formatted and indexed
5. **Custom file paths**: Use `--sr-vcf-path` and `--lr-vcf-pattern` for non-standard file locations

### Error Messages

- `Error: Reference file (X) not found`: Check that reference files exist in `benchmark_scripts/input/`
- `Error: Input VCF file not found`: Verify the SR VCF path is correct relative to base-dir
- `Error: Long-read VCF file not found`: Check LR VCF path and pattern settings
- `bcftools: command not found`: Install bcftools or ensure it's in your PATH

For additional help, run `./process_samples.sh --help` for detailed usage information 