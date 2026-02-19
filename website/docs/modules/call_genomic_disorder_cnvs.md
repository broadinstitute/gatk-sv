---
title: CallGenomicDisorderCNVs
description: Detect CNVs at known genomic disorder loci
sidebar_position: 21
slug: gd
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/CallGenomicDisorderCNVs.wdl)

Runs Bayesian copy-number variant inference at known genomic disorder loci using cohort-wide binned read-depth and
B-allele frequency data. This workflow operates as a standalone post-processing step after `GatherBatchEvidence`
completes, consuming the merged BAF and read-depth matrices. See the [gatk-sv-gd](https://github.com/broadinstitute/gatk-sv-gd) repository for more information.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gbe: GatherBatchEvidence
  gd: CallGenomicDisorderCNVs
  igdv: IntegrateGDVcf
  gbe --> gd
  gd --> igdv
  
  class gbe inModules
  class gd thisModule
  class igdv outModules
```

## Inputs

#### `batch`
An identifier for the batch; may only be alphanumeric with underscores.

#### `baf_matrix`
Path to the merged B-allele frequency matrix (TSV, gzipped) produced by `GatherBatchEvidence`. Must be indexed
(`.tbi` index file required at `<baf_matrix>.tbi`).

#### `high_res_rd_matrix`
Path to the merged high-resolution read-depth matrix (TSV, gzipped) produced by `GatherBatchEvidence`.

#### `gd_table`
Path to the genomic disorder regions table (TSV) defining known disorder loci with breakpoint annotations.

#### `segdup_bed`
Path to segmental duplication regions (BED) used for masking.

#### `centromere_bed`
Path to centromere regions (BED).

#### `acrocentric_arm_bed`
Path to acrocentric chromosome arm regions (BED).

#### `custom_mask_bed`
Path to custom masking regions (BED).

#### `hard_inclusion_bed`
Path to hard inclusion intervals (BED) that must be included in analysis regardless of other filters.

#### `par_bed`
Path to pseudoautosomal region intervals (BED).

#### `gaps_bed`
Path to assembly gap regions (BED).

#### `gtf`
Path to gene annotation file (GTF).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `flank_exclusion_intervals`
Array of BED files defining intervals to exclude from flank analysis. Typically set to the same file as `segdup_bed`.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `rebinned_interval_size`
Target interval size in base pairs for read-depth rebinning. Default: `10000`.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `truth_table`
Path to a truth set table for benchmarking.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `preprocess_args`, `infer_args`, `call_args`, `eval_args`, `plot_args`
Free-form string arguments passed to the respective `gatk-sv-gd` subcommands.

## Outputs

#### `gd_output_tarball`
Compressed tarball (`.tar.gz`) containing all GD analysis results, consumed by [IntegrateGDVcf](./igdv). See [gatk-sv-gd](https://github.com/broadinstitute/gatk-sv-gd) for details about the tarball's contents.
