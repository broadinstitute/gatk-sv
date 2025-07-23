---
title: TrainGCNV
description: Train gCNV
sidebar_position: 3
slug: gcnv
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/TrainGCNV.wdl)

[GATK-gCNV](https://www.nature.com/articles/s41588-023-01449-0)
is a method for detecting rare germline copy number variants (CNVs)
from short-read sequencing read-depth information.
The `TrainGCNV` module trains a gCNV model for use in the [GatherBatchEvidence](./gbe) workflow. 

The samples used for training should be homogeneous (concerning sequencing platform, 
coverage, library preparation, etc.) and similar to the samples on which the model will be applied.

For small, relatively homogeneous cohorts, a single gCNV model is usually sufficient. 
However, for larger cohorts, especially those with multiple data sources, 
we recommend training a separate model for each batch or group of batches (see 
[batching section](/docs/execution/joint/#batching) for details).
The model can be trained on all or a subset of the samples to which it will be applied. 
A subset of 100 randomly selected samples from the batch is a reasonable
input size for training the model; when the `n_samples_subsample` input is provided, 
the `TrainGCNV` workflow can automatically perform this random selection.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  batching: Batching, sample QC, and sex assignment
  t: TrainGCNV
  gbe: GatherBatchEvidence
  
  batching --> t
  t --> gbe 
  
  class t thisModule
  class batching inModules
  class gbe outModules
```

## Inputs

The majority of the optional inputs of the workflow map to the optional arguments of the 
tool the workflow uses, `GATK-GermlineCNVCaller`; hence, you may refer to the 
[documentation](https://gatk.broadinstitute.org/hc/en-us/articles/360040097712-GermlineCNVCaller) 
of the tool for a description on these optional inputs.  We recommend that most users use the defaults.

:::info
All array inputs of sample data must match in order. For example, the order of the `samples` array should match that 
of the `count_files` array.
:::

#### `samples`
Sample IDs

#### `count_files`
Per-sample binned read counts (`*.rd.txt.gz`) generated in the [GatherSampleEvidence](./gse#outputs) workflow.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `n_samples_subsample`, `sample_ids_training_subset`
Provide one of these inputs to subset the input batch. `n_samples_subsample` will randomly subset, while 
`sample_ids_training_subset` is for defining a predetermined subset. These options are provided for convenience in Terra.

## Outputs

#### `cohort_contig_ploidy_model_tar`
Contig ploidy model tarball.

#### `cohort_gcnv_model_tars`
CNV model tarballs scattered across genomic intervals.

#### `cohort_contig_ploidy_calls_tar`
Contig ploidy calls for the submitted batch.

#### `cohort_gcnv_calls_tars`
CNV call tarballs scattered by sample and genomic region prior to segmentation.

#### `cohort_genotyped_segments_vcfs`
Single-sample VCFs of CNV calls for the submitted batch.

#### `cohort_gcnv_tracking_tars`
Convergence tracking logs.

#### `cohort_genotyped_intervals_vcfs`
Single-sample VCFs for the submitted batch containing per-interval genotypes prior to segmentation.

#### `cohort_denoised_copy_ratios`
TSV files containing denoised copy ratios in each sample.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `annotated_intervals` 
The count files from [GatherSampleEvidence](./gse) with adjacent intervals combined into 
locus-sorted `DepthEvidence` files using `GATK CondenseDepthEvidence` tool, which are
annotated with GC content, mappability, and segmental-duplication content using 
[`GATK-AnnotateIntervals`](https://gatk.broadinstitute.org/hc/en-us/articles/360041416652-AnnotateIntervals)
tool. This output is generated if `do_explicit_gc_correction` is set to `True`. Disabled by default.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `filtered_intervals_cnv`, `filtered_intervals_ploidy`
Intervals of read count bins to be used for CNV and ploidy calling after filtering for problematic regions (e.g. 
high GC content). This output is generated if `filter_intervals` is set to `True`. Enabled by default.
