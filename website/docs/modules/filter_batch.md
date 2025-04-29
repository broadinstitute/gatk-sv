---
title: FilterBatch
description: Filter Batch
sidebar_position: 7
slug: fb
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/FilterBatch.wdl)

Filters poor quality variants and outlier samples. 
This workflow can be run all at once with the top-level WDL, 
or it can be run in two steps to enable tuning of outlier 
filtration cutoffs. The two subworkflows are:

1. [FilterBatchSites](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/FilterBatchSites.wdl): Per-batch variant filtration. 
    Visualize filtered SV counts per sample per type to help choose an IQR cutoff for outlier sample filtering, and preview
    outlier samples for a given cutoff.

2. [FilterBatchSamples](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/FilterBatchSamples.wdl): Per-batch outlier sample filtration; 
   provide an appropriate [outlier_cutoff_nIQR](#outlier_cutoff_niqr) based on the 
   SV count plots and outlier previews from step 2. Note 
   that not removing high outliers can result in increased 
   compute cost and a higher false positive rate in later steps.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gbm: GenerateBatchMetrics
  fb: FilterBatch
  mbs: MergeBatchSites
  gbm --> fb
  fb --> mbs
  
  class fb thisModule
  class gbm inModules
  class mbs outModules
```

### Inputs

#### `batch`
An identifier for the batch. Should match the name used in [GatherBatchEvidence](./gbe#batch).

#### `*_vcf`
Clustered VCFs from [ClusterBatch](./cb#clustered__vcf)

#### `evidence_metrics`
Metrics table [GenerateBatchMetrics](./gbm#metrics)

#### `evidence_metrics_common`
Common variant metrics table [GenerateBatchMetrics](./gbm#metrics_common)

#### `outlier_cutoff_nIQR`
Defines outlier sample cutoffs based on variant counts. Samples deviating from the batch median count by more than 
the given multiple of the interquartile range are hard filtered from the VCF. Recommended range is between `3` and `9`
depending on desired sensitivity (higher is less stringent), or disable with `10000`.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `outlier_cutoff_table`
A cutoff table to set permissible nIQR ranges for each SVTYPE. If provided, overrides `outlier_cutoff_nIQR`. Expected 
columns are: `algorithm`, `svtype`, `lower_cuff`, `higher_cff`. See the `outlier_cutoff_table` resource in 
[this json](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/ref_panel_1kg.json) for an example table.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `adjudicate_cutoffs`
A pre-computed table of random forest cutoffs that override the table that is generated in the workflow. It ensures
that the workflow does not re-run the random forest model training if cutoffs are already available, or if 
setting these in a manual manner is desired.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `adjudicate_scores`
A pre-computed table of variant scores that override the table that is generated in the workflow. It ensures
that the workflow does not re-run the random forest model evaluation if scores are already available, or if 
setting these in a manual manner is desired.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `adjudicate_rf_files`
A pre-computed series of interemediete random forest files that override what is generated in the workflow. 
It ensures that the workflow does not re-save these intermediete files if they are already available.

### Outputs

#### `filtered_depth_vcf`
Depth-based CNV caller VCFs after variant and sample filtering.

#### `filtered_pesr_vcf`
PE/SR (non-depth) caller VCFs after variant and sample filtering.

#### `cutoffs`
Variant metric cutoffs for genotyping.

#### `sv_counts`
Array of TSVs containing SV counts for each sample, i.e. `sample-svtype-count` triplets. Each file corresponds to
a different SV caller.

#### `sv_count_plots`
Array of images plotting SV counts stratified by SV type. Each file corresponds to a different SV caller.

#### `outlier_samples_excluded`
Array of sample IDs excluded by outlier analysis.

#### `outlier_samples_excluded_file`
Text file of sample IDs excluded by outlier analysis.

#### `batch_samples_postOutlierExclusion`
Array of remaining sample IDs after outlier exclusion.

#### `batch_samples_postOutlierExclusion_file`
Text file of remaining sample IDs after outlier exclusion.
