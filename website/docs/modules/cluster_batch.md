---
title: ClusterBatch
description: Cluster Batch
sidebar_position: 5
slug: cb
---

import { Highlight, HighlightOptionalArg } from "../../src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/ClusterBatch.wdl)

Clusters SV calls across a batch.

The following diagram illustrates the downstream workflows of the `GatherBatchEvidence` workflow
in the recommended invocation order.

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#00509d,color:#caf0f8
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gbe: GatherBatchEvidence
  cb: ClusterBatch
  gbm: GenerateBatchMetrics
  gbe --> cb
  cb --> gbm
  
  class cb thisModule
  class gbe inModules
  class gbm outModules
```

## Inputs

#### `batch`
An identifier for the batch. Should match the name used in [GatherBatchEvidence](./gbe#batch).

#### `*_vcf_tar`
Standardized VCF tarballs from [GatherBatchEvidence](./gbe#std__vcf_tar)

#### `del_bed`, `dup_bed`
Merged CNV call files (`.bed.gz`) from [GatherBatchEvidence](./gbe#merged_dels-merged_dels)

#### `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `N_IQR_cutoff_plotting`
If provided, plot SV counts per sample. This number is used as the cutoff of interquartile range multiples for flagging 
outlier samples. Example value: 4.

## Outputs

#### `clustered_*_vcf`
Clustered variants for each caller (`depth` corresponds to depth-based CNV callers `cnMOPS` and `GATK-gCNV`) in VCF format.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `clustered_sv_counts`, `clustered_sv_count_plots`, `clustered_outlier_samples_preview`, `clustered_outlier_samples_with_reason`, `clustered_num_outlier_samples`
SV count QC tables and plots. Enable by providing [N_IQR_cutoff_plotting](#optional--n_iqr_cutoff_plotting)
