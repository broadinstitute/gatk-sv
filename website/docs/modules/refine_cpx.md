---
title: RefineComplexVariants
description: Refines and filters complex variants
sidebar_position: 15
slug: refcv
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/RefineComplexVariants.wdl)

Refines complex SVs and translocations and filters based on discordant read pair and read depth evidence reassessment.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  cvcf: CleanVcf
  refcv: RefineComplexVariants
  svc: SVConcordance
  
  cvcf --> refcv
  refcv --> svc
  
  class refcv thisModule
  class cvcf inModules
  class svc outModules
```

### Inputs

:::info
All array inputs of batch data must match in order. For example, the order of the `batch_name_list` array should match 
that of `batch_sample_lists`, `PE_metrics`, etc.
:::

#### `vcf`
Input vcf, generated in [CleanVcf](./cvcf#cleaned_vcf).

#### `prefix`
Prefix for output VCF, such as the cohort name. May be alphanumeric with underscores.

#### `batch_name_list`
Array of batch names. These should be the same batch names used in [GatherBatchEvidence](./gbe#batch).

#### `batch_sample_lists`
Array of sample ID lists for all batches, generated in [FilterBatch](./fb#batch_samples_postoutlierexclusion). Order must match [batch_name_list](#batch_name_list).

#### `PE_metrics`
Array of PE metrics files for all batches, generated in [GatherBatchEvidence](./gbe#merged_pe). Order must match [batch_name_list](#batch_name_list).

#### `Depth_DEL_beds`, `Depth_DUP_beds`
Arrays of raw DEL and DUP depth calls for all batches, generated in [GatherBatchEvidence](./gbe#merged_dels-merged_dups). Order must match [batch_name_list](#batch_name_list).

#### `n_per_split`
Shard size for parallel computations. Decreasing this parameter may help reduce run time.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `min_pe_cpx`
Default: `3`. Minimum PE read count for complex variants (CPX).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `min_pe_ctx`
Default: `3`. Minimum PE read count for translocations (CTX).

### Outputs

#### `cpx_refined_vcf`
Output VCF.

#### `cpx_evidences`
Supplementary output table of complex variant evidence.
