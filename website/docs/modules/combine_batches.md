---
title: CombineBatches
description: Cross-batch variant clustering
sidebar_position: 11
slug: cmb
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/CombineBatches.wdl)

Merges variants across multiple batches. Variant merging uses similar methods and criteria as in [ClusterBatch](./cb), 
but in addition requires samples genotyped as non-reference to match sufficiently.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
    
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  rgc: RegenotypeCNVs
  cb: CombineBatches
  rcv: ResolveComplexVariants
  rgc --> cb
  cb --> rcv
  
  class cb thisModule
  class rgc inModules
  class rcv outModules
```

### Inputs

:::info
All array inputs of batch data must match in order. For example, the order of the `batches` array should match that of
`pesr_vcfs`, `depth_vcfs`, etc.
:::

#### `cohort_name`
Cohort name. The guidelines outlined in the [sample ID requirements](/docs/gs/inputs#sampleids) section apply here.

#### `batches`
Array of batch identifiers. Should match the name used in [GatherBatchEvidence](./gbe#batch). Order must match that of [depth_vcfs](#depth_vcfs).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `merge_vcfs`
Default: `false`. If true, merge contig-sharded VCFs into one genome-wide VCF. This may be used for convenience but cannot be used with 
downstream workflows.

#### `pesr_vcfs`
Array of genotyped depth caller variants for all batches, generated in [GenotypeBatch](./gb#genotyped_depth_vcf).

#### `depth_vcfs`
Array of re-genotyped depth caller variants for all batches, generated in [RegenotypeCNVs](./rgcnvs#regenotyped_depth_vcfs).

#### `raw_sr_bothside_pass_files`
Array of variant lists with bothside SR support for all batches, generated in [GenotypeBatch](./gb#sr_bothside_pass).

#### `raw_sr_background_fail_files`
Array of variant lists with low SR signal-to-noise ratio for all batches, generated in [GenotypeBatch](./gb#sr_background_fail).

#### `localize_shard_size`
Shard size for parallel computations. Decreasing this parameter may help reduce run time.

#### `min_sr_background_fail_batches`
Threshold fraction of batches with high SR background for a given variant required in order to assign this 
`HIGH_SR_BACKGROUND` flag. Most users should leave this at the default value.

### Outputs

#### `combined_vcfs`
Array of contig-sharded VCFs of combined PE/SR and depth calls.

#### `cluster_bothside_pass_lists`
Array of contig-sharded bothside SR support variant lists.

#### `cluster_background_fail_lists`
Array of contig-sharded high SR background variant lists.

#### `combine_batches_merged_vcf`
Genome-wide VCF of combined PE/SR and depth calls. Only generated if using [merge_vcfs](#optional--merge_vcfs).
