---
title: RegenotypeCNVs
description: Regenotype CNVs
sidebar_position: 10
slug: rgcnvs
---

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/RegenotypeCNVs.wdl)

Re-genotypes probable mosaic variants across multiple batches. This is a "cohort-level" workflow that operates on 
all batches.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gb: GenotypeBatch
  rgc: RegenotypeCNVs
  cb: CombineBatches
  gb --> rgc
  rgc --> cb
  
  class rgc thisModule
  class gb inModules
  class cb outModules
```

### Inputs

:::info
All array inputs of batch data must match in order. For example, the order of the `batches` array should match that of
`depth_vcfs`, `batch_depth_vcfs`, etc.
:::

#### `depth_vcfs`
Array of genotyped depth caller variants for all batches, generated in [GenotypeBatch](./gb#genotyped_depth_vcf).

#### `cohort_depth_vcf`
Merged depth caller variants for the cohort, generated in [MergeBatchSites](./msites#cohort_depth_vcf).

#### `batch_depth_vcfs`
Array of filtered depth caller variants for all batches, generated in [FilterBatch](./fb#filtered_depth_vcf). Order must match that of [depth_vcfs](#depth_vcfs).

#### `coveragefiles`
Array of merged RD evidence files for all batches from [GatherBatchEvidence](./gbe#merged_bincov). Order must match that of [depth_vcfs](#depth_vcfs).

#### `medianfiles`
Array of median coverage tables for all batches from [GatherBatchEvidence](./gbe#median_cov). Order must match that of [depth_vcfs](#depth_vcfs).

#### `RD_depth_sepcutoffs`
Array of "depth_depth" genotype cutoff files (depth evidence for depth-based calls) generated in 
[GenotypeBatch](./gb#trained_genotype___sepcutoff). Order must match that of [depth_vcfs](#depth_vcfs).

#### `n_per_split`
Records per shard when scattering variants. Decrease to increase parallelism if the workflow is running slowly.

#### `n_RD_genotype_bins`
Number of depth genotyping bins. Most users should leave this at the default value.

#### `batches`
Array of batch identifiers. Should match the name used in [GatherBatchEvidence](./gbe#batch). Order must match that of [depth_vcfs](#depth_vcfs).

#### `cohort`
Cohort name. May be alphanumeric with underscores.

#### `regeno_coverage_medians`
Array of regenotyping metrics generated in [GenotypeBatch](./gb#regeno_coverage_medians).

### Outputs

#### `regenotyped_depth_vcfs`
Array of batch depth VCFs after regenotyping.
