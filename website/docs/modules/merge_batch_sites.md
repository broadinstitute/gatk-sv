---
title: MergeBatchSites
description: Merge Batch Sites
sidebar_position: 8
slug: msites
---

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/MergeBatchSites.wdl)

Merges variants across batches. Variants are merged only if the following attributes match exactly:

- Contig
- Start position
- End position (`END` field)
- SV type (`SVTYPE` field)
- SV length (`SVLEN` field, if available)
- Strandedness (`STRANDS` field, if available)
- Second contig (`CHR2` field, if available)
- Second end (`END2` field, if available)

This is a "cohort-level" workflow, meaning that is aggregates data across all batches. This is in contrast to all previous 
modules, which are sample- or batch-level. Note that this workflow should still be run on cohorts consisting of 
a single batch.

:::info
Terra users must configure a "sample_set_set" in their data table before running this module. See the [Execution 
section on MergeBatchSites](/docs/execution/joint#09-mergebatchsites) for futher instructions.
:::

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  fb: FilterBatch
  mbs: MergeBatchSites
  gb: GenotypeBatch
  fb --> mbs
  mbs --> gb
  
  class mbs thisModule
  class fb inModules
  class gb outModules
```

### Inputs

#### `cohort`
An identifier for the cohort. The guidelines outlined in the [sample ID requirements](/docs/gs/inputs#sampleids) 
section apply here.

#### `depth_vcfs`
Array of filtered depth VCFs across batches generated in [FilterBatch](./fb#filtered_depth_vcf).

#### `pesr_vcfs`
Array of filtered PE/SR VCFs across batches generated in [FilterBatch](./fb#filtered_pesr_vcf).

### Outputs

#### `cohort_pesr_vcf`
Merged PE/SR caller VCF.

#### `cohort_depth_vcf`
Merged depth caller VCF.
