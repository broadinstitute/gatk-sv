---
title: GenotypeBatch
description: Genotype Batch
sidebar_position: 9
slug: gb
---

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/GenotypeBatch.wdl)

Genotypes a batch of samples across all variants in the cohort. Note that while the preceding step 
[MergeBatchSites](./msites) is a "cohort-level" module, genotyping is performed on one batch of samples at a time.

In brief, genotyping is performed by first training variant metric cutoffs on sites with clear evidence signatures, 
and then genotypes and genotype qualities are assigned based on parametric models tuned with these cutoffs. This is 
performed separately for PE/SR calls and depth-based calls.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  mbs: MergeBatchSites
  gb: GenotypeBatch
  rgc: RegenotypeCNVs
  mbs --> gb
  gb --> rgc
  
  class gb thisModule
  class mbs inModules
  class rgc outModules
```

### Inputs

:::info
A number of inputs to this module are only used in single-sample mode and therefore omitted here. In addition, some 
inputs marked as optional in the WDL are required for joint calling.
:::

#### `batch`
An identifier for the batch. Should match the name used in [GatherBatchEvidence](./gbe#batch).

#### `batch_pesr_vcf`
Batch PE/SR caller variants after filtering, generated in [FilterBatch](./fb#filtered_pesr_vcf).

#### `batch_depth_vcf`
Batch depth caller variants after filtering, generated in [FilterBatch](./fb#filtered_depth_vcf).

#### `cohort_pesr_vcf`
Merged PE/SR caller variants for the cohort, generated in [MergeBatchSites](./msites#cohort_pesr_vcf).

#### `cohort_depth_vcf`
Merged depth caller variants for the cohort, generated in [MergeBatchSites](./msites#cohort_depth_vcf).

#### `n_per_split`
Records per shard when scattering variants. Decrease to increase parallelism if the workflow is running slowly.

#### `coveragefile`
Merged RD evidence file from [GatherBatchEvidence](./gbe#merged_bincov).

#### `medianfile`
Merged median coverage table from [GatherBatchEvidence](./gbe#median_cov).

#### `rf_cutoffs`
Genotyping cutoffs trained with the random forest filtering model from [FilterBatch](./fb#cutoffs).

#### `seed_cutoffs`
See [here](/docs/resources#seed_cutoffs).

#### `n_RD_genotype_bins`
Number of depth genotyping bins. Most users should leave this at the default value.

#### `discfile`
Merged PE evidence file from [GatherBatchEvidence](./gbe#merged_pe).

#### `reference_build`
Reference build version. Only "hg38" is supported.

#### `sr_median_hom_ins`
Median normalized split read counts of homozygous insertions. Most users should leave this at the default value.

#### `sr_hom_cutoff_multiplier`
Cutoff multiplier for split read counts of homozygous insertions. Most users should leave this at the default value.

### Outputs

#### `sr_bothside_pass`
List of variant IDs with split reads found on both sides of the breakpoint.

#### `sr_background_fail`
List of variant IDs exhibiting low signal-to-noise ratio for split read evidence.

#### `trained_PE_metrics`
PE evidence genotyping metrics file.

#### `trained_SR_metrics`
SR evidence genotyping metrics file.

#### `trained_genotype_*_*_sepcutoff`
Trained genotyping cutoffs for variants called by PESR or depth when supported by PESR or depth evidence (4 total).

#### `genotyped_depth_vcf`
Genotyped depth call VCF.

#### `genotyped_pesr_vcf`
Genotyped PE/SR call VCF.

#### `regeno_coverage_medians`
Coverage metrics for downstream CNV regenotyping.
