---
title: GenotypeComplexVariants
description: Complex SV genotyping
sidebar_position: 13
slug: gcv
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/GenotypeComplexVariants.wdl)

Genotypes, filters, and classifies putative complex variants using depth evidence.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  rcv: ResolveComplexVariants
  gcv: GenotypeComplexVariants
  cvcf: CleanVcf
  rcv --> gcv
  gcv --> cvcf
  
  class gcv thisModule
  class rcv inModules
  class cvcf outModules
```

### Inputs

:::info
Some inputs of batch data must match in order. Specifically, the order of the `batches` array should match that of
`depth_vcfs`, `bincov_files`, `depth_gt_rd_sep_files`, and `median_coverage_files`.
:::

#### `cohort_name`
Cohort name. The guidelines outlined in the [sample ID requirements](/docs/gs/inputs#sampleids) section apply here.

#### `batches`
Array of batch identifiers. Should match the name used in [GatherBatchEvidence](./gbe#batch).

#### `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format).

#### `depth_vcfs`
Array of re-genotyped depth caller variants for all batches, generated in [RegenotypeCNVs](./rgcnvs#regenotyped_depth_vcfs).
Must match order of [batches](#batches).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `merge_vcfs`
Default: `false`. If true, merge contig-sharded VCFs into one genome-wide VCF. This may be used for convenience but cannot be used with
downstream workflows.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `localize_shard_size`
Default: `50000`. Shard size for parallel computations. Decreasing this parameter may help reduce run time.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `min_ddup_thresh`
Default: `5000`. Minimize size threshold used to classify a dispersed deletion.

#### `complex_resolve_vcfs`
Array of contig-sharded VCFs containing putative complex variants, generated in [ResolveComplexVariants](./rcv#complex_resolve_vcfs).

#### `bincov_files`
Array of RD evidence files for all batches from [GatherBatchEvidence](./gbe#counts). Must match order of [batches](#batches).

#### `depth_gt_rd_sep_files`
Array of "depth_depth" genotype cutoff files (depth evidence for depth-based calls) generated in
[GenotypeBatch](./gb#trained_genotype___sepcutoff). Order must match that of [batches](#batches).

#### `median_coverage_files`
Array of median coverage tables for all batches from [GatherBatchEvidence](./gbe#median_cov). Order must match that of [batches](#batches).

### Outputs

#### `complex_genotype_vcfs`
Array of contig-sharded VCFs containing fully resolved and genotyped complex variants.

#### `complex_genotype_merged_vcf`
Genome-wide output VCF. Only generated if using [merge_vcfs](#optional--merge_vcfs).
