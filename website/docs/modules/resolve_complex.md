---
title: ResolveComplexVariants
description: Complex SV discovery
sidebar_position: 12
slug: rcv
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/ResolveComplexVariants.wdl)

Identifies multi-breakpoint complex variants, which are annotated with the `CPX` value in the `SVTYPE` field. These 
variants are putative, as read depth evidence is not assessed at this stage.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  cb: CombineBatches
  rcv: ResolveComplexVariants
  gcv: GenotypeComplexVariants
  cb --> rcv
  rcv --> gcv
  
  class rcv thisModule
  class cb inModules
  class gcv outModules
```

### Inputs

:::info
Some inputs of batch data must match in order. Specifically, the order of the `disc_files` array should match that of
`rf_cutoff_files`.
:::

#### `cohort_name`
Cohort name. The guidelines outlined in the [sample ID requirements](/docs/gs/inputs#sampleids) section apply here.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `merge_vcfs`
Default: `false`. If true, merge contig-sharded VCFs into one genome-wide VCF. This may be used for convenience but cannot be used with
downstream workflows.

#### `cluster_vcfs`
Array of contig-sharded VCFs, generated in [CombineBatches](./cmb#combined_vcfs).

#### `cluster_bothside_pass_lists`
Array of variant lists with bothside SR support for all batches, generated in [CombineBatches](./cmb#cluster_bothside_pass_lists).

#### `cluster_background_fail_lists`
Array of variant lists with low SR signal-to-noise ratio for all batches, generated in [CombineBatches](./cmb#cluster_background_fail_lists).

#### `disc_files`
Array of PE evidence files for all batches from [GatherBatchEvidence](./gbe#merged_pe).

#### `rf_cutoffs`
Array of batch genotyping cutoff files trained with the random forest filtering model from [FilterBatch](./fb#cutoffs).
Must match the order of [disc_files](#disc_files).

### Outputs

#### `complex_resolve_vcfs`
Array of contig-sharded VCFs containing putative complex variants.

#### `complex_resolve_bothside_pass_list`
Array of contig-sharded bothside SR support variant lists.

#### `complex_resolve_background_fail_list`
Array of contig-sharded high SR background variant lists.

#### `breakpoint_overlap_dropped_record_vcfs`
Variants dropped due to exact overlap with another's breakpoint.

#### `complex_resolve_merged_vcf`
Genome-wide output VCF. Only generated if using [merge_vcfs](#optional--merge_vcfs).
