---
title: CleanVcf
description: VCF cleaning
sidebar_position: 14
slug: cvcf
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/CleanVcf.wdl)

Performs various VCF clean-up steps including:

- Adjusting genotypes on allosomal contigs
- Collapsing overlapping CNVs into multi-allelic CNVs
- Revising genotypes in overlapping CNVs
- Removing redundant CNVs
- Stitching large CNVs
- VCF formatting clean-up

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gcv: GenotypeComplexVariants
  cvcf: CleanVcf
  refcv: RefineComplexVariants
  
  gcv --> cvcf
  cvcf --> refcv
  
  class cvcf thisModule
  class gcv inModules
  class refcv outModules
```

### Inputs

#### `cohort_name`
Cohort name. The guidelines outlined in the [sample ID requirements](/docs/gs/inputs#sampleids) section apply here.

#### `complex_genotype_vcfs`
Array of contig-sharded VCFs containing genotyped complex variants, generated in [GenotypeComplexVariants](./gcv#complex_genotype_vcfs).

#### `complex_resolve_bothside_pass_list`
Array of variant lists with bothside SR support for all batches, generated in [ResolveComplexVariants](./rcv#complex_resolve_bothside_pass_list).

#### `complex_resolve_background_fail_list`
Array of variant lists with low SR signal-to-noise ratio for all batches, generated in [ResolveComplexVariants](./rcv#complex_resolve_background_fail_list).

#### `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format).

#### `max_shards_per_chrom_step1`, `min_records_per_shard_step1`, `samples_per_step2_shard`, `max_samples_per_shard_step3`, `clean_vcf1b_records_per_shard`, `clean_vcf5_records_per_shard`
These parameters control parallelism in scattered tasks. Please examine the 
[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/CleanVcf.wdl) to see how each is used.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg>  `outlier_samples_list`
Text file of samples IDs to exclude when identifying multi-allelic CNVs. Most users do not need this feature unless 
excessive multi-allelic CNVs driven by low-quality samples are observed.

### Outputs

#### `cleaned_vcf`
Genome-wide VCF of output.

