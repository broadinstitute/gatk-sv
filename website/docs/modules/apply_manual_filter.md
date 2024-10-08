---
title: ApplyManualVariantFilter
description: Complex SV genotyping
sidebar_position: 16
slug: amvf
---

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/ApplyManualVariantFilter.wdl)

This module hard-filters variants (dropping records) using [bcftools](https://github.com/samtools/bcftools). While the 
workflow is general-purpose, we recommend running it with default parameters to eliminate major sources of false 
positive variants:

1. Deletions called solely by `Wham`.
2. SVA MEIs called by `Scramble` with the `HIGH_SR_BACKGROUND` flag.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
    
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  refcv: RefineComplexVariants
  amvf: ApplyManualVariantFilter
  svc: SVConcordance
  refcv --> amvf
  amvf --> svc
  
  class amvf thisModule
  class refcv inModules
  class svc outModules
```

### Inputs

#### `prefix`
Prefix for the output VCF, such as the cohort name. May be alphanumeric with underscores.

#### `vcf`
Any VCF. Running on the [cleaned VCF](cvcf#cleaned_vcf) is recommended.

#### `filter_name`
A name for the filter, used for output file naming. May be alphanumeric with underscores.

#### `bcftools_filter`
[Bcftools EXPRESSION](https://samtools.github.io/bcftools/bcftools.html#expressions) to use for filtering. Variants 
matching this expression will be **excluded**, i.e. with the `-e` argument.

### Outputs

#### `manual_filtered_vcf`
Filtered VCF.
