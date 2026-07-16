---
title: IntegrateGDVcf
description: Integrate genomic disorder CNV calls into the cohort VCF
sidebar_position: 22
slug: igdv
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/IntegrateGDVcf.wdl)

Integrates the genomic disorder (GD) CNV calls produced by [CallGenomicDisorderCNVs](./gd) into the cohort's filtered
VCF. The workflow is scattered across chromosomes: the combined GD calls and per-batch ploidy tables are subset per
contig and merged into the VCF using the `gatk-sv-gd integrate` CLI, then the per-contig VCFs are concatenated back
into a single cohort VCF. See the [gatk-sv-gd](https://github.com/broadinstitute/gatk-sv-gd) repository for more
information.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gd: CallGenomicDisorderCNVs
  igdv: IntegrateGDVcf
  avcf: AnnotateVcf
  gd --> igdv
  igdv --> avcf
  
  class igdv thisModule
  class gd inModules
  class avcf outModules
```

## Inputs

#### `vcf`
Filtered VCF of the pipeline (cohort-level), and its index `vcf_index`.

#### `prefix`
Prefix for the output VCF, such as the cohort name.

#### `gd_output_tarballs`
Array of GD analysis tarballs, one per batch, produced by [CallGenomicDisorderCNVs](./gd#gd_output_tarball).

#### `ploidy_tables`
Array of ploidy tables, one per batch, generated in [JoinRawCalls](./jrc#ploidy_table).

#### `gd_table`
Path to the genomic disorder regions table (TSV), the same file used in [CallGenomicDisorderCNVs](./gd#gd_table).

#### `par_bed`
Path to pseudoautosomal region intervals (BED).

#### `contig_list`
List of contigs (one per line) to scatter the integration over.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `integrate_args`
Free-form string arguments passed to the `gatk-sv-gd integrate` subcommand.

## Outputs

#### `integrate_gd_vcf`
Cohort VCF with genomic disorder CNV calls integrated.
