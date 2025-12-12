---
title: AnnotateVcf
description: Annotate VCF
sidebar_position: 21
slug: av
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

Adds annotations, such as the inferred function and allele frequencies of variants, to a VCF.

Annotations methods include:
* Functional annotation - The GATK tool [SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/13832752531355-SVAnnotate) 
  is used to annotate SVs with inferred functional consequence on protein-coding regions, regulatory regions such as 
  UTR and promoters, and other non-coding elements.
* Allele Frequency (`AF`) annotation - annotate SVs with their allele frequencies across all samples, and samples of 
  specific sex, as well as specific subpopulations.
* Allele Frequency annotation with external callset - annotate SVs with the allele frequencies of their overlapping SVs
  in another callset, eg. the gnomAD-SV reference callset.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  fg: FilterGenotypes
  avcf: AnnotateVcf
  fg --> avcf
  
  class avcf thisModule
  class fg inModules
```

### Inputs

#### `vcf`
Any SV VCF. Running on the [genotype filtered VCF](./fg#filtered_vcf) is recommended.

#### `prefix`
Prefix for the output VCF, such as the cohort name. May be alphanumeric with underscores.

#### `protein_coding_gtf`
Coding transcript definitions, see [here](/docs/resources#protein_coding_gtf).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `noncoding_bed`
Non-coding reference intervals, see [here](//docs/resources#noncoding_bed).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `promoter_window`
Promoter window size. See [here](https://gatk.broadinstitute.org/hc/en-us/articles/27007964610331-SVAnnotate#--promoter-window-length).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `max_breakend_as_cnv_length`
Max size for treating `BND` records as CNVs. See [here](https://gatk.broadinstitute.org/hc/en-us/articles/27007964610331-SVAnnotate#--max-breakend-as-cnv-length).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `svannotate_additional_args`
Additional arguments for [GATK-SVAnnotate](https://gatk.broadinstitute.org/hc/en-us/articles/27007964610331-SVAnnotate).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `sample_pop_assignments`
Two-column file with sample ID & population assignment. "." for population will ignore the sample. If provided, 
annotates population-specific allele frequencies.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `sample_keep_list`
If provided, subset samples to this list in the output VCF.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format). 
If provided, sex-specific allele frequencies will be annotated.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `par_bed`
Pseudo-autosomal region (PAR) bed file. If provided, variants overlapping PARs will be annotated with the `PAR` field.

#### `sv_per_shard`
Shard size for parallel processing. Decreasing this may help if the workflow is running too slowly.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `external_af_ref_bed`
Reference SV set (see [here](/docs/resources#external_af_ref_bed)). If provided, annotates variants with allele frequencies 
from the reference population.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `external_af_ref_prefix`
External `AF` annotation prefix. Required if providing [external_af_ref_bed](#optional-external_af_ref_bed).

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `external_af_population`
Population names in the external SV reference set, e.g. "ALL", "AFR", "AMR", "EAS", "EUR". Required if providing 
[external_af_ref_bed](#optional-external_af_ref_bed) and must match the populations in the bed file.

### Outputs

#### `annotated_vcf`
Output VCF.
