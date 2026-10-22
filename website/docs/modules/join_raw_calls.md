---
title: JoinRawCalls
description: Clusters unfiltered variants across batches
sidebar_position: 17
slug: jrc
---

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/JoinRawCalls.wdl)

This module clusters raw unfiltered variants from [ClusterBatch](./cb) across all batches. Concordance between these 
genotypes and the joint call set usually can be indicative of variant quality and is used downstream for genotype 
filtering.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
    
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  cb: ClusterBatch
  jrc: JoinRawCalls
  svc: SVConcordance
  cb --> jrc
  jrc --> svc
  
  class jrc thisModule
  class cb inModules
  class svc outModules
```

### Inputs

#### `prefix`
Prefix for the output VCF, such as the cohort name. The guidelines outlined in the 
[sample ID requirements](/docs/gs/inputs#sampleids) section apply here.

#### `clustered_*_vcfs`
Array of clustered VCFs for all batches, generated in [ClusterBatch](./cb#clustered__vcf).

#### `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format).

### Outputs

#### `joined_raw_calls_vcf`
VCF containing all raw calls in the cohort.

#### `ploidy_table`
TSV of contig ploidies for all samples, assuming diploid autosome and sex assignments from the [ped_file](#ped_file).
