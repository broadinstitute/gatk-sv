---
title: EvidenceQC
description: Evidence QC
sidebar_position: 2
slug: eqc
---

import { Highlight, HighlightOptionalArg } from "../../src/components/highlight.js"

Runs ploidy estimation, dosage scoring, and optionally VCF QC. 
The results from this module can be used for QC and batching.

For large cohorts, this workflow can be run on arbitrary cohort 
partitions of up to about 500 samples. Afterwards, we recommend 
using the results to divide samples into smaller batches (~100-500 samples) 
with ~1:1 male:female ratio. Refer to the [Batching](/docs/run/joint#batching) section 
for further guidance on creating batches.

We also recommend using sex assignments generated from the ploidy 
estimates and incorporating them into the PED file, with sex = 0 for sex aneuploidies.

The following diagram illustrates the upstream and downstream workflows of the `EvidenceQC` workflow 
in the recommended invocation order. You may refer to 
[this diagram](https://github.com/broadinstitute/gatk-sv/blob/main/terra_pipeline_diagram.jpg) 
for the overall recommended invocation order.

<br/>

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#00509d,color:#caf0f8
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gse: GatherSampleEvidence
  eqc: EvidenceQC
  batching: Batching, sample QC, and sex assignment
  
  gse --> eqc
  eqc --> batching
  
  class eqc thisModule
  class gse inModules
  class batching outModules
```

<br/>


### Inputs

- Read count files (GatherSampleEvidence)
- (Optional) SV call VCFs (GatherSampleEvidence)

### Outputs

- Per-sample dosage scores with plots
- Median coverage per sample
- Ploidy estimates, sex assignments, with plots
- (Optional) Outlier samples detected by call counts

## Preliminary Sample QC

The purpose of sample filtering at this stage after EvidenceQC is to 
prevent very poor quality samples from interfering with the results for 
the rest of the callset. In general, samples that are borderline are 
okay to leave in, but you should choose filtering thresholds to suit 
the needs of your cohort and study. There will be future opportunities 
(as part of FilterBatch) for filtering before the joint genotyping 
stage if necessary. Here are a few of the basic QC checks that we recommend:

- Look at the X and Y ploidy plots, and check that sex assignments 
  match your expectations. If there are discrepancies, check for 
  sample swaps and update your PED file before proceeding.

- Look at the dosage score (WGD) distribution and check that 
  it is centered around 0 (the distribution of WGD for PCR- 
  samples is expected to be slightly lower than 0, and the distribution 
  of WGD for PCR+ samples is expected to be slightly greater than 0. 
  Refer to the gnomAD-SV paper for more information on WGD score). 
  Optionally filter outliers.

- Look at the low outliers for each SV caller (samples with 
  much lower than typical numbers of SV calls per contig for 
  each caller). An empty low outlier file means there were 
  no outliers below the median and no filtering is necessary. 
  Check that no samples had zero calls.

- Look at the high outliers for each SV caller and optionally 
  filter outliers; samples with many more SV calls than average may be poor quality.

- Remove samples with autosomal aneuploidies based on 
  the per-batch binned coverage plots of each chromosome.
