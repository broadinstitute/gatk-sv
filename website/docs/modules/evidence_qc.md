---
title: EvidenceQC
description: Evidence QC
sidebar_position: 2
slug: eqc
---

import { Highlight, HighlightOptionalArg } from "../../src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/EvidenceQC.wdl)

Runs ploidy estimation, dosage scoring, and optionally VCF QC. 
The results from this module can be used for QC and batching.

For large cohorts, this workflow can be run on arbitrary cohort 
partitions of up to about 500 samples. Afterward, we recommend 
using the results to divide samples into smaller batches (~100-500 samples) 
with ~1:1 male:female ratio. Refer to the [Batching](/docs/execution/joint#batching) section 
for further guidance on creating batches.

We also recommend using sex assignments generated from the ploidy 
estimates and incorporating them into the PED file, with sex = 0 for sex aneuploidies.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
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


### Inputs

All array inputs of sample data must match in order. For example, the order of the `samples` array should match that 
of the `counts` array.

#### `batch`
A name for the batch of samples being run. Can be alphanumeric with underscores.

#### `samples`
Sample IDs. Must match those used in [GatherSampleEvidence](./gse#outputs).

#### `counts`
Binned read counts (`.counts.tsv.gz`) from [GatherSampleEvidence](./gse#outputs)

#### `*_vcfs`
Raw SV call VCFs (`.vcf.gz`) from [GatherSampleEvidence](./gse#outputs). May be omitted in case a caller was not run.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `run_vcf_qc`
Default: `false`. Run raw call VCF QC analysis.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `run_ploidy`
Default: `true`. Run ploidy estimation.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `melt_insert_size`
Mean insert size for each sample. Produces QC tables and plots if available.


### Outputs

#### `WGD_*`
Per-sample whole-genome dosage scores with plots

#### `bincov_median`
Median coverage per sample

#### `bincov_matrix`
Binned read depth matrix for the submitted batch

#### `ploidy_*`
Ploidy estimates, sex assignments, with plots

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `*_qc_low`, `*_qc_high`
Outlier samples detected by call counts.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `qc_table`
QC summary table. Enable with [run_ploidy](#optional-run_ploidy).

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
