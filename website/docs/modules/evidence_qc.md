---
title: EvidenceQC
description: Evidence QC
sidebar_position: 2
slug: eqc
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/EvidenceQC.wdl)

Runs ploidy estimation, dosage scoring, and optionally VCF QC. 
The results from this module can be used for [QC](#preliminary-sample-qc) and [batching](#batching).

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

### Preliminary Sample QC

The purpose of sample filtering at this stage after EvidenceQC is to
prevent very poor quality samples from interfering with the results for
the rest of the callset. In general, samples that are borderline are
okay to leave in, but you should choose filtering thresholds to suit
the needs of your cohort and study. There will be future opportunities
(as part of [FilterBatch](/docs/modules/fb)) for filtering before the joint genotyping
stage if necessary. Here are a few of the basic QC checks that we recommend:

- Chromosome X and Y ploidy plots: check that sex assignments
  match your expectations. If there are discrepancies, check for
  sample swaps and update your PED file before proceeding.

- Whole-genome dosage score (WGD): examine distribution and check that
  it is centered around 0 (the distribution of WGD for PCR-
  samples is expected to be slightly lower than 0, and the distribution
  of WGD for PCR+ samples is expected to be slightly greater than 0.
  Refer to the gnomAD-SV paper for more information on WGD score).
  Optionally filter outliers.

- Low outliers for each SV caller: these are samples with
  much lower than typical numbers of SV calls per contig for
  each caller. An empty low outlier file means there were
  no outliers below the median and no filtering is necessary.
  Check that no samples had zero calls.

- High outliers for each SV caller: optionally
  filter outliers; samples with many more SV calls than average may be poor quality.

- Remove samples with autosomal aneuploidies based on
  the per-batch binned coverage plots of each chromosome.

In the joint calling mode Terra workspace, we provide a Jupyter notebook `SampleQC.ipynb`
for sample QC and filtering.


### Batching

For larger cohorts, samples should be split up into batches of about 100-500
samples with similar characteristics. We recommend batching based on overall
coverage and dosage score (WGD), which is generated in EvidenceQC.
You may also wish to batch samples based on other characteristics that could
impact SV calling, such as mean insert size or PCR status.
An example batching process is outlined below:

1. Divide the cohort by chromosome X ploidy (less than 2, greater than or equal to 2)
   based on copy ratio estimates from EvidenceQC. In this way, males and females will be
   batched separately before being merged back together for batches with equal sex balance
2. Partition the samples by median coverage from EvidenceQC,
   grouping samples with similar median coverage together
3. Partition the samples further by dosage score (WGD) from
   EvidenceQC, grouping samples with similar WGD score together
4. Optionally, partition the samples further by mean insert size if available,
   grouping samples with similar mean insert size together
5. Merge corresponding male and female partitions together to generate
   roughly equally sized batches of 100-500 samples with roughly equal sex balance

In the joint calling mode Terra workspace, we provide a Jupyter notebook `Batching.ipynb`
for batch creation.


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
