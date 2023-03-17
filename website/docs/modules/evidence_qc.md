---
title: EvidenceQC
description: Evidence QC
sidebar_position: 2
slug: eqc
---

Runs ploidy estimation, dosage scoring, and optionally VCF QC. 
The results from this module can be used for QC and batching.

For large cohorts, this workflow can be run on arbitrary cohort 
partitions of up to about 500 samples. Afterwards, we recommend 
using the results to divide samples into smaller batches (~100-500 samples) 
with ~1:1 male:female ratio. Refer to the [Batching](/docs/run/batching) section 
for further guidance on creating batches.

We also recommend using sex assignments generated from the ploidy 
estimates and incorporating them into the PED file, with sex = 0 for sex aneuploidies.

### Prerequisites

- [Gather Sample Evidence](./gse)

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
