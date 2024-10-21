---
title: Overview
description: Overview of the constituting components
sidebar_position: 0
---

The pipeline is written in [Workflow Description Language (WDL)](https://openwdl.org),
consisting of multiple modules to be executed in the following order. 

- **GatherSampleEvidence** SV evidence collection, including calls from a configurable set of 
  algorithms (Manta, MELT, and Wham), read depth (RD), split read positions (SR), 
  and discordant pair positions (PE).

- **EvidenceQC** Dosage bias scoring and ploidy estimation.

- **GatherBatchEvidence** Copy number variant calling using 
  `cn.MOPS` and `GATK gCNV`; B-allele frequency (BAF) generation; 
   call and evidence aggregation.

- **ClusterBatch** Variant clustering

- **GenerateBatchMetrics** Variant filtering metric generation

- **FilterBatch** Variant filtering; outlier exclusion

- **GenotypeBatch** Genotyping

- **MakeCohortVcf** Cross-batch integration; complex variant resolution and re-genotyping; vcf cleanup

- **Module 07 (in development)** Downstream filtering, including minGQ, batch effect check, 
  outlier samples removal and final recalibration;

- **AnnotateVCF** Annotations, including functional annotation, 
  allele frequency (AF) annotation and AF annotation with external population callsets;

- **Module 09 (in development)** Visualization, including scripts that generates IGV screenshots and rd plots.

- Additional modules to be added: de novo and mosaic scripts


## Pipeline Parameters

Several inputs are shared across different modules of the pipeline, which are explained in this section.

#### `ped_file`

A pedigree file describing the familial relationships between the samples in the cohort.
The file needs to be in the 
[PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).
Updated with [EvidenceQC](./eqc) sex assignments, including 
`sex = 0` for sex aneuploidies; 
genotypes on chrX and chrY for samples with `sex = 0` in the PED file will be set to 
`./.` and these samples will be excluded from sex-specific training steps.
