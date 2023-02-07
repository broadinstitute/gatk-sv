---
title: Overview
description: Overview of the constituting components
sidebar_position: 0
---

The pipeline is written in [Workflow Description Language (WDL)](https://openwdl.org),
consisting of multiple modules as listed in the following. 

- **Gather Sample Evidence** SV evidence collection, including calls from a configurable set of 
  algorithms (Manta, MELT, and Wham), read depth (RD), split read positions (SR), 
  and discordant pair positions (PE).

- **Evidence QC** Dosage bias scoring and ploidy estimation.

- **Train gCNV**

- **Gather Batch Evidence** Copy number variant calling using 
  `cn.MOPS` and `GATK gCNV`; B-allele frequency (BAF) generation; 
   call and evidence aggregation.

- **Cluster Batch** Variant clustering

- **Generate Batch Metrics** Variant filtering metric generation

- **Filter Batch** Variant filtering; outlier exclusion

- **Genotype Batch** Genotyping

- **Make Cohort Vcf** Cross-batch integration; complex variant resolution and re-genotyping; vcf cleanup

- **Module 07** Downstream filtering, including minGQ, batch effect check, 
  outlier samples removal and final recalibration;

- **Annotate Vcf** Annotations, including functional annotation, 
  allele frequency (AF) annotation and AF annotation with external population callsets;

- **Module 09** Visualization, including scripts that generates IGV screenshots and rd plots.

- Additional modules to be added: de novo and mosaic scripts