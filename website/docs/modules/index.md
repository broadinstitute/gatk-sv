---
title: Overview
description: Overview of the constituting components
sidebar_position: 0
---

The pipeline is written in [Workflow Description Language (WDL)](https://openwdl.org) and consists of multiple modules 
to be executed in the following order for joint calling. Note that while the single-sample mode calls many of these 
modules, it is implemented as a separate runnable workflow.

Below is a summary of the core modules in GATK-SV. Detailed information is documented in the following sections.
For details on running these modules for joint calling in Terra, refer to the 
[Workflow instructions](/docs/execution/joint#instructions) section.

- **GatherSampleEvidence** SV evidence collection, including calls from a configurable set of 
  algorithms (Manta, Scramble, and Wham), read depth (RD), split read positions (SR), 
  and discordant pair positions (PE).

- **EvidenceQC** Dosage bias scoring and ploidy estimation.

- **GatherBatchEvidence** Copy number variant calling using 
  `cn.MOPS` and `GATK gCNV`; B-allele frequency (BAF) generation; 
   call and evidence aggregation.

- **ClusterBatch** Intra-batch variant clustering

- **GenerateBatchMetrics** Variant filtering metric generation

- **FilterBatch** Variant filtering; outlier exclusion

- **GenotypeBatch** Joint genotyping

- **MakeCohortVcf** Cross-batch clustering; complex variant resolution and re-genotyping; VCF cleanup

- **AnnotateVCF** Annotations, including functional annotation, 
  allele frequency (AF) annotation and AF annotation with external population callsets
