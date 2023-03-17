---
title: ReGenotypeCNVs
description: Regenotype CNVs
sidebar_position: 10
slug: rgcnvs
---

Re-genotypes probable mosaic variants across multiple batches.

### Prerequisites
- Genotype batch

### Inputs

- Per-sample median coverage estimates (GatherBatchEvidence)
- Pre-genotyping depth VCFs (FilterBatch)
- Batch PED files (FilterBatch)
- Cohort depth VCF (MergeBatchSites)
- Genotyped depth VCFs (GenotypeBatch)
- Genotyped depth RD cutoffs file (GenotypeBatch)

### Outputs

- Re-genotyped depth VCFs.
