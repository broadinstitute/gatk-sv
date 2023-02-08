---
title: MakeCohortVcf
description: Make Cohort VCF
sidebar_position: 11
slug: cvcf
---

Combines variants across multiple batches, resolves complex variants, 
re-genotypes, and performs final VCF clean-up.

### Prerequisites

- GenotypeBatch
- (Optional) RegenotypeCNVs

### Inputs

- RD, PE and SR file URIs (GatherBatchEvidence)
- Batch filtered PED file URIs (FilterBatch)
- Genotyped PESR VCF URIs (GenotypeBatch)
- Genotyped depth VCF URIs (GenotypeBatch or RegenotypeCNVs)
- SR pass variant file URIs (GenotypeBatch)
- SR fail variant file URIs (GenotypeBatch)
- Genotyping cutoff file URIs (GenotypeBatch)
- Batch IDs
- Sample ID list URIs

### Outputs

- Finalized "cleaned" VCF and QC plots