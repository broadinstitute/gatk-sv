---
title: GenotypeBatch
description: Genotype Batch
sidebar_position: 9
slug: gb
---

Genotypes a batch of samples across unfiltered variants combined across all batches.

### Prerequisites
- Filter batch
- Merge Batch Sites

### Inputs

- Batch PESR and depth VCFs (FilterBatch)
- Cohort PESR and depth VCFs (MergeBatchSites)
- Batch read count, PE, and SR files (GatherBatchEvidence)

### Outputs

- Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
- Filtered depth-only call VCF with outlier samples excluded
- PED file with outlier samples excluded
- List of SR pass variants
- List of SR fail variants
- (Optional) Depth re-genotyping intervals list
