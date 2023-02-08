---
title: FilterBatch
description: Filter Batch
sidebar_position: 7
slug: fb
---

Filters poor quality variants and filters outlier samples. 
This workflow can be run all at once with the WDL at wdl/FilterBatch.wdl, 
or it can be run in three steps to enable tuning of outlier 
filtration cutoffs. The three subworkflows are:

1. FilterBatchSites: Per-batch variant filtration

2. PlotSVCountsPerSample: Visualize SV counts per 
   sample per type to help choose an IQR cutoff for 
   outlier filtering, and preview outlier samples for a given cutoff

3. FilterBatchSamples: Per-batch outlier sample filtration; 
   provide an appropriate outlier_cutoff_nIQR based on the 
   SV count plots and outlier previews from step 2. Note 
   that not removing high outliers can result in increased 
   compute cost and a higher false positive rate in later steps.

### Prerequisites

- Generate Batch Metrics

### Inputs

- Batch PED file
- Metrics file (GenerateBatchMetrics)
- Clustered SV and depth-only call VCFs (ClusterBatch)

### Outputs

- Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
- Filtered depth-only call VCF with outlier samples excluded
- Random forest cutoffs file
- PED file with outlier samples excluded