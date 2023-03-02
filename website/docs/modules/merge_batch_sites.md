---
title: MergeBatchSites
description: Merge Batch Sites
sidebar_position: 8
slug: msites
---

Combines filtered variants across batches. The WDL can be found at: /wdl/MergeBatchSites.wdl.

### Prerequisites
- Filter Batch

### Inputs

- List of filtered PESR VCFs (FilterBatch)
- List of filtered depth VCFs (FilterBatch)

### Outputs

- Combined cohort PESR and depth VCFs
