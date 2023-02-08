---
title: GenerateBatchMetrics
description: Generate Batch Metrics
sidebar_position: 6
slug: gbm
---

Generates variant metrics for filtering.

### Prerequisites

- Cluster batch

### Inputs

- Combined read count matrix, SR, PE, and BAF files (GatherBatchEvidence)
- Per-sample median coverage estimates (GatherBatchEvidence)
- Clustered SV VCFs (ClusterBatch)
- Clustered depth-only call VCF (ClusterBatch)

### Outputs

- Metrics file
