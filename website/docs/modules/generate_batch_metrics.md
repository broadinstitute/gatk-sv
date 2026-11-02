---
title: GenerateBatchMetrics
description: Generate Batch Metrics
sidebar_position: 6
slug: gbm
---

import { Highlight, HighlightOptionalArg } from "@site/src/components/highlight.js"

[WDL source code](https://github.com/broadinstitute/gatk-sv/blob/main/wdl/GenerateBatchMetrics.wdl)

Analyzes variants for RD, BAF, PE, and SR evidence and creates a table of metrics containing raw and statistical 
metrics. These results are used to assess variant quality in `FilterBatch` and for SR-based breakpoint refinement.

Modified tests are applied to common variants (carrier frequency at least 50%) and results are emitted in a separate table.

The following diagram illustrates the recommended invocation order:

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#caf0f8,color:#00509d
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  cb: ClusterBatch
  gbm: GenerateBatchMetrics
  fb: FilterBatch
  cb --> gbm
  gbm --> fb
  
  class gbm thisModule
  class cb inModules
  class fb outModules
```

## Inputs

#### `batch`
An identifier for the batch. Should match the name used in [GatherBatchEvidence](./gbe#batch).

#### `*_vcf`
Clustered VCFs from [ClusterBatch](./cb#clustered__vcf).

#### `baf_file`
Merged BAF evidence file from [GatherBatchEvidence](./gbe#merged_baf).

#### `pe_file`
Merged PE evidence file from [GatherBatchEvidence](./gbe#merged_pe).

#### `rd_file`
Merged RD evidence file from [GatherBatchEvidence](./gbe#merged_bincov).

#### `sr_file`
Merged SR evidence file from [GatherBatchEvidence](./gbe#merged_sr).

#### `median_file`
Merged median coverage table from [GatherBatchEvidence](./gbe#median_cov).

#### `ped_file`
Family structures and sex assignments determined in [EvidenceQC](./eqc). See [PED file format](/docs/gs/inputs#ped-format).

#### `records_per_shard_agg`
Variants per shard for evidence testing. Reduce defaults to increase parallelism if the workflow is too slow.

#### <HighlightOptionalArg>Optional</HighlightOptionalArg> `outlier_sample_ids`
Provide a file containing sample IDs, delimited by new lines, to exclude from the generation of batch metrics - 
this should be based on outlier samples identified following [EvidenceQC](./eqc) that are still retained in 
the cohort. If provided, the workflow discards any outlier samples from being considered as part of the called 
samples when calculating metrics at a given site, as long as non-outlier samples are also called for that site. 
It does the same for the set of background samples considered in the metric calculations for a given site. These 
outlier samples are *not* removed from joint calling in downstream modules.

## Outputs

#### `metrics`
TSV of variant metrics

#### `ploidy_table`
TSV of contig ploidies for all samples, assuming diploid autosome and sex assignments from the [ped_file](#ped_file).