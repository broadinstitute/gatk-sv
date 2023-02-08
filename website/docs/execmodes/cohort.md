---
sidebar_position: 2
---

# Cohort mode

A minimum cohort size of 100 is required, and a roughly equal 
number of males and females is recommended. For modest cohorts 
(~100-500 samples), the pipeline can be run as a single batch 
using `GATKSVPipelineBatch.wdl`.

For larger cohorts, samples should be split up into batches of 
about 100-500 samples. 

The pipeline should be executed as follows:

- Modules [GatherSampleEvidence](/docs/modules/gse) and 
  [EvidenceQC](/docs/modules/eqc) can be run on arbitrary cohort partitions

- Modules [GatherBatchEvidence](/docs/modules/gbe), 
  [ClusterBatch](/docs/modules/cb), [GenerateBatchMetrics](/docs/modules/gbm), 
  and [FilterBatch](/docs/modules/fb) are run separately per batch

- [GenotypeBatch](/docs/modules/gb) is run separately per batch, 
  using filtered variants ([FilterBatch](/docs/modules/fb) output) combined across all batches

- [MakeCohortVcf](/docs/modules/cvcf) and beyond are run on all batches together

Note: [GatherBatchEvidence](/docs/modules/gbe) requires a [trained gCNV model](/docs/modules/gcnv).
