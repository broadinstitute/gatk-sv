---
sidebar_position: 2
---

# Cohort mode

A minimum cohort size of 100 is required, and a roughly equal 
number of males and females is recommended. For modest cohorts 
(~100-500 samples), the pipeline can be run as a single batch 
using `GATKSVPipelineBatch.wdl`.

For larger cohorts, samples should be split up into batches of 
about 100-500 samples. Refer to the [Batching](/docs/run/batching) 
section for further guidance on creating batches.

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

#### <a name="batching">Batching</a>
For larger cohorts, samples should be split up into batches of about 100-500 
samples with similar characteristics. We recommend batching based on overall 
coverage and dosage score (WGD), which can be generated in [EvidenceQC](/docs/modules/eqc). 
An example batching process is outlined below:

1. Divide the cohort into PCR+ and PCR- samples
2. Partition the samples by median coverage from [EvidenceQC](/docs/modules/eqc), 
   grouping samples with similar median coverage together. The end goal is to 
   divide the cohort into roughly equal-sized batches of about 100-500 samples; 
   if your partitions based on coverage are larger or uneven, you can partition 
   the cohort further in the next step to obtain the final batches. 
3. Optionally, divide the samples further by dosage score (WGD) from 
   [EvidenceQC](/docs/modules/eqc), grouping samples with similar WGD score 
   together, to obtain roughly equal-sized batches of about 100-500 samples
4. Maintain a roughly equal sex balance within each batch, based on sex 
   assignments from [EvidenceQC](/docs/modules/eqc)