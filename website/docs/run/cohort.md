---
title: Cohort
description: Run the pipeline on a cohort
sidebar_position: 4
slug: cohort
---

## Batching

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