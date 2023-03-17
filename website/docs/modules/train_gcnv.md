---
title: TrainGCNV
description: Train gCNV
sidebar_position: 3
slug: gcnv
---

Trains a gCNV model for use in GatherBatchEvidence. 
The WDL can be found at /wdl/TrainGCNV.wdl. 

Both the cohort and single-sample modes use the 
GATK gCNV depth calling pipeline, which requires a 
trained model as input. The samples used for training 
should be technically homogeneous and similar to the 
samples to be processed (i.e. same sample type, 
library prep protocol, sequencer, sequencing center, etc.). 
The samples to be processed may comprise all or a subset 
of the training set. For small, relatively homogenous cohorts, 
a single gCNV model is usually sufficient. If a cohort 
contains multiple data sources, we recommend training a separate 
model for each batch or group of batches with similar dosage 
score (WGD). The model may be trained on all or a subset of 
the samples to which it will be applied; a reasonable default 
is 100 randomly-selected samples from the batch (the random 
selection can be done as part of the workflow by specifying 
a number of samples to the n_samples_subsample input 
parameter in /wdl/TrainGCNV.wdl).

### Prerequisites

- GatherSampleEvidence
- (Recommended) EvidenceQC

### Inputs

- Read count files (GatherSampleEvidence)

### Outputs

- Contig ploidy model tarball
- gCNV model tarballs