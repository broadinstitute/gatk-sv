---
title: GatherBatchEvidence
description: Gather Batch Evidence
sidebar_position: 4
slug: gbe
---

Runs CNV callers (cnMOPs, GATK gCNV) and combines single-sample 
raw evidence into a batch. See above for more information on batching.

### Prerequisites

- GatherSampleEvidence
- (Recommended) EvidenceQC
- gCNV training. 

### Inputs
- PED file (updated with EvidenceQC sex assignments, including sex = 0 
  for sex aneuploidies. Calls will not be made on sex chromosomes 
  when sex = 0 in order to avoid generating many confusing calls 
  or upsetting normalized copy numbers for the batch.)
- Read count, BAF, PE, SD, and SR files (GatherSampleEvidence)
- Caller VCFs (GatherSampleEvidence)
- Contig ploidy model and gCNV model files (gCNV training)

### Outputs

- Combined read count matrix, SR, PE, and BAF files
- Standardized call VCFs
- Depth-only (DEL/DUP) calls
- Per-sample median coverage estimates
- (Optional) Evidence QC plots
