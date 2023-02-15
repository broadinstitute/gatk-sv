---
title: GatherSampleEvidence 
description: Gather Sample Evidence
sidebar_position: 1
slug: gse
---

Runs raw evidence collection on each sample with the following SV callers: 
Manta, Wham, and/or MELT. For guidance on pre-filtering prior to GatherSampleEvidence, 
refer to the Sample Exclusion section.

Note: a list of sample IDs must be provided. Refer to the sample ID 
requirements for specifications of allowable sample IDs. 
IDs that do not meet these requirements may cause errors.

### Inputs

- Per-sample BAM or CRAM files aligned to hg38. Index files (.bai) must be provided if using BAMs.

### Outputs

- Caller VCFs (Manta, MELT, and/or Wham)
- Binned read counts file
- Split reads (SR) file
- Discordant read pairs (PE) file
