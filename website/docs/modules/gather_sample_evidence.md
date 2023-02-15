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

### Sample ID requirements
#### Sample IDs must

- Be unique within the cohort
- Contain only alphanumeric characters and underscores (no dashes, whitespace, or special characters)

Sample IDs should not:

- Contain only numeric characters
- Be a substring of another sample ID in the same cohort
- Contain any of the following substrings: `chr`, `name`, `DEL`, `DUP`, `CPX`, `CHROM`

The same requirements apply to family IDs in the PED file, 
as well as batch IDs and the cohort ID provided as workflow inputs.

Sample IDs are provided to GatherSampleEvidence directly and 
need not match sample names from the BAM/CRAM headers. 
`GetSampleID.wdl` can be used to fetch BAM sample IDs and 
also generates a set of alternate IDs that are considered 
safe for this pipeline; alternatively, [this script](https://github.com/talkowski-lab/gnomad_sv_v3/blob/master/sample_id/convert_sample_ids.py)
transforms a list of sample IDs to fit these requirements. 
Currently, sample IDs can be replaced again in [GatherBatchEvidence](./gbe).

The following inputs will need to be updated with the transformed sample IDs:

- Sample ID list for [GatherSampleEvidence](./gse) or [GatherBatchEvidence](./gbe)
- PED file

### Inputs

- Per-sample BAM or CRAM files aligned to hg38. Index files (.bai) must be provided if using BAMs.

### Outputs

- Caller VCFs (Manta, MELT, and/or Wham)
- Binned read counts file
- Split reads (SR) file
- Discordant read pairs (PE) file
