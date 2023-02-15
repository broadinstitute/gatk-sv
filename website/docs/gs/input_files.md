---
title: Input Data
description: Supported input and reference data.
sidebar_position: 3
slug: ./inputs
---

GATK-SV requires the following input data:

- Illumina short-read whole-genome CRAMs or BAMs, aligned to hg38 with [bwa-mem](https://github.com/lh3/bwa). 
  BAMs must also be indexed.

- Family structure definitions file in 
  [PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format).  
  **Note**: Sex aneuploidies (detected in EvidenceQC) should be entered as sex = 0.

### Sample Exclusion
We recommend filtering out samples with a high percentage 
of improperly paired reads (>10% or an outlier for your data) 
as technical outliers prior to running [GatherSampleEvidence](/docs/modules/gse). 
A high percentage of improperly paired reads may indicate issues 
with library prep, degradation, or contamination. Artifactual 
improperly paired reads could cause incorrect SV calls, and 
these samples have been observed to have longer runtimes and 
higher compute costs for [GatherSampleEvidence](/docs/modules/gse).


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
