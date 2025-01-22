---
title: Input data
description: Supported input and reference data.
sidebar_position: 2
slug: ./inputs
---

GATK-SV requires the following input data:

1. Sequencing alignments in BAM or CRAM format that are:
   - Short-read, paired-end Illumina (e.g. Novaseq)
   - Deep whole-genome coverage (~30x); RNA-seq and targeted (exome) libraries are not supported
   - Indexed (have a companion `.bai` or `.crai` file)
   - Aligned to 
     [hg38](https://gatk.broadinstitute.org/hc/en-us/articles/360035890951-Human-genome-reference-builds-GRCh38-or-hg38-b37-hg19) 
     with either [GATK Best Practices](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows) 
     and [bwa-mem](https://github.com/lh3/bwa), 
     or [Illumina DRAGEN](https://www.illumina.com/products/by-type/informatics-products/dragen-secondary-analysis.html) v3.4.12 or v3.7.8
2. (Joint calling mode only) Family structure definitions file in [PED format](/docs/gs/inputs#ped-format). This file is required even if your dataset does not contain related individuals.

Note that the supported alignment pipeline versions have been extensively tested for robustness and accuracy. While other 
versions of DRAGEN may work as well, they have not been validated with GATK-SV. We do not recommend mixing aligners within call sets. 

### Sample Exclusion
We recommend filtering out samples with a high percentage 
of improperly paired or chimeric reads 
as technical outliers prior to running [GatherSampleEvidence](/docs/modules/gse). 
Samples with high rates of anomalous reads may indicate issues 
with library preparation, degradation, or contamination and can lead to poor variant set quality. 
Samples failing these criteria often require longer run times and higher compute costs.


### Sample IDs {#sampleids}

GATK-SV imposes certain restrictions on sample names (IDs) in order to avoid certain parsing errors (e.g. with the 
use of the `grep` command). While future releases will obviate some of these restrictions, users must modify 
their sample IDs according to the following requirements.

#### Sample IDs must:
- Be unique within the cohort
- Contain only alphanumeric characters and underscores (no dashes, whitespace, or special characters)

#### Sample IDs should not:
- Contain only numeric characters, e.g. `10004928`
- Be a substring of another sample ID in the same cohort
- Contain any of the following substrings: `chr`, `name`, `DEL`, `DUP`, `CPX`, `CHROM`

The same requirements apply to family IDs in the PED file, as well as batch IDs and the cohort ID provided as workflow 
inputs. [This script](https://github.com/broadinstitute/gatk-sv/blob/main/scripts/inputs/convert_sample_ids.py)
can be used to transform a list of sample IDs to meet safe ID requirements.

Users should assign sample IDs in [GatherSampleEvidence](/docs/modules/gse) with the `sample_id` input, which needs not 
match the sample name defined in the BAM/CRAM header. Alternatively, sample IDs can be replaced again in 
[GatherBatchEvidence](/docs/modules/gbe) by setting the parameter `rename_samples = True` and providing updated 
sample IDs via the `samples` parameter.

Note that following inputs will need to be updated with the transformed sample IDs:

- Sample ID list for [GatherSampleEvidence](/docs/modules/gse) or [GatherBatchEvidence](/docs/modules/gbe)
- PED file


### PED file format {#ped-format}
The PED file format is described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format). Note that GATK-SV imposes additional requirements:
* The file must be tab-delimited.
* The sex column must only contain 0, 1, or 2: 1=Male, 2=Female, 0=Other/Unknown. Sex chromosome aneuploidies (detected in [EvidenceQC](/docs/modules/eqc)) should be entered as sex = 0.
* All family, individual, and parental IDs must conform to the [sample ID requirements](/docs/gs/inputs#sampleids).
* Missing parental IDs should be entered as 0.
* Header lines are allowed if they begin with a # character.
* To validate the PED file, you may use `src/sv-pipeline/scripts/validate_ped.py -p pedigree.ped -s samples.list`.
