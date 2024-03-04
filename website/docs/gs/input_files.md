---
title: Input Data
description: Supported input and reference data.
sidebar_position: 5
slug: ./inputs
---

GATK-SV requires the following input data:

- Illumina short-read whole-genome CRAMs or BAMs, aligned to hg38 with [bwa-mem](https://github.com/lh3/bwa). 
  BAMs must also be indexed.

- Family structure definitions file in 
  [PED format](/docs/gs/inputs#ped-format).

### PED file format {#ped-format}
The PED file format is described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format). Note that GATK-SV imposes additional requirements:
* The file must be tab-delimited.
* The sex column must only contain 0, 1, or 2: 1=Male, 2=Female, 0=Other/Unknown. Sex chromosome aneuploidies (detected in [EvidenceQC](/docs/modules/eqc)) should be entered as sex = 0.
* All family, individual, and parental IDs must conform to the [sample ID requirements](/docs/gs/inputs#sampleids).
* Missing parental IDs should be entered as 0.
* Header lines are allowed if they begin with a # character.
To validate the PED file, you may use `src/sv-pipeline/scripts/validate_ped.py -p pedigree.ped -s samples.list`.

### Sample Exclusion
We recommend filtering out samples with a high percentage 
of improperly paired reads (>10% or an outlier for your data) 
as technical outliers prior to running [GatherSampleEvidence](/docs/modules/gse). 
A high percentage of improperly paired reads may indicate issues 
with library prep, degradation, or contamination. Artifactual 
improperly paired reads could cause incorrect SV calls, and 
these samples have been observed to have longer runtimes and 
higher compute costs for [GatherSampleEvidence](/docs/modules/gse).


### Sample ID requirements {#sampleids}
#### Sample IDs must

- Be unique within the cohort
- Contain only alphanumeric characters and underscores (no dashes, whitespace, or special characters)

#### Sample IDs should not

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
Currently, sample IDs can be replaced again in [GatherBatchEvidence](/docs/modules/gbe).

The following inputs will need to be updated with the transformed sample IDs:

- Sample ID list for [GatherSampleEvidence](/docs/modules/gse) or [GatherBatchEvidence](/docs/modules/gbe)
- PED file


### Generating a reference panel (single-sample mode only)
New reference panels can be generated easily from a single run of the `GATKSVPipelineBatch` workflow. If using a Cromwell server, we recommend copying the outputs to a permanent location by adding the following option to the workflow configuration file:
```
  "final_workflow_outputs_dir" : "gs://my-outputs-bucket",
  "use_relative_output_paths": false,
```
Here is an example of how to generate workflow input jsons from `GATKSVPipelineBatch` workflow metadata:
```
> cromshell -t60 metadata 38c65ca4-2a07-4805-86b6-214696075fef > metadata.json
> python scripts/inputs/create_test_batch.py \
    --execution-bucket gs://my-exec-bucket \
    --final-workflow-outputs-dir gs://my-outputs-bucket \
    metadata.json \
    > inputs/values/my_ref_panel.json
> # Define your google project id (for Cromwell inputs) and Terra billing project (for workspace inputs)
> echo '{ "google_project_id": "my-google-project-id", "terra_billing_project_id": "my-terra-billing-project" }' > inputs/values/google_cloud.my_project.json
> # Build test files for batched workflows (google cloud project id required)
> python scripts/inputs/build_inputs.py \
    inputs/values \
    inputs/templates/test \
    inputs/build/my_ref_panel/test \
    -a '{ "test_batch" : "ref_panel_1kg", "cloud_env": "google_cloud.my_project" }'
> # Build test files for the single-sample workflow
> python scripts/inputs/build_inputs.py \
    inputs/values \
    inputs/templates/test/GATKSVPipelineSingleSample \
    inputs/build/NA19240/test_my_ref_panel \
    -a '{ "single_sample" : "test_single_sample_NA19240", "ref_panel" : "my_ref_panel" }'
> # Build files for a Terra workspace
> python scripts/inputs/build_inputs.py \
    inputs/values \
    inputs/templates/terra_workspaces/single_sample \
    inputs/build/NA12878/terra_my_ref_panel \
    -a '{ "single_sample" : "test_single_sample_NA12878", "ref_panel" : "my_ref_panel" }'
```
Note that the inputs to `GATKSVPipelineBatch` may be used as resources for the reference panel and therefore should also be in a permanent location.
