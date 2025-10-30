---
title: Joint calling
description: Run the pipeline on a cohort
sidebar_position: 4
slug: joint
---

## Terra workspace
Users should clone the [Terra joint calling workspace](https://app.terra.bio/#workspaces/broad-firecloud-dsde-methods/GATK-Structural-Variants-Joint-Calling)
which is configured with a demo sample set. 
Refer to the following sections for instructions on how to run the pipeline on your data using this workspace.

### Default data
The demonstration data in this workspace is 156 publicly-available 1000 Genomes Project samples from the 
[NYGC/AnVIL high coverage data set](https://app.terra.bio/#workspaces/anvil-datastorage/1000G-high-coverage-2019).

## Pipeline Expectations
### What does it do?
This pipeline performs structural variation discovery from CRAMs, joint genotyping, and variant resolution on a cohort 
of samples.

### Required inputs
Refer to the [Input Data section](/docs/gs/inputs) for details on file formats, sample QC, and sample ID restrictions.

The following inputs must be provided for each sample in the cohort, via the sample table described in **Workspace 
Setup** step 2:

|Input Type|Input Name|Description|
|---------|--------|--------------|
|`String`|`sample_id`|Case sample identifier|
|`File`|`bam_or_cram_file`|Path to the GCS location of the input CRAM or BAM file.|

The following cohort-level or batch-level inputs are also required:

|Input Type|Input Name|Description|
|---------|--------|--------------|
|`String`|`sample_set_id`|Batch identifier|
|`String`|`sample_set_set_id`|Cohort identifier|
|`File`|`cohort_ped_file`|Path to the GCS location of a family structure definitions file in [PED format](/docs/gs/inputs#ped-format).|

### Pipeline outputs

The following are the main pipeline outputs. For more information on the outputs of each module, refer to the 
[Modules section](/docs/category/modules).

|Output Type|Output Name|Description|
|---------|--------|--------------|
|`File`|`annotated_vcf`|Annotated SV VCF for the cohort***|
|`File`|`annotated_vcf_idx`|Index for `annotated_vcf`|
|`File`|`sv_vcf_qc_output`|QC plots (bundled in a .tar.gz file)|

***Note that this VCF is not filtered

### Pipeline overview

![Pipeline Diagram](https://media.githubusercontent.com/media/broadinstitute/gatk-sv/refs/tags/v1.0/terra_pipeline_diagram.jpg)

The following workflows and Jupyter notebooks are included in this workspace, to be executed in this order:

1. `01-GatherSampleEvidence`: Per-sample SV evidence collection, including calls from a configurable set of 
algorithms (Manta, MELT, and Wham), read depth (RD), split read positions (SR), and discordant pair positions (PE).
2. `02-EvidenceQC`: Dosage bias scoring and ploidy estimation, run on preliminary batches
3. [Notebook] `SampleQC.ipynb`: Interactively perform sample QC and filtering using outputs from `02-EvidenceQC`
4. [Notebook] `Batching.ipynb`: Create batches for subsequent steps. For cohorts >500 samples or smaller heterogeneous cohorts
5. `03-TrainGCNV`: Per-batch training of a gCNV model for use in `04-GatherBatchEvidence`
6. `04-GatherBatchEvidence`: Per-batch copy number variant calling using cn.MOPS and GATK gCNV; B-allele frequency (BAF) 
generation; call and evidence aggregation
7. `05-ClusterBatch`: Per-batch variant clustering
8. `06-GenerateBatchMetrics`: Per-batch metric generation
9. `07-FilterBatchSites`: Per-batch variant filtering and plot SV counts per sample per SV type to enable choice of IQR 
cutoff for outlier filtration in `08-FilterBatchSamples`
10. `08-FilterBatchSamples`: Per-batch outlier sample filtration
11. `09-MergeBatchSites`: Site merging of SVs discovered across batches, run on a cohort-level `sample_set_set`
12. `10-GenotypeBatch`: Per-batch genotyping of all sites in the cohort
13. `11-RegenotypeCNVs`: Cohort-level genotype refinement of some depth calls
14. `12-CombineBatches`: Cohort-level cross-batch integration and clustering
15. `13-ResolveComplexVariants`: Complex variant resolution
16. `14-GenotypeComplexVariants`: Complex variant re-genotyping
17. `15-CleanVcf`: VCF cleanup
18. `16-RefineComplexVariants`: Complex variant filtering and refinement
19. `17-JoinRawCalls`: Raw call aggregation
20. `18-SVConcordance`: Annotate genotype concordance with raw calls
21. `19-TrainGenotypeFilteringModel`: Training GQ recalibrator model for genotype filtering 
21. `19-FilterGenotypes`: Apply genotype filtering
22. `20-AnnotateVcf`: Cohort VCF annotations, including functional annotation, allele frequency (AF) annotation, and 
AF annotation with external population callsets

Extra workflows (Not part of canonical pipeline, but included for your convenience. May require manual configuration):
* `MainVcfQc`: Generate detailed call set QC plots
* `PlotSVCountsPerSample`: Plot SV counts per sample per SV type. Recommended to run before `FilterOutlierSamples` 
  (configured with the single VCF you want to filter) to enable IQR cutoff choice.
* `FilterOutlierSamples`: Filter outlier samples (in terms of SV counts) from a single VCF.
* `VisualizeCnvs`: Plot multi-sample depth profiles for CNVs

For detailed instructions on running the pipeline in Terra, see [workflow instructions](#instructions) below.

#### What is the maximum number of samples the pipeline can handle?

In Terra, we have tested batch sizes of up to 500 samples and cohort sizes (consisting of multiple batches) of up to 
11,000 samples (and 98,000 samples with the final steps split by chromosome). On a dedicated Cromwell server, we have 
tested the pipeline on cohorts of up to ~140,000 samples.


### Time and cost estimates

The following estimates pertain to the 1000 Genomes sample data in this workspace. They represent aggregated run time 
and cost across modules for the whole pipeline. For workflows run multiple times (on each sample or on each batch), 
the longest individual runtime was used. Call caching may affect some of this information.

|Number of samples|Time|Total run cost|Per-sample run cost|
|--------------|--------|----------|----------|
|312|~76 hours|~$675|~$2.16/sample|

Please note that sample characteristics, cohort size, and level of filtering may influence pipeline compute costs, 
with average costs ranging between $2-$3 per sample. For instance, PCR+ samples and samples with a high percentage 
of improperly paired reads have been observed to cost more. Consider 
[excluding low-quality samples](/docs/gs/inputs#sample-exclusion) prior to processing to keep costs low.


### Workspace setup

1. Clone this workspace into a Terra project to which you have access. Select `us-central1` for the workspace region.
   If you must use a different region, you will need to copy all GATK-SV docker images to the other region
   before running the pipeline. See the [docker images section](/docs/gs/dockers#regions-important) for details.

2. In your new workspace, delete the example data. To do this, go to the *Data* tab of the workspace. Delete the data 
   tables in this order: `sample_set_set`, `sample_set`, and `sample`. For each table, click the 3 dots icon to the 
   right of the table name and click "Delete table". Confirm when prompted.
   <img alt="deleting data tables" title="How to delete the sample data table" src="https://i.imgur.com/43M51WH.png" width="300" height="420" />

3. Create and upload a new sample data table for your samples. This should be a tab-separated file (.tsv) with one line 
   per sample, as well as a header (first) line. It should contain the columns `entity:sample_id` (first column) and 
   `bam_or_cram_file` at minimum. See the **Required inputs** section above for more information on these inputs. For 
   an example sample data table, refer to the sample data table for the 1000 Genomes samples in this workspace 
   [here in the GATK-SV GitHub repository](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/templates/terra_workspaces/cohort_mode/samples_1kgp_156.tsv.tmpl). 
   To upload the TSV file, navigate to the *Data* tab of the workspace, click the `Import Data` button on the top left, 
   and select "Upload TSV".
   <img alt="uploading a TSV data table" title="How to upload a TSV data table" src="https://i.imgur.com/1ZtwseH.png" width="300" height="250" />

4. Edit the `cohort_ped_file` item in the Workspace Data table (as shown in the screenshot below) to provide the Google 
   URI to the PED file for your cohort (make sure to share it with your Terra proxy account!).
   <img alt="editing cohort_ped_file" title="How to edit the cohort_ped_file attribute" src="https://i.imgur.com/IFwc0gs.png" width="800" height="400" />


### Creating sample_sets

To create batches (in the `sample_set` table), we recommend using the `Batching.ipynb` notebook (see [batching](#batching)). 
To create batches manually, the easiest way is to upload a tab-separated sample set membership file. 
This file should have one line per sample, plus a header (first) line. The first column should be 
`membership:sample_set_id` (containing the `sample_set_id` for the sample in question), and the second should be 
`sample` (containing the sample IDs). Recall that batch IDs (`sample_set_id`) should follow the 
[sample ID requirements](/docs/gs/inputs#sampleids). For an example sample set membership file, refer to
[this one in the GATK-SV GitHub repository](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/templates/terra_workspaces/cohort_mode/sample_set_membership_1kgp.tsv.tmpl).

## Workflow instructions {#instructions}

### General recommendations

* It is recommended to run each workflow first on one sample/batch to check that the method is properly configured 
before you attempt to process all of your data.
* We recommend enabling call-caching (on by default in each workflow configuration).
* We recommend enabling automatic intermediate file deletion by checking the box labeled "Delete intermediate outputs" 
at the top of the workflow launch page every time you start a workflow. With this option enabled, intermediate files 
(those not present in the Terra data table, and not needed for any further GATK-SV processing) will be deleted 
automatically if the workflow succeeds. If the workflow fails, the outputs will be retained to enable a re-run to 
pick up where it left off with call-caching. However, call-caching will not be possible for workflows that have 
succeeded. For more information on this option, see 
[this article](https://terra.bio/delete-intermediates-option-now-available-for-workflows-in-terra/). For guidance on 
managing intermediate storage from failed workflows, or from workflows without the delete intermediates option enabled, 
see the next bullet point.
* There are cases when you may need to manage storage in other ways: for workflows that failed (only delete files from 
a failed workflow after a version has succeeded, to avoid disabling call-caching), for workflows without intermediate 
file deletion enabled, or once you are done processing and want to delete files from earlier steps in the pipeline 
that you no longer need.
   * One option is to manually delete large files, or directories containing failed workflow intermediates (after 
  re-running the workflow successfully to take advantage of call-caching) with the command 
  `gsutil -m rm gs://path/to/workflow/directory/**file_extension_to_delete` to delete all files with the given extension 
  for that workflow, or `gsutil -m rm -r gs://path/to/workflow/directory/` to delete an entire workflow directory 
  (only after you are done with all the files!). Note that this can take a very long time for larger workflows, which 
  may contain thousands of files.
   * Another option is to use the `fiss mop` API call to delete all files that do not appear in one of the Terra data 
  tables (intermediate files). Always ensure that you are completely done with a step and you will not need to return 
  before using this option, as it will break call-caching. See 
  [this blog post](https://terra.bio/deleting-intermediate-workflow-outputs/) for more details. This can also be done 
  [via the command line](https://github.com/broadinstitute/fiss/wiki/MOP:-reducing-your-cloud-storage-footprint).
* If your workflow fails, check the job manager for the error message. Most issues can be resolved by increasing the 
memory or disk. Do not delete workflow log files until you are done troubleshooting. If call-caching is enabled, do not 
delete any files from the failed workflow until you have run it successfully.
* To display run costs, see [this article](https://support.terra.bio/hc/en-us/articles/360037862771#h_01EX5ED53HAZ59M29DRCG24CXY) 
for one-time setup instructions for non-Broad users.

### 01-GatherSampleEvidence

Read the full GatherSampleEvidence documentation [here](/docs/modules/gse).
* This workflow runs on a per-sample level, but you can launch many (a few hundred) samples at once, in arbitrary 
partitions. Make sure to try just one sample first though!
* Refer to the [Input Data section](/docs/gs/inputs) for details on file formats, sample QC, and sample ID restrictions.
* It is normal for a few samples in a cohort to run out of memory during Wham SV calling, so we recommend enabling 
auto-retry for out-of-memory errors for `01-GatherSampleEvidence` only. Before you launch the workflow, click the 
checkbox reading "Retry with more memory" and set the memory retry factor to 2. This action must be performed each 
time you launch a `01-GatherSampleEvidence` job.
* If you enable "Delete intermediate outputs" whenever you launch this workflow (recommended), BAM files will be 
deleted for successful runs; but BAM files will not be deleted if the run fails or if intermediate file deletion is 
not enabled. Since BAM files are large, we recommend deleting them to save on storage costs, but only after fixing and 
re-running the failed workflow, so that it will call-cache.


### 02-EvidenceQC

Read the full EvidenceQC documentation [here](/docs/modules/eqc).
* `02-EvidenceQC` is run on arbitrary cohort partitions of up to 500 samples.
* The outputs from `02-EvidenceQC` can be used for [sample QC](#sample-qc) and 
[batching](#batching) before moving on to [TrainGCNV](#traingcnv).


### Sample QC {#sample-qc}
Read the documentation on preliminary sample QC [here](/docs/modules/eqc#preliminary-sample-qc).
Follow the `SampleQC.ipynb` notebook step-by-step to evaluate sample data quality and remove low-quality samples as needed.
The notebook will produce a table of passing samples to use for [batching](#batching).


### Batching {#batching}
Read the documentation on batching [here](/docs/modules/eqc#batching).
If necessary, follow the `Batching.ipynb` notebook step-by-step to divide samples into batches
and create corresponding `sample_sets` for use in `03-TrainGCNV` and beyond.


### 03-TrainGCNV {#traingcnv}

Read the full TrainGCNV documentation [here](/docs/modules/gcnv).
* Before running this workflow, create the batches (~100-500 samples) you will use for the rest of the pipeline according 
to the [batching](#batching) instructions. These will likely not be the same as the batches you used for `02-EvidenceQC`.
* By default, `03-TrainGCNV` is configured to be run once per `sample_set` on 100 randomly-chosen samples from that 
set to create a gCNV model for each batch. To modify this behavior, you can set the `n_samples_subsample` parameter 
to the number of samples to use for training.

### 04-GatherBatchEvidence

Read the full GatherBatchEvidence documentation [here](/docs/modules/gbe).
* Use the same `sample_set` definitions you used for `03-TrainGCNV`.
* Before running this workflow, ensure that you have updated the `cohort_ped_file` attribute in Workspace Data with 
your cohort's PED file, with sex assignments updated based on ploidy detection from `02-EvidenceQC`.

### Steps 05-06

Read the full documentation for [ClusterBatch](/docs/modules/cb) and [GenerateBatchMetrics](/docs/modules/gbm).
* Use the same `sample_set` definitions you used for `03-TrainGCNV` and `04-GatherBatchEvidence`.


### Steps 07-08

These two workflows make up FilterBatch; they are subdivided in this workspace to enable tuning of outlier filtration 
cutoffs. Read the full FilterBatch documentation [here](/docs/modules/fb).
* Use the same `sample_set` definitions you used for `03-TrainGCNV` through `06-GenerateBatchMetrics`.
* `07-FilterBatchSites` produces SV count plots and files, as well as a preview of the outlier samples to be filtered. 
The input `N_IQR_cutoff_plotting` is used to visualize filtration thresholds on the SV count plots and preview the 
samples to be filtered; the default value is set to 6. You can adjust this value depending on your needs, and you can 
re-run the workflow with new `N_IQR_cutoff_plotting` values until the plots and outlier sample lists suit the purposes 
of your study. Once you have chosen an IQR cutoff, provide it to the `N_IQR_cutoff` input in `08-FilterBatchSamples` to 
filter the VCFs using the chosen cutoff.
* `08-FilterBatchSamples` performs outlier sample filtration, removing samples with an abnormal number of SV calls of 
at least one SV type. To tune the filtering threshold to your needs, edit the `N_IQR_cutoff` input value based on the 
plots and outlier sample preview lists from `07-FilterBatchSites`. The default value for `N_IQR_cutoff` in this step 
is 10000, which essentially means that no samples are filtered.

### 09-MergeBatchSites

Read the full MergeBatchSites documentation [here](/docs/modules/msites).
* `09-MergeBatchSites` is a cohort-level workflow, so it is run on a `sample_set_set` containing all the batches 
in the cohort. Navigate to the Data tab of your workspace. If there is no `sample_set_set` data table, you will need 
to create it. To do this, select the `sample_set` data table, then select (with the checkboxes) all the batches 
(`sample_set`) in your cohort. These should be the `sample_sets` that you used to run steps `03-TrainGCNV` through 
`08-FilterBatchSamples`. Then click the "Edit" icon above the table and choose "Save selection as set." Enter a name 
that follows the **Sample ID requirements**. This will create a new `sample_set_set` containing all of the `sample_sets` 
in your cohort. When you launch MergeBatchSites, you can now select this `sample_set_set`.

<img alt="selecting batches" title="Selecting sample_sets in the data table" src="https://i.imgur.com/E5x3qqk.png" width="400" height="200" />
<img alt="creating a new set" title="Creating a new sample_set_set" src="https://i.imgur.com/pizOtX9.png" width="400" height="200" />
* If there is already a `sample_set_set` data table in your workspace, you can create this `sample_set_set` while you 
are launching the `09-MergeBatchSites` workflow: click "Select Data", choose "Create new sample_set_set [...]", check 
all the batches to include (all the ones used in `03-TrainGCNV` through `08-FilterBatchSamples`), and give it a name 
that follows the [sample ID requirements](/docs/gs/inputs#sampleids).

<img alt="creating a cohort sample_set_set" title="How to create a cohort sample_set_set" src="https://i.imgur.com/zKEtSbe.png" width="377" height="363" />

### 10-GenotypeBatch

Read the full GenotypeBatch documentation [here](/docs/modules/gb).
* Use the same `sample_set` definitions you used for `03-TrainGCNV` through `08-FilterBatchSamples`.

### Steps 11-20

Read the full documentation for [RegenotypeCNVs](/docs/modules/rgcnvs), [CombineBatches](/docs/modules/cmb), 
[ResolveComplexVariants](/docs/modules/rcv), [GenotypeComplexVariants](/docs/modules/gcv), [CleanVcf](/docs/modules/cvcf), 
[RefineComplexVariants](/docs/modules/refcv), [JoinRawCalls](/docs/modules/jrc), [SVConcordance](/docs/modules/svc), 
[FilterGenotypes](/docs/modules/fg), and [AnnotateVcf](/docs/modules/av).
* Use the same cohort `sample_set_set` you created and used for `09-MergeBatchSites`.
