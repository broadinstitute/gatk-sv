# GATK-SV

A structural variation discovery pipeline for Illumina short-read whole-genome sequencing (WGS) data.

## Table of Contents
* [Requirements](#requirements)
* [Citation](#citation)
* [Quickstart](#quickstart)
* [Pipeline Overview](#overview)
    * [Cohort mode](#cohort-mode)
    * [Single-sample mode](#single-sample-mode)
    * [gCNV model](#gcnv-training-overview)
    * [Generating a reference panel](#reference-panel-generation)
* [Module Descriptions](#descriptions)
    * [GatherSampleEvidence](#gather-sample-evidence) - Raw callers and evidence collection
    * [EvidenceQC](#evidence-qc) - Batch QC
    * [TrainGCNV](#gcnv-training) - gCNV model creation
    * [GatherBatchEvidence](#gather-batch-evidence) - Batch evidence merging, BAF generation, and depth callers
    * [ClusterBatch](#cluster-batch) - Site clustering
    * [GenerateBatchMetrics](#generate-batch-metrics) - Site metrics
    * [FilterBatch](#filter-batch) - Filtering
    * [MergeBatchSites](#merge-batch-sites) - Cross-batch site merging
    * [GenotypeBatch](#genotype-batch) - Genotyping
    * [RegenotypeCNVs](#regenotype-cnvs) - Genotype refinement (optional)
    * [MakeCohortVcf](#make-cohort-vcf) - Cross-batch integration, complex event resolution, and VCF cleanup
    * [Module 07](#module07) - Downstream Filtering
    * [AnnotateVcf](#annotate-vcf) - Annotation
    * [Module 09](#module09) - QC and Visualization
    * Additional modules - Mosaic and de novo
* [CI/CD](#cicd)
* [Troubleshooting](#troubleshooting)


## <a name="requirements">Requirements</a>


### Deployment and execution:
* A [Google Cloud](https://cloud.google.com/) account.
* A workflow execution system supporting the [Workflow Description Language](https://openwdl.org/) (WDL), either:
  * [Cromwell](https://github.com/broadinstitute/cromwell) (v36 or higher). A dedicated server is highly recommended.
  * or [Terra](https://terra.bio/) (note preconfigured GATK-SV workflows are not yet available for this platform)
* Recommended: [MELT](https://melt.igs.umaryland.edu/). Due to licensing restrictions, we cannot provide a public docker image or reference panel VCFs for this algorithm.
* Recommended: [cromshell](https://github.com/broadinstitute/cromshell) for interacting with a dedicated Cromwell server.
* Recommended: [WOMtool](https://cromwell.readthedocs.io/en/stable/WOMtool/) for validating WDL/json files.

#### Alternative backends
Because GATK-SV has been tested only on the Google Cloud Platform (GCP), we are unable to provide specific guidance or support for other execution platforms including HPC clusters and AWS. Contributions from the community to improve portability between backends will be considered on a case-by-case-basis. We ask contributors to please adhere to the following guidelines when submitting issues and pull requests:

1. Code changes must be functionally equivalent on GCP backends, i.e. not result in changed output
2. Increases to cost and runtime on GCP backends should be minimal
3. Avoid adding new inputs and tasks to workflows. Simpler changes are more likely to be approved, e.g. small in-line changes to scripts or WDL task command sections
4. Avoid introducing new code paths, e.g. conditional statements
5. Additional backend-specific scripts, workflows, tests, and Dockerfiles will not be approved
6. Changes to Dockerfiles may require extensive testing before approval

We still encourage members of the community to adapt GATK-SV for non-GCP backends and share code on forked repositories. Here are a some considerations:
* Refer to Cromwell's [documentation](https://cromwell.readthedocs.io/en/stable/backends/Backends/) for configuration instructions.
* The handling and ordering of `glob` commands may differ between platforms.
* Shell commands that are potentially destructive to input files (e.g. `rm`, `mv`, `tabix`) can cause unexpected behavior on shared filesystems. Enabling [copy localization](https://cromwell.readthedocs.io/en/stable/Configuring/#local-filesystem-options) may help to more closely replicate the behavior on GCP.
* For clusters that do not support Docker, Singularity is an alternative. See [Cromwell documentation on Singularity(https://cromwell.readthedocs.io/en/stable/tutorials/Containers/#singularity).
* The GATK-SV pipeline takes advantage of the massive parallelization possible in the cloud. Local backends may not have the resources to execute all of the workflows. Workflows that use fewer resources or that are less parallelized may be more successful. For instance, some users have been able to run [GatherSampleEvidence](#gather-sample-evidence) on a SLURM cluster.

### Data:
* Illumina short-read whole-genome CRAMs or BAMs, aligned to hg38 with [bwa-mem](https://github.com/lh3/bwa). BAMs must also be indexed.
* Indexed GVCFs produced by GATK HaplotypeCaller, or a jointly genotyped VCF.
* Family structure definitions file in [PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format). Sex aneuploidies (detected in [EvidenceQC](#evidence-qc)) should be entered as sex = 0.

#### <a name="sample-exclusion">Sample Exclusion</a>
We recommend filtering out samples with a high percentage of improperly paired reads (>10% or an outlier for your data) as technical outliers prior to running [GatherSampleEvidence](#gather-sample-evidence). A high percentage of improperly paired reads may indicate issues with library prep, degradation, or contamination. Artifactual improperly paired reads could cause incorrect SV calls, and these samples have been observed to have longer runtimes and higher compute costs for [GatherSampleEvidence](#gather-sample-evidence).

#### <a name="sampleids">Sample ID requirements:</a>

Sample IDs must:
* Be unique within the cohort
* Contain only alphanumeric characters and underscores (no dashes, whitespace, or special characters)

Sample IDs should not:
* Contain only numeric characters
* Be a substring of another sample ID in the same cohort
* Contain any of the following substrings: `chr`, `name`, `DEL`, `DUP`, `CPX`, `CHROM`

The same requirements apply to family IDs in the PED file, as well as batch IDs and the cohort ID provided as workflow inputs.

Sample IDs are provided to [GatherSampleEvidence](#gather-sample-evidence) directly and need not match sample names from the BAM/CRAM headers or GVCFs. `GetSampleID.wdl` can be used to fetch BAM sample IDs and also generates a set of alternate IDs that are considered safe for this pipeline; alternatively, [this script](https://github.com/talkowski-lab/gnomad_sv_v3/blob/master/sample_id/convert_sample_ids.py) transforms a list of sample IDs to fit these requirements. Currently, sample IDs can be replaced again in [GatherBatchEvidence](#gather-batch-evidence). 

The following inputs will need to be updated with the transformed sample IDs:
* Sample ID list for [GatherSampleEvidence](#gather-sample-evidence) or [GatherBatchEvidence](#gather-batch-evidence)
* PED file

If using a SNP VCF in [GatherBatchEvidence](#gather-batch-evidence), it does not need to be re-headered; simply provide the `vcf_samples` argument.


## <a name="citation">Citation</a>
Please cite the following publication:
[Collins, Brand, et al. 2020. "A structural variation reference for medical and population genetics." Nature 581, 444-451.](https://doi.org/10.1038/s41586-020-2287-8)

Additional references:
[Werling et al. 2018. "An analytical framework for whole-genome sequence association studies and its implications for autism spectrum disorder." Nature genetics 50.5, 727-736.](http://dx.doi.org/10.1038/s41588-018-0107-y)


## <a name="quickstart">Quickstart</a>

#### WDLs
There are two scripts for running the full pipeline:
* `wdl/GATKSVPipelineBatch.wdl`: Runs GATK-SV on a batch of samples.
* `wdl/GATKSVPipelineSingleSample.wdl`: Runs GATK-SV on a single sample, given a reference panel

#### Building inputs
Example workflow inputs can be found in `/inputs`. Build using `scripts/inputs/build_default_inputs.sh`, which 
generates input jsons in `/inputs/build`. Except the MELT docker image, all required resources are available in public 
Google buckets. 

Some workflows require a Google Cloud Project ID to be defined in a cloud environment parameter group. Workspace builds 
require a Terra billing project ID as well. An example is  provided at `/inputs/values/google_cloud.json` but should 
not be used, as modifying this file will cause tracked changes in the repository. Instead, create a copy in the same 
directory with the format `google_cloud.my_project.json` and modify as necessary.

Note that these inputs are required only when certain data are located in requester pays buckets. If this does not 
apply, users may use placeholder values for the cloud configuration and simply delete the inputs manually.

#### MELT
**Important**: The example input files contain MELT inputs that are NOT public (see [Requirements](#requirements)). These include:

* `GATKSVPipelineSingleSample.melt_docker` and `GATKSVPipelineBatch.melt_docker` - MELT docker URI (see [Docker readme](https://github.com/talkowski-lab/gatk-sv-v1/blob/master/dockerfiles/README.md))
* `GATKSVPipelineSingleSample.ref_std_melt_vcfs` - Standardized MELT VCFs ([GatherBatchEvidence](#gather-batch-evidence))

The input values are provided only as an example and are not publicly accessible. In order to include MELT, these values must be provided by the user. MELT can be disabled by deleting these inputs and setting `GATKSVPipelineBatch.use_melt` to `false`.

#### Requester pays buckets
**Important**: The following parameters must be set when certain input data is in requester pays (RP) buckets:

* `GATKSVPipelineSingleSample.requester_pays_cram` and `GATKSVPipelineBatch.GatherSampleEvidenceBatch.requester_pays_crams` - set to `True` if inputs are CRAM format and in an RP bucket, otherwise `False`.
* `GATKSVPipelineBatch.GATKSVPipelinePhase1.gcs_project_for_requester_pays` - set to your Google Cloud Project ID if gVCFs are in an RP bucket, otherwise omit this parameter.

#### Execution
We recommend running the pipeline on a dedicated [Cromwell](https://github.com/broadinstitute/cromwell) server with a [cromshell](https://github.com/broadinstitute/cromshell) client. A batch run can be started with the following commands:

```
> mkdir gatksv_run && cd gatksv_run
> mkdir wdl && cd wdl
> cp $GATK_SV_ROOT/wdl/*.wdl .
> zip dep.zip *.wdl
> cd ..
> echo '{ "google_project_id": "my-google-project-id", "terra_billing_project_id": "my-terra-billing-project" }' > inputs/values/google_cloud.my_project.json
> bash scripts/inputs/build_default_inputs.sh -d $GATK_SV_ROOT -c google_cloud.my_project
> cp $GATK_SV_ROOT/inputs/build/ref_panel_1kg/test/GATKSVPipelineBatch/GATKSVPipelineBatch.json GATKSVPipelineBatch.my_run.json
> cromshell submit wdl/GATKSVPipelineBatch.wdl GATKSVPipelineBatch.my_run.json cromwell_config.json wdl/dep.zip
```

where `cromwell_config.json` is a Cromwell [workflow options file](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/). Note users will need to re-populate batch/sample-specific parameters (e.g. BAMs and sample IDs).

## <a name="overview">Pipeline Overview</a>
The pipeline consists of a series of modules that perform the following:
* [GatherSampleEvidence](#gather-sample-evidence): SV evidence collection, including calls from a configurable set of algorithms (Delly, Manta, MELT, and Wham), read depth (RD), split read positions (SR), and discordant pair positions (PE).
* [EvidenceQC](#evidence-qc): Dosage bias scoring and ploidy estimation
* [GatherBatchEvidence](#gather-batch-evidence): Copy number variant calling using cn.MOPS and GATK gCNV; B-allele frequency (BAF) generation; call and evidence aggregation
* [ClusterBatch](#cluster-batch): Variant clustering
* [GenerateBatchMetrics](#generate-batch-metrics): Variant filtering metric generation
* [FilterBatch](#filter-batch): Variant filtering; outlier exclusion
* [GenotypeBatch](#genotype-batch): Genotyping
* [MakeCohortVcf](#make-cohort-vcf): Cross-batch integration; complex variant resolution and re-genotyping; vcf cleanup
* [Module 07](#module07): Downstream filtering, including minGQ, batch effect check, outlier samples removal and final recalibration;
* [AnnotateVcf](#annotate-vcf): Annotations, including functional annotation, allele frequency (AF) annotation and AF annotation with external population callsets;
* [Module 09](#module09): Visualization, including scripts that generates IGV screenshots and rd plots.
* Additional modules to be added: de novo and mosaic scripts 


Repository structure:
* `/dockerfiles`: Resources for building pipeline docker images
* `/inputs`: files for generating workflow inputs
  * `/templates`: Input json file templates
  * `/values`: Input values used to populate templates
* `/wdl`: WDLs running the pipeline. There is a master WDL for running each module, e.g. `ClusterBatch.wdl`.
* `/scripts`: scripts for running tests, building dockers, and analyzing cromwell metadata files
* `/src`: main pipeline scripts
  * `/RdTest`: scripts for depth testing
  * `/sv-pipeline`: various scripts and packages used throughout the pipeline
  * `/svqc`: Python module for checking that pipeline metrics fall within acceptable limits
  * `/svtest`: Python module for generating various summary metrics from module outputs
  * `/svtk`: Python module of tools for SV-related datafile parsing and analysis
  * `/WGD`: whole-genome dosage scoring scripts


## <a name="cohort-mode">Cohort mode</a>
A minimum cohort size of 100 with roughly equal number of males and females is recommended. For modest cohorts (~100-500 samples), the pipeline can be run as a single batch using `GATKSVPipelineBatch.wdl`.

For larger cohorts, samples should be split up into batches of about 100-500 samples. Refer to the [Batching](#batching) section for further guidance on creating batches.

The pipeline should be executed as follows:
* Modules [GatherSampleEvidence](#gather-sample-evidence) and [EvidenceQC](#evidence-qc) can be run on arbitrary cohort partitions
* Modules [GatherBatchEvidence](#gather-batch-evidence), [ClusterBatch](#cluster-batch), [GenerateBatchMetrics](#generate-batch-metrics), and [FilterBatch](#filter-batch) are run separately per batch
* [GenotypeBatch](#genotype-batch) is run separately per batch, using filtered variants ([FilterBatch](#filter-batch) output) combined across all batches
* [MakeCohortVcf](#make-cohort-vcf) and beyond are run on all batches together

Note: [GatherBatchEvidence](#gather-batch-evidence) requires a [trained gCNV model](#gcnv-training).

#### <a name="batching">Batching</a>
For larger cohorts, samples should be split up into batches of about 100-500 samples with similar characteristics. We recommend batching based on overall coverage and dosage score (WGD), which can be generated in [EvidenceQC](#evidence-qc). An example batching process is outlined below:
1. Divide the cohort into PCR+ and PCR- samples
2. Partition the samples by median coverage from [EvidenceQC](#evidence-qc), grouping samples with similar median coverage together. The end goal is to divide the cohort into roughly equal-sized batches of about 100-500 samples; if your partitions based on coverage are larger or uneven, you can partition the cohort further in the next step to obtain the final batches. 
3. Optionally, divide the samples further by dosage score (WGD) from [EvidenceQC](#evidence-qc), grouping samples with similar WGD score together, to obtain roughly equal-sized batches of about 100-500 samples
4. Maintain a roughly equal sex balance within each batch, based on sex assignments from [EvidenceQC](#evidence-qc)


## <a name="sample-sample-mode">Single-sample mode</a>
`GATKSVPipelineSingleSample.wdl` runs the pipeline on a single sample using a fixed reference panel. An example run with reference panel containing 156 samples from the [NYGC 1000G Terra workspace](https://app.terra.bio/#workspaces/anvil-datastorage/1000G-high-coverage-2019) can be found in `inputs/build/NA12878/test` after [building inputs](#Building inputs)).

## <a name="gcnv-training-overview">gCNV Training</a>
Both the cohort and single-sample modes use the GATK gCNV depth calling pipeline, which requires a [trained model](#gcnv-training) as input. The samples used for training should be technically homogeneous and similar to the samples to be processed (i.e. same sample type, library prep protocol, sequencer, sequencing center, etc.). The samples to be processed may comprise all or a subset of the training set. For small, relatively homogenous cohorts, a single gCNV model is usually sufficient. If a cohort contains multiple data sources, we recommend training a separate model for each [batch](#batching) or group of batches with similar dosage score (WGD). The model may be trained on all or a subset of the samples to which it will be applied; a reasonable default is 100 randomly-selected samples from the batch (the random selection can be done as part of the workflow by specifying a number of samples to the `n_samples_subsample` input parameter in `/wdl/TrainGCNV.wdl`).

## <a name="reference-panel-generation">Generating a reference panel</a>
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

## <a name="descriptions">Module Descriptions</a>
The following sections briefly describe each module and highlights inter-dependent input/output files. Note that input/output mappings can also be gleaned from `GATKSVPipelineBatch.wdl`, and example input templates for each module can be found in `/inputs/templates/test`.

## <a name="gather-sample-evidence">GatherSampleEvidence</a>
*Formerly Module00a*

Runs raw evidence collection on each sample with the following SV callers: Manta, Wham, and/or MELT. Delly can be enabled but is no longer officially supported. For guidance on pre-filtering prior to `GatherSampleEvidence`, refer to the [Sample Exclusion](#sample-exclusion) section.

Note: a list of sample IDs must be provided. Refer to the [sample ID requirements](#sampleids) for specifications of allowable sample IDs. IDs that do not meet these requirements may cause errors.

#### Inputs:
* Per-sample BAM or CRAM files aligned to hg38. Index files (`.bai`) must be provided if using BAMs.

#### Outputs:
* Caller VCFs (Delly, Manta, MELT, and/or Wham)
* Binned read counts file
* Split reads (SR) file
* Discordant read pairs (PE) file

## <a name="evidence-qc">EvidenceQC</a>
*Formerly Module00b*

Runs ploidy estimation, dosage scoring, and optionally VCF QC. The results from this module can be used for QC and batching.

For large cohorts, this workflow can be run on arbitrary cohort partitions of up to about 500 samples. Afterwards, we recommend using the results to divide samples into smaller batches (~100-500 samples) with ~1:1 male:female ratio. Refer to the [Batching](#batching) section for further guidance on creating batches.

We also recommend using sex assignments generated from the ploidy estimates and incorporating them into the PED file, with sex = 0 for sex aneuploidies.

#### Prerequisites:
* [GatherSampleEvidence](#gather-sample-evidence)

#### Inputs:
* Read count files ([GatherSampleEvidence](#gather-sample-evidence))
* (Optional) SV call VCFs ([GatherSampleEvidence](#gather-sample-evidence))

#### Outputs:
* Per-sample dosage scores with plots
* Median coverage per sample
* Ploidy estimates, sex assignments, with plots
* (Optional) Outlier samples detected by call counts

#### <a name="prelim-sample-qc">Preliminary Sample QC</a>
The purpose of sample filtering at this stage after EvidenceQC is to prevent very poor quality samples from interfering with the results for the rest of the callset. In general, samples that are borderline are okay to leave in, but you should choose filtering thresholds to suit the needs of your cohort and study. There will be future opportunities (as part of [FilterBatch](#filter-batch)) for filtering before the joint genotyping stage if necessary. Here are a few of the basic QC checks that we recommend:
* Look at the X and Y ploidy plots, and check that sex assignments match your expectations. If there are discrepancies, check for sample swaps and update your PED file before proceeding.
* Look at the dosage score (WGD) distribution and check that it is centered around 0 (the distribution of WGD for PCR- samples is expected to be slightly lower than 0, and the distribution of WGD for PCR+ samples is expected to be slightly greater than 0. Refer to the [gnomAD-SV paper](https://doi.org/10.1038/s41586-020-2287-8) for more information on WGD score). Optionally filter outliers.
* Look at the low outliers for each SV caller (samples with much lower than typical numbers of SV calls per contig for each caller). An empty low outlier file means there were no outliers below the median and no filtering is necessary. Check that no samples had zero calls.
* Look at the high outliers for each SV caller and optionally filter outliers; samples with many more SV calls than average may be poor quality.


## <a name="gcnv-training">TrainGCNV</a>
Trains a gCNV model for use in [GatherBatchEvidence](#gather-batch-evidence). The WDL can be found at `/wdl/TrainGCNV.wdl`. See the [gCNV training overview](#gcnv-training-overview) for more information.

#### Prerequisites:
* [GatherSampleEvidence](#gather-sample-evidence)
* (Recommended) [EvidenceQC](#evidence-qc)

#### Inputs:
* Read count files ([GatherSampleEvidence](#gather-sample-evidence))

#### Outputs:
* Contig ploidy model tarball
* gCNV model tarballs


## <a name="gather-batch-evidence">GatherBatchEvidence</a>
*Formerly Module00c*

Runs CNV callers (cnMOPs, GATK gCNV) and combines single-sample raw evidence into a batch. See [above](#cohort-mode) for more information on batching.

#### Prerequisites:
* [GatherSampleEvidence](#gather-sample-evidence)
* (Recommended) [EvidenceQC](#evidence-qc)
* [gCNV training](#gcnv-training)

#### Inputs:
* PED file (updated with [EvidenceQC](#evidence-qc) sex assignments, including sex = 0 for sex aneuploidies. Calls will not be made on sex chromosomes when sex = 0 in order to avoid generating many confusing calls or upsetting normalized copy numbers for the batch.)
* Per-sample indexed GVCFs generated with HaplotypeCaller (`gvcfs` input), or a jointly-genotyped VCF (position-sharded, `snp_vcfs` input or `snp_vcfs_shard_list` input). The jointly-genotyped VCF may contain multi-allelic sites and indels, but only biallelic SNVs will be used by the pipeline. We recommend shards of 10 GB or less to lower compute time and resources.
* Read count, BAF, PE, and SR files ([GatherSampleEvidence](#gather-sample-evidence))
* Caller VCFs ([GatherSampleEvidence](#gather-sample-evidence))
* Contig ploidy model and gCNV model files ([gCNV training](#gcnv-training))

#### Outputs:
* Combined read count matrix, SR, PE, and BAF files
* Standardized call VCFs
* Depth-only (DEL/DUP) calls
* Per-sample median coverage estimates
* (Optional) Evidence QC plots


## <a name="cluster-batch">ClusterBatch</a>
*Formerly Module01*

Clusters SV calls across a batch.

#### Prerequisites:
* [GatherBatchEvidence](#gather-batch-evidence)

#### Inputs:
* Standardized call VCFs ([GatherBatchEvidence](#gather-batch-evidence))
* Depth-only (DEL/DUP) calls ([GatherBatchEvidence](#gather-batch-evidence))

#### Outputs:
* Clustered SV VCFs
* Clustered depth-only call VCF


## <a name="generate-batch-metrics">GenerateBatchMetrics</a>
*Formerly Module02*

Generates variant metrics for filtering.

#### Prerequisites:
* [ClusterBatch](#cluster-batch)

#### Inputs:
* Combined read count matrix, SR, PE, and BAF files ([GatherBatchEvidence](#gather-batch-evidence))
* Per-sample median coverage estimates ([GatherBatchEvidence](#gather-batch-evidence))
* Clustered SV VCFs ([ClusterBatch](#cluster-batch))
* Clustered depth-only call VCF ([ClusterBatch](#cluster-batch))

#### Outputs:
* Metrics file


## <a name="generate-batch-metrics">FilterBatch</a>
*Formerly Module03*

Filters poor quality variants and filters outlier samples. This workflow can be run all at once with the WDL at `wdl/FilterBatch.wdl`, or it can be run in three steps to enable tuning of outlier filtration cutoffs. The three subworkflows are:
1. FilterBatchSites: Per-batch variant filtration
2. PlotSVCountsPerSample: Visualize SV counts per sample per type to help choose an IQR cutoff for outlier filtering, and preview outlier samples for a given cutoff
3. FilterBatchSamples: Per-batch outlier sample filtration; provide an appropriate `outlier_cutoff_nIQR` based on the SV count plots and outlier previews from step 2.

#### Prerequisites:
* [GenerateBatchMetrics](#generate-batch-metrics)

#### Inputs:
* Batch PED file
* Metrics file ([GenerateBatchMetrics](#generate-batch-metrics))
* Clustered SV and depth-only call VCFs ([ClusterBatch](#cluster-batch))

#### Outputs:
* Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
* Filtered depth-only call VCF with outlier samples excluded
* Random forest cutoffs file
* PED file with outlier samples excluded


## <a name="merge-batch-sites">MergeBatchSites</a>
*Formerly MergeCohortVcfs*

Combines filtered variants across batches. The WDL can be found at: `/wdl/MergeBatchSites.wdl`.

#### Prerequisites:
* [FilterBatch](#filter-batch)

#### Inputs:
* List of filtered PESR VCFs ([FilterBatch](#filter-batch))
* List of filtered depth VCFs ([FilterBatch](#filter-batch))

#### Outputs:
* Combined cohort PESR and depth VCFs


## <a name="genotype-batch">GenotypeBatch</a>
*Formerly Module04*

Genotypes a batch of samples across unfiltered variants combined across all batches.

#### Prerequisites:
* [FilterBatch](#filter-batch)
* [MergeBatchSites](#merge-batch-sites)

#### Inputs:
* Batch PESR and depth VCFs ([FilterBatch](#filter-batch))
* Cohort PESR and depth VCFs ([MergeBatchSites](#merge-batch-sites))
* Batch read count, PE, and SR files ([GatherBatchEvidence](#gather-batch-evidence)) 

#### Outputs:
* Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
* Filtered depth-only call VCF with outlier samples excluded
* PED file with outlier samples excluded
* List of SR pass variants
* List of SR fail variants
* (Optional) Depth re-genotyping intervals list


## <a name="regenotype-cnvs">RegenotypeCNVs</a>
*Formerly Module04b*

Re-genotypes probable mosaic variants across multiple batches.

#### Prerequisites:
* [GenotypeBatch](#genotype-batch)

#### Inputs:
* Per-sample median coverage estimates ([GatherBatchEvidence](#gather-batch-evidence))
* Pre-genotyping depth VCFs ([FilterBatch](#filter-batch))
* Batch PED files ([FilterBatch](#filter-batch))
* Cohort depth VCF ([MergeBatchSites](#merge-batch-sites))
* Genotyped depth VCFs ([GenotypeBatch](#genotype-batch))
* Genotyped depth RD cutoffs file ([GenotypeBatch](#genotype-batch))

#### Outputs:
* Re-genotyped depth VCFs


## <a name="make-cohort-vcf">MakeCohortVcf</a>
*Formerly Module0506*

Combines variants across multiple batches, resolves complex variants, re-genotypes, and performs final VCF clean-up.

#### Prerequisites:
* [GenotypeBatch](#genotype-batch)
* (Optional) [RegenotypeCNVs](#regenotype-cnvs)

#### Inputs:
* RD, PE and SR file URIs ([GatherBatchEvidence](#gather-batch-evidence))
* Batch filtered PED file URIs ([FilterBatch](#filter-batch))
* Genotyped PESR VCF URIs ([GenotypeBatch](#genotype-batch))
* Genotyped depth VCF URIs ([GenotypeBatch](#genotype-batch) or [RegenotypeCNVs](#regenotype-cnvs))
* SR pass variant file URIs ([GenotypeBatch](#genotype-batch))
* SR fail variant file URIs ([GenotypeBatch](#genotype-batch))
* Genotyping cutoff file URIs ([GenotypeBatch](#genotype-batch))
* Batch IDs
* Sample ID list URIs

#### Outputs:
* Finalized "cleaned" VCF and QC plots

## <a name="module07">Module 07</a> (in development)
Apply downstream filtering steps to the cleaned vcf to further control the false discovery rate; all steps are optional and users should decide based on the specific purpose of their projects.

Filterings methods include:
* minGQ - remove variants based on the genotype quality across populations.
Note: Trio families are required to build the minGQ filtering model in this step. We provide tables pre-trained with the 1000 genomes samples at different FDR thresholds for projects that lack family structures, and they can be found at the paths below.  These tables assume that GQ has a scale of [0,999], so they will not work with newer VCFs where GQ has a scale of [0,99].
```
gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.10perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.1perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.5perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
```

* BatchEffect - remove variants that show significant discrepancies in allele frequencies across batches
* FilterOutlierSamplesPostMinGQ - remove outlier samples with unusually high or low number of SVs
* FilterCleanupQualRecalibration - sanitize filter columns and recalibrate variant QUAL scores for easier interpretation

## <a name="annotate-vcf">AnnotateVcf</a> (in development)
*Formerly Module08Annotation*

Add annotations, such as the inferred function and allele frequencies of variants, to final vcf.

Annotations methods include:
* Functional annotation - annotate SVs with inferred function on protein coding regions, regulatory regions such as UTR and Promoters and other non coding elements;
* Allele Frequency annotation - annotate SVs with their allele frequencies across all samples, and samples of specific sex, as well as specific sub-populations.
* Allele Frequency annotation with external callset - annotate SVs with the allele frequencies of their overlapping SVs in another callset, eg. gnomad SV callset.

## <a name="module09">Module 09</a> (in development)
Visualize SVs with [IGV](http://software.broadinstitute.org/software/igv/) screenshots and read depth plots.

Visualization methods include:
* RD Visualization - generate RD plots across all samples, ideal for visualizing large CNVs.
* IGV Visualization - generate IGV plots of each SV for individual sample, ideal for visualizing de novo small SVs.
* Module09.visualize.wdl - generate RD plots and IGV plots, and combine them for easy review.

## CI/CD
This repository is maintained following the norms of 
continuous integration (CI) and continuous delivery (CD). 
GATK-SV CI/CD is developed as a set of Github Actions
workflows that are available under the `.github/workflows`
directory. Please refer to the [workflow's README](.github/workflows/README.md) 
for their current coverage and setup. 

## <a name="troubleshooting">Troubleshooting</a>

### VM runs out of memory or disk
* Default pipeline settings are tuned for batches of 100 samples. Larger batches or cohorts may require additional VM resources. Most runtime attributes can be modified through the `RuntimeAttr` inputs. These are formatted like this in the json:
```
"MyWorkflow.runtime_attr_override": {
  "disk_gb": 100,
  "mem_gb": 16
},
```
Note that a subset of the struct attributes can be specified. See `wdl/Structs.wdl` for available attributes.


### Calculated read length causes error in MELT workflow

Example error message from `GatherSampleEvidence.MELT.GetWgsMetrics`:
```
Exception in thread "main" java.lang.ArrayIndexOutOfBoundsException: The requested index 701766 is out of counter bounds. Possible cause of exception can be wrong READ_LENGTH parameter (much smaller than actual read length)
```

This error message was observed for a sample with an average read length of 117, but for which half the reads were of length 90 and half were of length 151. As a workaround, override the calculated read length by providing a `read_length` input of 151 (or the expected read length for the sample in question) to `GatherSampleEvidence`.
