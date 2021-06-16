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
* [Module Descriptions](#descriptions)
    * [Module 00a](#module00a) - Raw callers and evidence collection
    * [Module 00b](#module00b) - Batch QC
    * [gCNV training](#gcnv-training) - gCNV model creation
    * [Module 00c](#module00c) - Batch evidence merging, BAF generation, and depth callers
    * [Module 01](#module01) - Site clustering
    * [Module 02](#module02) - Site metrics
    * [Module 03](#module03) - Filtering
    * [Gather Cohort VCFs](#gather-vcfs) - Cross-batch site merging
    * [Module 04](#module04) - Genotyping
    * [Module 04b](#module04b) - Genotype refinement (optional)
    * [Module 05/06](#module0506) - Cross-batch integration, complex event resolution, and VCF cleanup
    * [Module 07](#module07) - Downstream Filtering
    * [Module 08](#module08) - Annotation
    * [Module 09](#module09) - QC and Visualization
    * Additional modules - Mosaic and de novo
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

### Data:
* Illumina short-read whole-genome CRAMs or BAMs, aligned to hg38 with [bwa-mem](https://github.com/lh3/bwa). BAMs must also be indexed.
* Indexed GVCFs produced by GATK HaplotypeCaller, or a jointly genotyped VCF.
* Family structure definitions file in [PED format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format). Sex aneuploidies (detected in [Module 00b](#module00b)) should be entered as sex = 0.

#### <a name="sampleids">Sample ID requirements:</a>

Sample IDs must:
* Be unique within the cohort
* Contain only alphanumeric characters and underscores (no dashes, whitespace, or special characters)

Sample IDs should not:
* Contain only numeric characters
* Be a substring of another sample ID in the same cohort
* Contain any of the following substrings: `chr`, `name`, `DEL`, `DUP`, `CPX`, `CHROM`

The same requirements apply to family IDs in the PED file, as well as batch IDs and the cohort ID provided as workflow inputs.

Sample IDs are provided to [Module00a](#module00a) directly and need not match sample names from the BAM/CRAM headers or GVCFs. `GetSampleID.wdl` can be used to fetch BAM sample IDs and also generates a set of alternate IDs that are considered safe for this pipeline; alternatively, [this script](https://github.com/talkowski-lab/gnomad_sv_v3/blob/master/sample_id/convert_sample_ids.py) transforms a list of sample IDs to fit these requirements. Currently, sample IDs can be replaced again in [Module 00c](#module00c). 

The following inputs will need to be updated with the transformed sample IDs:
* Sample ID list for [Module00a](#module00a) or [Module 00c](#module00c)
* PED file

If using a SNP VCF in [Module 00c](#module00c), it does not need to be re-headered; simply provide the `vcf_samples` argument.


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

#### Inputs
Example workflow inputs can be found in `/inputs`. All required resources are available in public Google buckets.

#### MELT
**Important**: The example input files contain MELT inputs that are NOT public (see [Requirements](#requirements)). These include:

* `GATKSVPipelineSingleSample.melt_docker` and `GATKSVPipelineBatch.melt_docker` - MELT docker URI (see [Docker readme](https://github.com/talkowski-lab/gatk-sv-v1/blob/master/dockerfiles/README.md))
* `GATKSVPipelineSingleSample.ref_std_melt_vcfs` - Standardized MELT VCFs ([Module00c](#module00c))

The input values are provided only as an example and are not publicly accessible. In order to include MELT, these values must be provided by the user. MELT can be disabled by deleting these inputs and setting `GATKSVPipelineBatch.use_melt` to `false`.

#### Requester pays buckets
**Important**: The following parameters must be set when certain input data is in requester pays (RP) buckets:

* `GATKSVPipelineSingleSample.requester_pays_cram` and `GATKSVPipelineBatch.Module00aBatch.requester_pays_crams` - set to `True` if inputs are CRAM format and in an RP bucket, otherwise `False`.
* `GATKSVPipelineBatch.GATKSVPipelinePhase1.gcs_project_for_requester_pays` - set to your Google Cloud Project ID if gVCFs are in an RP bucket, otherwise omit this parameter.

#### Execution
We recommend running the pipeline on a dedicated [Cromwell](https://github.com/broadinstitute/cromwell) server with a [cromshell](https://github.com/broadinstitute/cromshell) client. A batch run can be started with the following commands:

```
> mkdir gatksv_run && cd gatksv_run
> mkdir wdl && cd wdl
> cp $GATK_SV_V1_ROOT/wdl/*.wdl .
> zip dep.zip *.wdl
> cd ..
> cp $GATK_SV_V1_ROOT/inputs/GATKSVPipelineBatch.ref_panel_1kg.json GATKSVPipelineBatch.my_run.json
> cromshell submit wdl/GATKSVPipelineBatch.wdl GATKSVPipelineBatch.my_run.json cromwell_config.json wdl/dep.zip
```

where `cromwell_config.json` is a Cromwell [workflow options file](https://cromwell.readthedocs.io/en/stable/wf_options/Overview/). Note users will need to re-populate batch/sample-specific parameters (e.g. BAMs and sample IDs).

## <a name="overview">Pipeline Overview</a>
The pipeline consists of a series of modules that perform the following:
* [Module 00a](#module00a): SV evidence collection, including calls from a configurable set of algorithms (Delly, Manta, MELT, and Wham), read depth (RD), split read positions (SR), and discordant pair positions (PE).
* [Module 00b](#module00b): Dosage bias scoring and ploidy estimation
* [Module 00c](#module00c): Copy number variant calling using cn.MOPS and GATK gCNV; B-allele frequency (BAF) generation; call and evidence aggregation
* [Module 01](#module01): Variant clustering
* [Module 02](#module02): Variant filtering metric generation
* [Module 03](#module03): Variant filtering; outlier exclusion
* [Module 04](#module04): Genotyping
* [Module 05/06](#module0506): Cross-batch integration; complex variant resolution and re-genotyping; vcf cleanup
* [Module 07](#module07): Downstream filtering, including minGQ, batch effect check, outlier samples removal and final recalibration;
* [Module 08](#module08): Annotations, including functional annotation, allele frequency (AF) annotation and AF annotation with external population callsets;
* [Module 09](#module09): Visualization, including scripts that generates IGV screenshots and rd plots.
* Additional modules to be added: de novo and mosaic scripts 


Repository structure:
* `/inputs`: Example workflow parameter files for running gCNV training, GATK-SV batch mode, and GATK-SV single-sample mode
* `/dockerfiles`: Resources for building pipeline docker images (see [readme](https://github.com/talkowski-lab/gatk-sv-v1/blob/master/dockerfiles/README.md))
* `/wdl`: WDLs running the pipeline. There is a master WDL for running each module, e.g. `Module01.wdl`.
* `/scripts`: scripts for running tests, building dockers, and analyzing cromwell metadata files
* `/src`: main pipeline scripts
  * `/RdTest`: scripts for depth testing
  * `/sv-pipeline`: various scripts and packages used throughout the pipeline
  * `/svqc`: Python module for checking that pipeline metrics fall within acceptable limits
  * `/svtest`: Python module for generating various summary metrics from module outputs
  * `/svtk`: Python module of tools for SV-related datafile parsing and analysis
  * `/WGD`: whole-genome dosage scoring scripts
* `/test`: WDL test parameter files. Please note that file inputs may not be publicly available.


## <a name="cohort-mode">Cohort mode</a>
A minimum cohort size of 100 with roughly equal number of males and females is recommended. For modest cohorts (~100-500 samples), the pipeline can be run as a single batch using `GATKSVPipelineBatch.wdl`.

For larger cohorts, samples should be split up into batches of ~100-500 samples. We recommend batching based on overall coverage and dosage score (WGD), which can be generated in [Module 00b](#module00b). 

The pipeline should be executed as follows:
* Modules [00a](#module00a) and [00b](#module00b) can be run on arbitrary cohort partitions
* Modules [00c](#module00c), [01](#module01), [02](#module02), and [03](#module03) are run separately per batch
* [Module 04](#module04) is run separately per batch, using filtered variants ([Module 03](#module03) output) combined across all batches
* [Module 05/06](#module0506) and beyond are run on all batches together

Note: [Module 00c](#module00c) requires a [trained gCNV model](#gcnv-training).


## <a name="sample-sample-mode">Single-sample mode</a>
`GATKSVPipelineSingleSample.wdl` runs the pipeline on a single sample using a fixed reference panel. An example reference panel containing 156 samples from the [NYGC 1000G Terra workspace](https://app.terra.bio/#workspaces/anvil-datastorage/1000G-high-coverage-2019) is provided with `inputs/GATKSVPipelineSingleSample.ref_panel_1kg.na12878.json`. 

Custom reference panels can be generated by running `GATKSVPipelineBatch.wdl` and `trainGCNV.wdl` and using the outputs to replace the following single-sample workflow inputs:

* `GATKSVPipelineSingleSample.ref_ped_file` : `batch.ped` - Manually created (see [data requirements](#requirements))
* `GATKSVPipelineSingleSample.contig_ploidy_model_tar` : `batch-contig-ploidy-model.tar.gz` - gCNV contig ploidy model ([gCNV training](#gcnv-training))
* `GATKSVPipelineSingleSample.gcnv_model_tars` : `batch-model-files-*.tar.gz` - gCNV model tarballs ([gCNV training](#gcnv-training))
* `GATKSVPipelineSingleSample.ref_pesr_disc_files` - `sample.disc.txt.gz` - Paired-end evidence files ([Module 00a](#module00a))
* `GATKSVPipelineSingleSample.ref_pesr_split_files` - `sample.split.txt.gz` - Split read evidence files ([Module 00a](#module00a))
* `GATKSVPipelineSingleSample.ref_panel_bincov_matrix`: `batch.RD.txt.gz` - Read counts matrix ([Module 00c](#module00c))
* `GATKSVPipelineSingleSample.ref_panel_del_bed` : `batch.DEL.bed.gz` - Depth deletion calls ([Module 00c](#module00c))
* `GATKSVPipelineSingleSample.ref_panel_dup_bed` : `batch.DUP.bed.gz` - Depth duplication calls ([Module 00c](#module00c))
* `GATKSVPipelineSingleSample.ref_samples` - Reference panel sample IDs
* `GATKSVPipelineSingleSample.ref_std_manta_vcfs` - `std_XXX.manta.sample.vcf.gz` - Standardized Manta VCFs ([Module 00c](#module00c))
* `GATKSVPipelineSingleSample.ref_std_melt_vcfs` - `std_XXX.melt.sample.vcf.gz` - Standardized Melt VCFs ([Module 00c](#module00c))
* `GATKSVPipelineSingleSample.ref_std_wham_vcfs` - `std_XXX.wham.sample.vcf.gz` - Standardized Wham VCFs ([Module 00c](#module00c))
* `GATKSVPipelineSingleSample.cutoffs` : `batch.cutoffs` - Filtering cutoffs ([Module 03](#module03))
* `GATKSVPipelineSingleSample.genotype_pesr_pesr_sepcutoff` : `genotype_pesr.pesr_sepcutoff.txt` - Genotyping cutoffs ([Module 04](#module04))
* `GATKSVPipelineSingleSample.genotype_pesr_depth_sepcutoff` : `genotype_pesr.depth_sepcutoff.txt` - Genotyping cutoffs ([Module 04](#module04))
* `GATKSVPipelineSingleSample.genotype_depth_pesr_sepcutoff` : `genotype_depth.pesr_sepcutoff.txt` - Genotyping cutoffs ([Module 04](#module04))
* `GATKSVPipelineSingleSample.genotype_depth_depth_sepcutoff` : `genotype_depth.depth_sepcutoff.txt` - Genotyping cutoffs ([Module 04](#module04))
* `GATKSVPipelineSingleSample.PE_metrics` : `pe_metric_file.txt` - Paired-end evidence genotyping metrics ([Module 04](#module04))
* `GATKSVPipelineSingleSample.SR_metrics` : `sr_metric_file.txt` - Split read evidence genotyping metrics ([Module 04](#module04))
* `GATKSVPipelineSingleSample.ref_panel_vcf` : `batch.cleaned.vcf.gz` - Final output VCF ([Module 05/06](#module0506))


## <a name="gcnv-training-overview">gCNV Training</a>
Both the cohort and single-sample modes use the GATK gCNV depth calling pipeline, which requires a [trained model](#gcnv-training) as input. The samples used for training should be technically homogeneous and similar to the samples to be processed (i.e. same sample type, library prep protocol, sequencer, sequencing center, etc.). The samples to be processed may comprise all or a subset of the training set. For small cohorts, a single gCNV model is usually sufficient. If a cohort contains multiple data sources, we recommend clustering them using the dosage score, and training a separate model for each cluster.


## <a name="descriptions">Module Descriptions</a>
The following sections briefly describe each module and highlights inter-dependent input/output files. Note that input/output mappings can also be gleaned from `GATKSVPipelineBatch.wdl`, and example input files for each module can be found in `/test`.

## <a name="module00a">Module 00a</a>
Runs raw evidence collection on each sample.

Note: a list of sample IDs must be provided. Refer to the [sample ID requirements](#sampleids) for specifications of allowable sample IDs. IDs that do not meet these requirements may cause errors.

#### Inputs:
* Per-sample BAM or CRAM files aligned to hg38. Index files (`.bai`) must be provided if using BAMs.

#### Outputs:
* Caller VCFs (Delly, Manta, MELT, and/or Wham)
* Binned read counts file
* Split reads (SR) file
* Discordant read pairs (PE) file
* B-allele fraction (BAF) file


## <a name="module00b">Module 00b</a>
Runs ploidy estimation, dosage scoring, and optionally VCF QC. The results from this module can be used for QC and batching.

For large cohorts, we recommend dividing samples into smaller batches (~500 samples) with ~1:1 male:female ratio.

We also recommend using sex assignments generated from the ploidy estimates and incorporating them into the PED file.

#### Prerequisites:
* [Module 00a](#module00a)

#### Inputs:
* Read count files ([Module 00a](#module00a))
* (Optional) SV call VCFs ([Module 00a](#module00a))

#### Outputs:
* Per-sample dosage scores with plots
* Ploidy estimates, sex assignments, with plots
* (Optional) Outlier samples detected by call counts


## <a name="gcnv-training">gCNV Training</a>
Trains a gCNV model for use in [Module 00c](#module00c). The WDL can be found at `/gcnv/trainGCNV.wdl`.

#### Prerequisites:
* [Module 00a](#module00a)
* (Recommended) [Module 00b](#module00b)

#### Inputs:
* Read count files ([Module 00a](#module00a))

#### Outputs:
* Contig ploidy model tarball
* gCNV model tarballs


## <a name="module00c">Module 00c</a>
Runs CNV callers (cnMOPs, GATK gCNV) and combines single-sample raw evidence into a batch. See [above]("#cohort-mode") for more information on batching.

#### Prerequisites:
* [Module 00a](#module00a)
* (Recommended) [Module 00b](#module00b)
* gCNV training

#### Inputs:
* PED file (updated with [Module 00b](#module00b) sex assignments, including sex = 0 for sex aneuploidies. Calls will not be made on sex chromosomes when sex = 0 in order to avoid generating many confusing calls or upsetting normalized copy numbers for the batch.)
* Per-sample GVCFs generated with HaplotypeCaller (`gvcfs` input), or a jointly-genotyped VCF (position-sharded, `snp_vcfs` input)
* Read count, BAF, PE, and SR files ([Module 00a](#module00a))
* Caller VCFs ([Module 00a](#module00a))
* Contig ploidy model and gCNV model files (gCNV training)

#### Outputs:
* Combined read count matrix, SR, PE, and BAF files
* Standardized call VCFs
* Depth-only (DEL/DUP) calls
* Per-sample median coverage estimates
* (Optional) Evidence QC plots


## <a name="module01">Module 01</a>
Clusters SV calls across a batch.

#### Prerequisites:
* [Module 00c](#module00c)

#### Inputs:
* Standardized call VCFs ([Module 00c](#module00c))
* Depth-only (DEL/DUP) calls ([Module 00c](#module00c))

#### Outputs:
* Clustered SV VCFs
* Clustered depth-only call VCF


## <a name="module02">Module 02</a>
Generates variant metrics for filtering.

#### Prerequisites:
* [Module 01](#module01)

#### Inputs:
* Combined read count matrix, SR, PE, and BAF files ([Module 00c](#module00c))
* Per-sample median coverage estimates ([Module 00c](#module00c))
* Clustered SV VCFs ([Module 01](#module01))
* Clustered depth-only call VCF ([Module 01](#module01))

#### Outputs:
* Metrics file


## <a name="module02">Module 03</a>
Filters poor quality variants and filters outlier samples.

#### Prerequisites:
* [Module 02](#module02)

#### Inputs:
* Batch PED file
* Metrics file ([Module 02](#module02))
* Clustered SV and depth-only call VCFs ([Module 01](#module01))

#### Outputs:
* Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
* Filtered depth-only call VCF with outlier samples excluded
* Random forest cutoffs file
* PED file with outlier samples excluded


## <a name="module04">Merge Cohort VCFs</a>
Combines filtered variants across batches. The WDL can be found at: `/wdl/MergeCohortVcfs.wdl`.

#### Prerequisites:
* [Module 03](#module03)

#### Inputs:
* List of filtered PESR VCFs ([Module 03](#module03))
* List of filtered depth VCFs ([Module 03](#module03))

#### Outputs:
* Combined cohort PESR and depth VCFs
* Cohort and clustered depth variant BED files


## <a name="module04">Module 04</a>
Genotypes a batch of samples across unfiltered variants combined across all batches.

#### Prerequisites:
* [Module 03](#module03)
* Merge Cohort VCFs

#### Inputs:
* Batch PESR and depth VCFs ([Module 03](#module03))
* Cohort PESR and depth VCFs (Merge Cohort VCFs)
* Batch read count, PE, and SR files ([Module 00c](#module00c)) 

#### Outputs:
* Filtered SV (non-depth-only a.k.a. "PESR") VCF with outlier samples excluded
* Filtered depth-only call VCF with outlier samples excluded
* PED file with outlier samples excluded
* List of SR pass variants
* List of SR fail variants
* (Optional) Depth re-genotyping intervals list


## <a name="module04b">Module 04b</a>
Re-genotypes probable mosaic variants across multiple batches.

#### Prerequisites:
* [Module 04](#module04)

#### Inputs:
* Per-sample median coverage estimates ([Module 00c](#module00c))
* Pre-genotyping depth VCFs ([Module 03](#module03))
* Batch PED files ([Module 03](#module03))
* Clustered depth variant BED file (Merge Cohort VCFs)
* Cohort depth VCF (Merge Cohort VCFs)
* Genotyped depth VCFs ([Module 04](#module04))
* Genotyped depth RD cutoffs file ([Module 04](#module04))

#### Outputs:
* Re-genotyped depth VCFs


## <a name="module0506">Module 05/06</a>
Combines variants across multiple batches, resolves complex variants, re-genotypes, and performs final VCF clean-up.

#### Prerequisites:
* [Module 04](#module04)
* (Optional) [Module 04b](#module04b)

#### Inputs:
* RD, PE and SR file URIs ([Module 00c](#module00c))
* Batch filtered PED file URIs ([Module 03](#module03))
* Genotyped PESR VCF URIs ([Module 04](#module04))
* Genotyped depth VCF URIs ([Module 04](#module04) or [04b](#module04b))
* SR pass variant file URIs ([Module 04](#module04))
* SR fail variant file URIs ([Module 04](#module04))
* Genotyping cutoff file URIs ([Module 04](#module04))
* Batch IDs
* Sample ID list URIs

#### Outputs:
* Finalized "cleaned" VCF and QC plots

## <a name="module07">Module 07</a> (in development)
Apply downstream filtering steps to the cleaned vcf to further control the false discovery rate; all steps are optional and users should decide based on the specific purpose of their projects.

Filterings methods include:
* minGQ - remove variants based on the genotype quality across populations.
Note: Trio families are required to build the minGQ filtering model in this step. We provide tables pre-trained with the 1000 genomes samples at different FDR thresholds for projects that lack family structures, and they can be found here: 
```
gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.10perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.1perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
gs://gatk-sv-resources-public/hg38/v0/sv-resources/ref-panel/1KG/v2/mingq/1KGP_2504_and_698_with_GIAB.5perc_fdr.PCRMINUS.minGQ.filter_lookup_table.txt
```

* BatchEffect - remove variants that show significant discrepancies in allele frequencies across batches
* FilterOutlierSamples - remove outlier samples with unusually high or low number of SVs
* FilterCleanupQualRecalibration - sanitize filter columns and recalibrate variant QUAL scores for easier interpretation

## <a name="module08">Module 08</a> (in development)
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


## <a name="troubleshooting">Troubleshooting</a>

### VM runs out of memory or disk
* Default pipeline settings are tuned for batches of 100 samples. Larger batches or cohorts may require additional VM resources. Most runtime attributes can be modified through the `RuntimeAttr` inputs. These are formatted like this in the json:
```
"ModuleX.runtime_attr_override": {
  "disk_gb": 100,
  "mem_gb": 16
},
```
Note that a subset of the struct attributes can be specified. See `wdl/Structs.wdl` for available attributes.
