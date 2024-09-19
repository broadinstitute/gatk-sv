---
title: GatherBatchEvidence
description: Gather Batch Evidence
sidebar_position: 4
slug: gbe
---

Runs CNV callers ([cn.MOPS](https://academic.oup.com/nar/article/40/9/e69/1136601), GATK-gCNV) 
and combines single-sample raw evidence into a batch.

The following diagram illustrates the downstream workflows of the `GatherBatchEvidence` workflow 
in the recommended invocation order. You may refer to 
[this diagram](https://github.com/broadinstitute/gatk-sv/blob/main/terra_pipeline_diagram.jpg) 
for the overall recommended invocation order.

```mermaid

stateDiagram
  direction LR
  
  classDef inModules stroke-width:0px,fill:#00509d,color:#caf0f8
  classDef thisModule font-weight:bold,stroke-width:0px,fill:#ff9900,color:white
  classDef outModules stroke-width:0px,fill:#caf0f8,color:#00509d

  gbe: GatherBatchEvidence
  t: TrainGCNV
  cb: ClusterBatch
  t --> gbe
  gbe --> cb
  
  class gbe thisModule
  class t inModules
  class cb outModules
```

## Inputs
This workflow takes as input the read counts, BAF, PE, SD, SR, and per-caller VCF files 
produced in the GatherSampleEvidence workflow, and contig ploidy and gCNV models from 
the TrainGCNV workflow.
The following is the list of the inputs the GatherBatchEvidence workflow takes.


#### `batch`
An identifier for the batch.


#### `samples`
Sets the list of sample IDs. 


#### `counts`
Set to the [`GatherSampleEvidence.coverage_counts`](./gse#coverage-counts) output.


#### Raw calls

The following inputs set the per-caller raw SV calls, and should be set 
if the caller was run in the [`GatherSampleEvidence`](./gse) workflow.
You may set each of the following inputs to the linked output from 
the GatherSampleEvidence workflow.


- `manta_vcfs`: [`GatherSampleEvidence.manta_vcf`](./gse#manta-vcf);
- `melt_vcfs`: [`GatherSampleEvidence.melt_vcf`](./gse#melt-vcf);
- `scramble_vcfs`: [`GatherSampleEvidence.scramble_vcf`](./gse#scramble-vcf);
- `wham_vcfs`: [`GatherSampleEvidence.wham_vcf`](./gse#wham-vcf).

#### `PE_files`
Set to the [`GatherSampleEvidence.pesr_disc`](./gse#pesr-disc) output.

#### `SR_files`
Set to the [`GatherSampleEvidence.pesr_split`](./gse#pesr-split)


#### `SD_files`
Set to the [`GatherSampleEvidence.pesr_sd`](./gse#pesr-sd)


#### `matrix_qc_distance`
You may refer to [this file](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/templates/terra_workspaces/cohort_mode/workflow_configurations/GatherBatchEvidence.json.tmpl)
for an example value. 


#### `min_svsize`
Sets the minimum size of SVs to include.


#### `ped_file`
A pedigree file describing the familial relationshipts between the samples in the cohort.
Please refer to [this section](./#ped_file) for details. 


#### `run_matrix_qc`
Enables or disables running optional QC tasks. 


#### `gcnv_qs_cutoff`
You may refer to [this file](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/templates/terra_workspaces/cohort_mode/workflow_configurations/GatherBatchEvidence.json.tmpl)
for an example value. 

#### cn.MOPS files
The workflow needs the following cn.MOPS files.

- `cnmops_chrom_file` and `cnmops_allo_file`: FASTA index files (`.fai`) for respectively 
  non-sex chromosomes (autosomes) and chromosomes X and Y (allosomes). 
  The file format is explained [on this page](https://www.htslib.org/doc/faidx.html).

  You may use the following files for these fields:

  ```json
  "cnmops_chrom_file": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/autosome.fai"
  "cnmops_allo_file": "gs://gcp-public-data--broad-references/hg38/v0/sv-resources/resources/v1/allosome.fai"
  ```
  
- `cnmops_exclude_list`: 
  You may use [this file](https://github.com/broadinstitute/gatk-sv/blob/d66f760865a89f30dbce456a3f720dec8b70705c/inputs/values/resources_hg38.json#L10)
  for this field.

#### GATK-gCNV inputs

The following inputs are configured based on the outputs generated in the [`TrainGCNV`](./gcnv) workflow.

- `contig_ploidy_model_tar`: [`TrainGCNV.cohort_contig_ploidy_model_tar`](./gcnv#contig-ploidy-model-tarball)
- `gcnv_model_tars`: [`TrainGCNV.cohort_gcnv_model_tars`](./gcnv#model-tarballs)


The workflow also enables setting a few optional arguments of gCNV.
The arguments and their default values are provided 
[here](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/templates/terra_workspaces/cohort_mode/workflow_configurations/GatherBatchEvidence.json.tmpl) 
as the following, and each argument is documented on 
[this page](https://gatk.broadinstitute.org/hc/en-us/articles/360037593411-PostprocessGermlineCNVCalls)
and
[this page](https://gatk.broadinstitute.org/hc/en-us/articles/360047217671-GermlineCNVCaller).


#### Docker images

The workflow needs the following Docker images, the latest versions of which are in 
[this file](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/dockers.json).

  - `cnmops_docker`;
  - `condense_counts_docker`;
  - `linux_docker`;
  - `sv_base_docker`;
  - `sv_base_mini_docker`;
  - `sv_pipeline_docker`;
  - `sv_pipeline_qc_docker`;
  - `gcnv_gatk_docker`;
  - `gatk_docker`.

#### Static inputs

You may refer to [this reference file](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json)
for values of the following inputs.

 - `primary_contigs_fai`;
 - `cytoband`;
 - `ref_dict`;
 - `mei_bed`;
 - `genome_file`;
 - `sd_locs_vcf`.


#### Optional Inputs
The following is the list of a few optional inputs of the 
workflow, with an example of possible values. 

- `"allosomal_contigs": [["chrX", "chrY"]]`
- `"ploidy_sample_psi_scale": 0.001`





## Outputs

- Combined read count matrix, SR, PE, and BAF files
- Standardized call VCFs
- Depth-only (DEL/DUP) calls
- Per-sample median coverage estimates
- (Optional) Evidence QC plots
