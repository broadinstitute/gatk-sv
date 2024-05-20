---
title: DeNovoSV
description: Downstream filtering (work in progress)
sidebar_position: 15
slug: denovo
---

The de-novo workflow operates on the annotated multi-sample VCF file created by 
the [AnnotateVCF](./av) workflow. 


### Inputs

- `vcf_file`: output of [AnnotateVcf](./av) called output_vcf.
  Note thatAll families in the vcf file must be included in the pedigree file


- `ped_input`: Must have a header as follows:

  | FamID | IndividualID | FatherID | MotherID | Gender    | Affected |
  |-|-|-|-|-|-|

- `genomic_disorder_input`: a file in BED format that contains regions of genomic disorder; 
   variants that overlap these regions will not be removed from the input VCF file. 

- `contigs`: Should be set to the following list.

  ```
  [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX" ]`
  ```

- `python_config`: a text file as the following.

  ```text
  large_cnv_size: '1000'
  gnomad_col: 'gnomAD_V2_AF'
  alt_gnomad_col: 'gnomad_v2.1_sv_AF'
  gnomad_AF: '0.01'
  parents_AF: '0.05'
  large_raw_overlap: '0.5'
  small_raw_overlap: '0.5'
  cohort_AF: '0.05'
  coverage_cutoff: '10'
  depth_only_size: '10000'
  parents_overlap: '0.5'
  gq_min: '0'
  ```
  
  Note that you value may increase the value of `cohort_AF` if the cohort is small.

- `batch_raw_file`: 
  a txt file with first column as batch and second column raw file generated from 
  module05-ClusterBatch for all callers except depth (clustered_manta_vcf, clustered_melt_vcf, clustered_wham_vcf).

  - Must match batch names in sample_batches (see below)
  - These batches and the samples contained in them are relevant in regards to the bincov matrices and raw files
  - If you are using the fam_ids input to only run the de novo script on certain samples, you can still include all raw files for this cohort and the script will only localize that files containing your samples of interest

- `batch_depth_raw_file`:
  a txt file with first column as batch and second column raw file generated from module05-ClusterBatch for depth caller (clustered_depth_vcf)
 
  - Must match batch names in sample_batches (see below)
  - These batches and the samples contained in them are relevant in regards to the bincov matrices and raw files
  - If you are using the fam_ids input to only run the de novo script on certain samples, you can still include all raw files for this cohort and the script will only localize that files containing your samples of interest

- `batch_bincov_index`: 
  a txt file with first column as batch name, second column as merged_bincov (output of module04-GatherBatchEvidence), and third column as merged_bincov_index (output of module04-GatherBatchEvidence) for each batch
  
  - Must match batch names in sample_batches (see below)
  - These batches and the samples contained in them are relevant in regards to the bincov matrices and raw files
  - If you are using the fam_ids input to only run the de novo script on certain samples, you can still include all bincov matrcies for this cohort and the script will only localize that files containing your samples of interest

- `sample_batches`
  txt file with samples in first column and batch in second column

  - Must match batch names in batch_bincov, batch_raw_file, and batch_depth_raw_file
  - These batches and the samples contained in them are relevant in regards to the bincov matrices and raw files

- `prefix`: choose any prefix which will become the prefix of output files

- `records_per_shard`: number of variants per sharded vcf file; 4000 variants/shard runs in ~1 hr

- `fam_ids`: optional input file of family ids that you want to run the script on
If you are running large cohorts (with one large vcf file) this option can be useful to break the vcf up into x samples per vcf file (and then the de novo calling will only be run on the vcf with those samples). You should try to keep batches together so that you are not localizing the same inputs in multiple runs.





