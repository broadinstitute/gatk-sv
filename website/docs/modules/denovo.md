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


## Identifying de-novo variants

The de-novo workflow identifies a variant as de-novo based the criteria defined 
on site and sample specific metrics of a variant.
These criteria are discussed in details in the following sections. 


#### Site-specific filters

- Exclude multi-allelic copy number variants (mCNV) and breakends (BND).
- Exclude variants with gnomAD `AF > 0.01`, cohort `AF > 0.05` or parents `AF > 0.05` (except for GDs).
- Remove Wham-only calls with `GQ = 1`.

#### Sample-specific filters

- De novo SV (call absent in parents)
- CNVs filters:
  - Large CNVs (>5,000bp) in autosomal chromosomes: 
    - Remove if read depth copy number (RD_CN) of parents is equal to RD_CN of proband.
    - Remove if the variant overlaps with a CNV in the parents, with minimum overlap of 0.5 (using bedtools coverage).
    - Check that variant has raw depth support (uses bedtools coverage for >5000bp, and intersect for 1000-5000bp).
    - Remove if variant overlaps with parental raw depth CNV (uses bedtools coverage for >5000bp, and intersect for 1000-5000bp).
  - Small CNVs (<=5000bp):
    - If the CNV is small, ignore RD evidence if “RD,SR” (and use only SR)
    - Check that variant has raw evidence (bedtools intersect)
    - Remove if variant overlaps with parental raw evidence (bedtools intersect)
    - Remove if variant is SR only but doesn’t have “bothsides support”
  - Remove depth only calls that are <10000bp
  - Remove deletions that are >500bp and RD_CN=2 and PE only evidence
  - Flag if parents GQ is smaller than min_gq (does not apply, we have been using 0)
  - Flag INS that are either manta or melt and SR only and have high SR background (filter removed)
  - Flag if median coverage in parents is <= 10
- Blacklist: remove variants that overlap >0.5 with “blacklist” (repetitive and bad regions)
- INS filters:
  - Requires raw evidence
  - Remove variant if it’s in parental raw evidence file
- No filtering is applied on complex SVs, translocations and inversions
- Reformatting is applied at the end to remove duplicated CPX SVs.


