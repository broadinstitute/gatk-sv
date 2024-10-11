---
title: Resource files
description: Resource files manifest
sidebar_position: 5
---

This page contains descriptions of fixed resource files used in the pipeline. All required files are publicly available 
in Google Cloud Storage Buckets. URIs for these files are available in the [hg38 resources json](https://github.com/broadinstitute/gatk-sv/blob/main/inputs/values/resources_hg38.json).


:::info
Some resources contains sensitive data but are maintained for development purposes. None of these files are required 
to run GATK-SV. These files are located in the bucket `gs://gatk-sv-resources-secure`.
:::

#### allosome_file
Reference [fasta index file](https://www.htslib.org/doc/faidx.html) containing only allosomal contigs.

#### autosome_file
Reference [fasta index file](https://www.htslib.org/doc/faidx.html) containing only autosomal contigs.

#### bin_exclude
[Bed file](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) of intervals to exclude in the call set.

#### cnmops_exclude_list
Bed file of N-masked regions in the reference.

#### contig_ploidy_priors
TSV of prior probabilities contig ploidies used by GATK-gCNV.

#### cytobands
Bed file of cytoband intervals.

#### sd_locs_vcf
SNP sites at which to collect site depth (SD) evidence.

#### depth_exclude_list
Bed file of intervals over which to exclude overlapping depth-only calls.

#### empty_file
Empty file used to satisfy some workflow code paths.

#### exclude_intervals_for_gcnv_filter_intervals
Bed file of intervals to exclude from GATK-gCNV.

#### external_af_ref_bed


#### genome_file

#### manta_region_bed

#### mei_bed

#### melt_std_vcf_header

#### noncoding_bed

#### par_bed

#### pesr_exclude_list

#### preprocessed_intervals

#### primary_contigs_fai

#### primary_contigs_list

#### contigs_header

#### protein_coding_gtf

#### reference_dict

#### reference_fasta

#### reference_index

#### rmsk

#### segdups

#### seed_cutoffs

#### single_sample_qc_definitions

#### wgd_scoring_mask

#### wham_include_list_bed_file

#### aou_recalibrate_gq_model_file

#### hgdp_recalibrate_gq_model_file

#### recalibrate_gq_genome_tracks

#### ccdg_abel_site_level_benchmarking_dataset

#### gnomad_v2_collins_sample_level_benchmarking_dataset

#### gnomad_v2_collins_site_level_benchmarking_dataset

#### hgsv_byrska_bishop_sample_level_benchmarking_dataset

#### hgsv_byrska_bishop_sample_renaming_tsv

#### hgsv_byrska_bishop_site_level_benchmarking_dataset

#### hgsv_ebert_sample_level_benchmarking_dataset

#### ssc_belyeu_sample_level_benchmarking_dataset

#### ssc_belyeu_site_level_benchmarking_dataset

#### ssc_sanders_sample_level_benchmarking_dataset

#### thousand_genomes_site_level_benchmarking_dataset

#### asc_site_level_benchmarking_dataset

#### hgsv_site_level_benchmarking_dataset

#### collins_2017_sample_level_benchmarking_dataset

#### sanders_2015_sample_level_benchmarking_dataset

#### werling_2018_sample_level_benchmarking_dataset

