version 1.0

import "Structs.wdl"
import "GenotypeGenomicDisorderRegionsBatch.wdl" as gdr_batch
import "VcfClusterSingleChromsome.wdl" as vcf_cluster
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow GenotypeGenomicDisorderRegions {
  input {
    String output_prefix
    String variant_prefix
    Array[String] batch_names
    Array[File] rd_files
    Array[File] median_files
    Array[File] depth_sepcutoff_files

    Array[File] cohort_vcfs
    File contig_list

    File ploidy_table
    File ped_file

    Float? min_gdr_overlap_frac_plotting

    File genomic_disorder_regions_bed
    File par_bed
    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    # Required only if running a subset of batches. Must include all samples in the cohort.
    File? vcf_header

    File? revise_script

    String linux_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_attr_preprocess
    RuntimeAttr? runtime_attr_subset_by_samples
    RuntimeAttr? runtime_override_concat_batch
    RuntimeAttr? runtime_gdr_overlapping_variants
    RuntimeAttr? runtime_rdtest_full
    RuntimeAttr? runtime_rdtest_subdiv
    RuntimeAttr? runtime_revise_vcf_batch
    RuntimeAttr? runtime_vcf2bed_new_records
    RuntimeAttr? runtime_vcf2bed_original_invalid
    RuntimeAttr? runtime_vcf2bed_subracted_invalid
    RuntimeAttr? runtime_rdtest_new_records
    RuntimeAttr? runtime_rdtest_original_invalid
    RuntimeAttr? runtime_rdtest_subtracted_invalid

    RuntimeAttr? runtime_cat_subtracted_genotypes
    RuntimeAttr? runtime_subtract_genotypes
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_revise_vcf_cohort
    RuntimeAttr? runtime_override_concat_revised_vcfs
    RuntimeAttr? runtime_override_concat_new_records
  }

  call PreprocessGenomicDisorderIntervals {
    input:
      prefix = "~{output_prefix}.preprocess_gdr",
      bed = genomic_disorder_regions_bed,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_preprocess
  }

  scatter (i in range(length(batch_names))) {
    call gdr_batch.GenotypeGenomicDisorderRegionsBatch {
      input:
        output_prefix = "~{output_prefix}.~{batch_names[i]}",
        batch_name = batch_names[i],
        rd_file = rd_files[i],
        median_file = median_files[i],
        depth_sepcutoff_file = depth_sepcutoff_files[i],
        cohort_vcfs = cohort_vcfs,
        ped_file = ped_file,
        preprocessed_genomic_disorder_regions_bed = PreprocessGenomicDisorderIntervals.out,
        genomic_disorder_regions_bed = genomic_disorder_regions_bed,
        par_bed = par_bed,
        min_gdr_overlap_frac_plotting = min_gdr_overlap_frac_plotting,
        revise_script = revise_script,
        linux_docker = linux_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_override_ids_from_median = runtime_override_ids_from_median,
        runtime_attr_subset_by_samples = runtime_attr_subset_by_samples,
        runtime_override_concat_batch = runtime_override_concat_batch,
        runtime_gdr_overlapping_variants = runtime_gdr_overlapping_variants,
        runtime_rdtest_full = runtime_rdtest_full,
        runtime_rdtest_subdiv = runtime_rdtest_subdiv,
        runtime_revise_vcf_batch = runtime_revise_vcf_batch,
        runtime_vcf2bed_new_records = runtime_vcf2bed_new_records,
        runtime_vcf2bed_original_invalid = runtime_vcf2bed_original_invalid,
        runtime_vcf2bed_subracted_invalid = runtime_vcf2bed_subracted_invalid,
        runtime_rdtest_new_records = runtime_rdtest_new_records,
        runtime_rdtest_original_invalid = runtime_rdtest_original_invalid,
        runtime_rdtest_subtracted_invalid = runtime_rdtest_subtracted_invalid
    }
  }

  # Note these tsv's are gzipped but can be concatenated normally
  call tasks_cohort.CatUncompressedFiles {
    input:
        shards = GenotypeGenomicDisorderRegionsBatch.batch_gdr_genotypes_tsv,
        outfile_name = "~{output_prefix}.gdr_subtracted_genotypes",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_cat_subtracted_genotypes
  }

  Array[String] contigs = transpose(read_tsv(select_first([contig_list])))[0]
  scatter (i in range(length(cohort_vcfs))) {
    call SubtractGenotypes {
      input:
        prefix="~{output_prefix}.~{contigs[i]}.subtract_genotypes",
        vcf = cohort_vcfs[i],
        genotype_tsv = CatUncompressedFiles.outfile,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_subtract_genotypes
    }
    call tasks_cluster.SVCluster {
      input:
        vcfs=select_all(flatten([[vcf_header], GenotypeGenomicDisorderRegionsBatch.batch_new_gdr_records_vcf])),
        ploidy_table=ploidy_table,
        output_prefix="~{output_prefix}.~{contigs[i]}.clustered_new_records",
        contig=contigs[i],
        fast_mode=false,
        pesr_sample_overlap=0,
        pesr_interval_overlap=1,
        pesr_breakend_window=0,
        depth_sample_overlap=0,
        depth_interval_overlap=1,
        depth_breakend_window=0,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        variant_prefix="~{variant_prefix}_~{contigs[i]}_GDR_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
    call tasks_cohort.ConcatVcfs as ConcatVcfsFinal {
      input:
        vcfs = [SubtractGenotypes.out, SVCluster.out],
        vcfs_idx = [SubtractGenotypes.out_index, SVCluster.out_index],
        allow_overlaps = true,
        outfile_prefix = "~{output_prefix}.concat_revise_gdr",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_concat_revised_vcfs
    }
  }
  call tasks_cohort.ConcatVcfs as ConcatNewRecords {
    input:
      vcfs = SVCluster.out,
      vcfs_idx = SVCluster.out_index,
      naive = true,
      outfile_prefix = "~{output_prefix}.concat_new_gdr_records",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat_new_records
  }
  output{
    # Cohort VCF outputs
    Array[File] cohort_gdr_revised_vcf = ConcatVcfsFinal.concat_vcf
    Array[File] cohort_gdr_revised_vcf_index = ConcatVcfsFinal.concat_vcf_idx

    File cohort_gdr_subtracted_genotypes = CatUncompressedFiles.outfile
    File cohort_gdr_new_records_vcf = ConcatNewRecords.concat_vcf
    File cohort_gdr_new_records_vcf_index = ConcatNewRecords.concat_vcf_idx

    # Batch RdTest outputs

    # Plots of input variants that overlap one or more GDRs, with carriers shown
    Array[File] batch_variants_overlapping_gdr_out = GenotypeGenomicDisorderRegionsBatch.batch_variants_overlapping_gdr_out
    # Plots of GDRs that overlap one or more input variants, with carriers shown
    Array[File] batch_gdr_overlapping_variants_out = GenotypeGenomicDisorderRegionsBatch.batch_gdr_overlapping_variants_out
    # Plots of all GDRs, carriers not shown (random sample is the carrier)
    Array[File] batch_rdtest_gdr_full_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_full_out
    # Plots and genotyping of GDR subdivisions (default 10 per region), carriers not shown (random sample is the carrier)
    Array[File] batch_rdtest_gdr_subdiv_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_subdiv_out

    # Plots of new variants with novel breakpoints generated by variant revision, with carriers shown
    Array[File] batch_rdtest_gdr_new_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_new_out
    # Plots of revised input variants with their original genotypes, with carriers shown
    Array[File] batch_rdtest_gdr_orig_invalid_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_orig_invalid_out
    # Plots of revised input variants after subtracting filtered genotypes, with carriers shown
    Array[File] batch_rdtest_gdr_subtracted_invalid_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_subtracted_invalid_out

    # Batch VCF outputs

    Array[File] batch_new_gdr_records_vcf = GenotypeGenomicDisorderRegionsBatch.batch_new_gdr_records_vcf
    Array[File] batch_new_gdr_records_index = GenotypeGenomicDisorderRegionsBatch.batch_new_gdr_records_index
    Array[File] batch_original_invalidated_gdr_records_vcf = GenotypeGenomicDisorderRegionsBatch.batch_original_invalidated_gdr_records_vcf
    Array[File] batch_original_invalidated_gdr_records_index = GenotypeGenomicDisorderRegionsBatch.batch_original_invalidated_gdr_records_index
    Array[File] batch_subtracted_invalidated_gdr_records_vcf = GenotypeGenomicDisorderRegionsBatch.batch_subtracted_invalidated_gdr_records_vcf
    Array[File] batch_subtracted_invalidated_gdr_records_index = GenotypeGenomicDisorderRegionsBatch.batch_subtracted_invalidated_gdr_records_index
    Array[File] batch_gdr_subtracted_genotypes_tsv = GenotypeGenomicDisorderRegionsBatch.batch_gdr_genotypes_tsv

    Array[File] batch_gdr_subsetted_vcf = GenotypeGenomicDisorderRegionsBatch.batch_subsetted_vcf
    Array[File] batch_gdr_subsetted_index = GenotypeGenomicDisorderRegionsBatch.batch_subsetted_index
  }
}

task PreprocessGenomicDisorderIntervals {
  input{
    String prefix
    File bed
    Int? subdivisions
    Int? min_size
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    python ~{default="/opt/src/sv-pipeline/scripts/preprocess_genomic_disorder_regions.py" script} \
      --input ~{bed} \
      --out ~{prefix}.bed \
      ~{"--sudivisions " + subdivisions} \
      ~{"--min-size " + min_size}
  >>>
  output{
    File out = "~{prefix}.bed"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task SubtractGenotypes {
  input{
    String prefix
    File vcf
    File genotype_tsv
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50.0 + size(vcf, "GiB") * 2 + size(genotype_tsv, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    python ~{default="/opt/src/sv-pipeline/scripts/reset_del_dup_genotypes.py" script} \
      --vcf ~{vcf} \
      --genotype-tsv ~{genotype_tsv} \
      --out ~{prefix}.vcf.gz
  >>>
  output{
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
