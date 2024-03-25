version 1.0

import "Structs.wdl"
import "GenotypeGenomicDisorderRegionsBatch.wdl" as gdr_batch
import "VcfClusterSingleChromsome.wdl" as vcf_cluster
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as util

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

    File? preprocess_intervals_script
    String? revise_args
    File? revise_script
    File? reset_genotypes_script

    String linux_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_attr_subset_by_samples
    RuntimeAttr? runtime_override_concat_batch
    RuntimeAttr? runtime_gdr_overlapping_variants
    RuntimeAttr? runtime_rdtest_full
    RuntimeAttr? runtime_rdtest_subdiv
    RuntimeAttr? runtime_revise_vcf_batch
    RuntimeAttr? runtime_vcf2bed_before_revise
    RuntimeAttr? runtime_vcf2bed_after_revise
    RuntimeAttr? runtime_rdtest_before_revise
    RuntimeAttr? runtime_rdtest_after_revise

    RuntimeAttr? runtime_attr_preprocess
    RuntimeAttr? runtime_cat_subtracted_genotypes
    RuntimeAttr? runtime_get_vcf_header
    RuntimeAttr? runtime_subtract_genotypes
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_reorder_samples
    RuntimeAttr? runtime_set_missing_formats
    RuntimeAttr? runtime_override_concat_revised_vcfs
    RuntimeAttr? runtime_override_concat_new_records
  }

  call PreprocessGenomicDisorderIntervals {
    input:
      prefix = "~{output_prefix}.preprocess_gdr",
      bed = genomic_disorder_regions_bed,
      script = preprocess_intervals_script,
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
        revise_args = revise_args,
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
        runtime_vcf2bed_before_revise = runtime_vcf2bed_before_revise,
        runtime_vcf2bed_after_revise = runtime_vcf2bed_after_revise,
        runtime_rdtest_before_revise = runtime_rdtest_before_revise,
        runtime_rdtest_after_revise = runtime_rdtest_after_revise
    }
  }

  # Note these tsv's are gzipped but can be concatenated normally
  call tasks_cohort.CatUncompressedFiles {
    input:
        shards = GenotypeGenomicDisorderRegionsBatch.batch_gdr_revised_genotypes_tsv,
        outfile_name = "~{output_prefix}.gdr_subtracted_genotypes.tsv.gz",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_cat_subtracted_genotypes
  }
  # Use this make sure all samples are included in the clustering output
  # Needed when running on a subset of batches
  call util.GetVcfHeader {
    input:
      vcf = cohort_vcfs[0],
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_get_vcf_header
  }

  Array[String] contigs = transpose(read_tsv(select_first([contig_list])))[0]
  scatter (i in range(length(cohort_vcfs))) {
    call SubtractGenotypes {
      input:
        prefix="~{output_prefix}.~{contigs[i]}.subtract_genotypes",
        vcf = cohort_vcfs[i],
        vcf_index = cohort_vcfs[i] + ".tbi",
        genotype_tsv = CatUncompressedFiles.outfile,
        contig = contigs[i],
        script = reset_genotypes_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_subtract_genotypes
    }
    call tasks_cluster.SVCluster {
      input:
        vcfs=flatten([[GetVcfHeader.out], GenotypeGenomicDisorderRegionsBatch.batch_gdr_revised_after_update_vcf]),
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
    call SetMissingGenotypingFormatFields {
      input:
        prefix="~{output_prefix}.~{contigs[i]}.set_missing_fields",
        vcf = SVCluster.out,
        vcf_index = SVCluster.out_index,
        script = reset_genotypes_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_set_missing_formats
    }
    call util.ReorderVcfSamples {
      input:
        vcf = SetMissingGenotypingFormatFields.out,
        sample_ordered_vcf = GetVcfHeader.out,
        output_prefix = "~{output_prefix}.~{contigs[i]}.clustered_new_records.reorder_samples",
        generate_index = true,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_reorder_samples
    }
    call tasks_cohort.ConcatVcfs as ConcatVcfsFinal {
      input:
        vcfs = [SubtractGenotypes.out, ReorderVcfSamples.out],
        vcfs_idx = [SubtractGenotypes.out_index, ReorderVcfSamples.out_index],
        allow_overlaps = true,
        outfile_prefix = "~{output_prefix}.concat_revise_gdr",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_concat_revised_vcfs
    }
  }
  call tasks_cohort.ConcatVcfs as ConcatNewRecords {
    input:
      vcfs = ReorderVcfSamples.out,
      vcfs_idx = ReorderVcfSamples.out_index,
      naive = true,
      outfile_prefix = "~{output_prefix}.concat_new_gdr_records",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat_new_records
  }
  output{
    # Cohort VCF outputs
    Array[File] cohort_gdr_revised_vcf = ConcatVcfsFinal.concat_vcf
    Array[File] cohort_gdr_revised_index = ConcatVcfsFinal.concat_vcf_idx

    File cohort_gdr_revised_genotypes = CatUncompressedFiles.outfile
    File cohort_gdr_revised_record_subset_vcf = ConcatNewRecords.concat_vcf
    File cohort_gdr_revised_record_subset_index = ConcatNewRecords.concat_vcf_idx

    # Batch RdTest outputs

    # Plots of input variants that overlap one or more GDRs, with carriers shown
    Array[File] batch_rdtest_variants_overlapping_gdr = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_variants_overlapping_gdr
    # Plots of GDRs that overlap one or more input variants, with carriers shown
    Array[File] batch_rdtest_gdr_overlapping_variants = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_overlapping_variants
    # Plots of all GDRs, carriers not shown (random sample is the carrier)
    Array[File] batch_rdtest_gdr_full = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_full
    # Plots and genotyping of GDR subdivisions (default 10 per region), carriers not shown (random sample is the carrier)
    Array[File] batch_rdtest_gdr_subdiv = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_subdiv

    # Plots of revised input variants with their original genotypes, with carriers shown
    Array[File] batch_rdtest_gdr_before_revise = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_before_revise
    # Plots of revised input variants after revising and new records of rescued false negatives, with carriers shown
    Array[File] batch_rdtest_gdr_after_revise = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_after_revise

    # Batch VCF outputs
    Array[File] batch_gdr_revised_before_update_vcf = GenotypeGenomicDisorderRegionsBatch.batch_gdr_revised_before_update_vcf
    Array[File] batch_gdr_revised_before_update_index = GenotypeGenomicDisorderRegionsBatch.batch_gdr_revised_before_update_index
    Array[File] batch_gdr_revised_after_update_vcf = GenotypeGenomicDisorderRegionsBatch.batch_gdr_revised_after_update_vcf
    Array[File] batch_gdr_revised_after_update_index = GenotypeGenomicDisorderRegionsBatch.batch_gdr_revised_after_update_index
    Array[File] batch_gdr_revised_genotypes_tsv = GenotypeGenomicDisorderRegionsBatch.batch_gdr_revised_genotypes_tsv

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
    File vcf_index
    File genotype_tsv
    String contig
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
    # Filter down to this contig
    zcat ~{genotype_tsv} \
      | awk -F'\t' -v OFS='\t' -v chrom=~{contig} '$1==chrom' \
      | gzip > genotypes.tsv.gz
    N_RECORDS=$(zcat genotypes.tsv.gz | wc -l)
    # If there are no changes on this contig, just copy the input
    if [ "$N_RECORDS" -eq "0" ]; then
      cp ~{vcf} ~{prefix}.vcf.gz
      cp ~{vcf_index} ~{prefix}.vcf.gz.tbi
    else
      python ~{default="/opt/src/sv-pipeline/scripts/reset_del_dup_genotypes.py" script} \
        --vcf ~{vcf} \
        --genotype-tsv ~{genotype_tsv} \
        --out ~{prefix}.vcf.gz
    fi
  >>>
  output {
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


task SetMissingGenotypingFormatFields {
  input{
    String prefix
    File vcf
    File vcf_index
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50.0 + size(vcf, "GiB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    python ~{default="/opt/src/sv-pipeline/scripts/reset_del_dup_genotypes.py" script} \
      --vcf ~{vcf} \
      --reset-rd-genotype \
      --out ~{prefix}.vcf.gz
  >>>
  output {
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