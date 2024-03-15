version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as util

workflow GenotypeGenomicDisorderRegionsBatch {
  input {
    String output_prefix
    String batch_name
    File batch_sample_list
    File rd_file
    File median_file
    File depth_sepcutoff_file

    Array[File] cohort_vcfs

    File ped_file
    File preprocessed_genomic_disorder_regions_bed
    File genomic_disorder_regions_bed
    File par_bed

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_subset_by_samples
    RuntimeAttr? runtime_override_concat_batch
    RuntimeAttr? runtime_rdtest_full
    RuntimeAttr? runtime_rdtest_subdiv
    RuntimeAttr? runtime_revise_vcf
    RuntimeAttr? runtime_vcf2bed_new_records
    RuntimeAttr? runtime_vcf2bed_original_invalid
    RuntimeAttr? runtime_vcf2bed_subracted_invalid
    RuntimeAttr? runtime_rdtest_new_records
    RuntimeAttr? runtime_rdtest_original_invalid
    RuntimeAttr? runtime_rdtest_subtracted_invalid
  }

  scatter (i in range(length(cohort_vcfs))) {
    call util.SubsetVcfBySamplesList {
      input:
        vcf = cohort_vcfs[i],
        list_of_samples = batch_sample_list,
        outfile_name = "~{output_prefix}.shard_{i}",
        remove_samples = false,
        remove_private_sites = true,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_subset_by_samples
    }
  }
  call tasks_cohort.ConcatVcfs {
    input:
      vcfs = SubsetVcfBySamplesList.vcf_subset,
      vcfs_idx = SubsetVcfBySamplesList.vcf_subset_index,
      naive = true,
      outfile_prefix = "~{output_prefix}.concat",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat_batch
  }
  # Run RdTest and generate plots over the full GD regions
  call RunRdTest as RunRdTestFullRegions {
    input:
      output_prefix = "rdtest_full_~{batch_name}",
      rdtest_bed = genomic_disorder_regions_bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_full
  }
  # Run RdTest and generate plots over subdivided GD regions
  call RunRdTest as RunRdTestSubdivision {
    input:
      output_prefix = "rdtest_subdiv_~{batch_name}",
      rdtest_bed = preprocessed_genomic_disorder_regions_bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      depth_sepcutoff = depth_sepcutoff_file,
      do_plot = true,
      do_genotyping = true,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_subdiv
  }
  call ReviseGenomicDisorderRegions {
    input:
      prefix = "~{output_prefix}.revise_gdr",
      rdtest_tars = [RunRdTestSubdivision.out],
      vcf = ConcatVcfs.concat_vcf,
      ped_file = ped_file,
      genomic_disorder_regions_bed = genomic_disorder_regions_bed,
      par_bed = par_bed,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_revise_vcf
  }
  call VcfToBed as VcfToBedNewRecords {
    input:
      vcf = ReviseGenomicDisorderRegions.new_records_vcf,
      prefix = "~{output_prefix}.new_records",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_vcf2bed_new_records
  }
  call VcfToBed as VcfToBedOriginalInvalidatedRecords {
    input:
      vcf = ReviseGenomicDisorderRegions.original_invalidated_records_vcf,
      prefix = "~{output_prefix}.original_invalidated_records",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_vcf2bed_original_invalid
  }
  call VcfToBed as VcfToBedSubtractedInvalidatedRecords {
    input:
      vcf = ReviseGenomicDisorderRegions.subtracted_invalidated_records_vcf,
      prefix = "~{output_prefix}.subtracted_invalidated_records",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_vcf2bed_subracted_invalid
  }
  # Run RdTest to visualize results
  call RunRdTest as RunRdTestNewRecords {
    input:
      output_prefix = "rdtest_new_~{batch_name}",
      rdtest_bed = VcfToBedNewRecords.bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_new_records
  }
  call RunRdTest as RunRdTestOriginalInvalid {
    input:
      output_prefix = "rdtest_orig_invalid_~{batch_name}",
      rdtest_bed = VcfToBedOriginalInvalidatedRecords.bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_original_invalid
  }
  call RunRdTest as RunRdTestSubtractedInvalid {
    input:
      output_prefix = "rdtest_subtract_invalid_~{batch_name}",
      rdtest_bed = VcfToBedOriginalInvalidatedRecords.bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_subtracted_invalid
  }
  output{
    File batch_rdtest_gdr_full_out = RunRdTestFullRegions.out
    File batch_rdtest_gdr_subdiv_out = RunRdTestFullRegions.out

    File batch_rdtest_gdr_new_out = RunRdTestNewRecords.out
    File batch_rdtest_gdr_orig_invalid_out = RunRdTestOriginalInvalid.out
    File batch_rdtest_gdr_subtracted_invalid_out = RunRdTestSubtractedInvalid.out

    File batch_new_gdr_records_vcf = ReviseGenomicDisorderRegions.new_records_vcf
    File batch_new_gdr_records_index = ReviseGenomicDisorderRegions.new_records_index

    File batch_original_invalidated_gdr_records_vcf = ReviseGenomicDisorderRegions.original_invalidated_records_index
    File batch_original_invalidated_gdr_records_index = ReviseGenomicDisorderRegions.original_invalidated_records_index

    File batch_subtracted_invalidated_gdr_records_vcf = ReviseGenomicDisorderRegions.subtracted_invalidated_records_vcf
    File batch_subtracted_invalidated_gdr_records_index = ReviseGenomicDisorderRegions.subtracted_invalidated_records_index

    File batch_gdr_subtracted_vcf = ReviseGenomicDisorderRegions.subtracted_vcf
    File batch_gdr_subtracted_index = ReviseGenomicDisorderRegions.subtracted_index
  }
}

task RunRdTest {
  input{
    String output_prefix
    File rdtest_bed
    File rd_file
    File rd_index
    File median_file
    File? depth_sepcutoff  # Required if do_genotyping = true
    Boolean do_plot
    Boolean do_genotyping
    Int large_size_cutoff = 1000000
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(40.0 + size(rd_file, "GiB") * 4 + size(rdtest_bed, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    # Inject one sample from the batch into the 5th column
    SAMPLE=$(awk -F'\t' '{ if (NR==1) {print $1} }' ~{median_file})
    awk -F'\t' -v OFS='\t' -v s="$SAMPLE" '{print $1,$2,$3,$4,s,$5}' ~{rdtest_bed} > intervals.bed
    mkdir ~{output_prefix}/
    Rscript /opt/RdTest/RdTest.R \
      ~{if do_genotyping then "-g TRUE -v TRUE" else ""}
      ~{if do_plot then "-p TRUE" else ""} \
      ~{"-r " + depth_sepcutoff} \
      -b intervals.bed \
      -c ~{rd_file} \
      -m ~{median_file} \
      -n ~{output_prefix} \
      -s ~{large_size_cutoff} \
      -o ~{output_prefix}
    tar czvf ~{output_prefix}.tar.gz ~{output_prefix}/
  >>>
  output{
    File out = "~{output_prefix}.tar.gz"
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

task ReviseGenomicDisorderRegions {
  input{
    String prefix
    Array[File] rdtest_tars
    File vcf
    File ped_file
    File genomic_disorder_regions_bed
    File par_bed
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50.0 + size(vcf, "GiB") * 3 + size(rdtest_tars, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    mkdir rdtest/
    while read -r FILE; do
      tar xzf $FILE -C rdtest/
    done < ~{write_lines(rdtest_tars)}
    ls rdtest/*/*.median_geno > median_geno_files.list
    python ~{default="/opt/src/sv-pipeline/scripts/revise_genomic_disorder_regions.py" script} \
      --vcf ~{vcf} \
      --median-geno-list median_geno_files.list\
      --ped-file ~{ped_file} \
      --region-bed ~{genomic_disorder_regions_bed} \
      --par-bed ~{par_bed} \
      --out ~{prefix}
  >>>
  output{
    File new_records_vcf = "~{prefix}.new_records.vcf.gz"
    File new_records_index = "~{prefix}.new_records.vcf.gz.tbi"
    File original_invalidated_records_vcf = "~{prefix}.original_invalidated_records.vcf.gz"
    File original_invalidated_records_index = "~{prefix}.original_invalidated_records.vcf.gz.tbi"
    File subtracted_invalidated_records_vcf = "~{prefix}.subtracted_invalidated_records.vcf.gz"
    File subtracted_invalidated_records_index = "~{prefix}.subtracted_invalidated_records.vcf.gz.tbi"
    File subtracted_vcf = "~{prefix}.subtracted.vcf.gz"
    File subtracted_index = "~{prefix}.subtracted.vcf.gz.tbi"
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


task VcfToBed {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 50 + size(vcf, "GiB") * 2,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bed = "~{prefix}.bed"
  }
  command <<<
    set -euxo pipefail
    svtk vcf2bed ~{vcf} intervals.bed
    # Swap last two columns for RdTest
    awk -F"\t" -v OFS="\t" '{print $1, $2, $3, $4, $6, $5} intervals.bed > ~{prefix}.bed
  >>>
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