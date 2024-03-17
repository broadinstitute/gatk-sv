version 1.0

import "Structs.wdl"
import "GenotypeGenomicDisorderRegionsBatch.wdl" as gdr_batch

workflow GenotypeGenomicDisorderRegions {
  input {
    String output_prefix
    Array[String] batch_names
    Array[File] rd_files
    Array[File] median_files
    Array[File] depth_sepcutoff_files

    Array[File] cohort_vcfs
    File ped_file
    File genomic_disorder_regions_bed
    File par_bed

    File? revise_script

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_attr_preprocess
    RuntimeAttr? runtime_attr_subset_by_samples
    RuntimeAttr? runtime_override_concat_batch
    RuntimeAttr? runtime_gdr_overlapping_variants
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
        runtime_revise_vcf = runtime_revise_vcf,
        runtime_vcf2bed_new_records = runtime_vcf2bed_new_records,
        runtime_vcf2bed_original_invalid = runtime_vcf2bed_original_invalid,
        runtime_vcf2bed_subracted_invalid = runtime_vcf2bed_subracted_invalid,
        runtime_rdtest_new_records = runtime_rdtest_new_records,
        runtime_rdtest_original_invalid = runtime_rdtest_original_invalid,
        runtime_rdtest_subtracted_invalid = runtime_rdtest_subtracted_invalid
    }
  }
  output{
    Array[File] batch_variants_overlapping_gdr_out = GenotypeGenomicDisorderRegionsBatch.batch_variants_overlapping_gdr_out
    Array[File] batch_gdr_overlapping_variants_out = GenotypeGenomicDisorderRegionsBatch.batch_gdr_overlapping_variants_out
    Array[File] batch_rdtest_gdr_full_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_full_out
    Array[File] batch_rdtest_gdr_subdiv_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_subdiv_out

    Array[File] batch_rdtest_gdr_new_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_new_out
    Array[File] batch_rdtest_gdr_orig_invalid_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_orig_invalid_out
    Array[File] batch_rdtest_gdr_subtracted_invalid_out = GenotypeGenomicDisorderRegionsBatch.batch_rdtest_gdr_subtracted_invalid_out

    Array[File] batch_new_gdr_records_vcf = GenotypeGenomicDisorderRegionsBatch.batch_new_gdr_records_vcf
    Array[File] batch_new_gdr_records_index = GenotypeGenomicDisorderRegionsBatch.batch_new_gdr_records_index
    Array[File] batch_original_invalidated_gdr_records_vcf = GenotypeGenomicDisorderRegionsBatch.batch_original_invalidated_gdr_records_vcf
    Array[File] batch_original_invalidated_gdr_records_index = GenotypeGenomicDisorderRegionsBatch.batch_original_invalidated_gdr_records_index
    Array[File] batch_subtracted_invalidated_gdr_records_vcf = GenotypeGenomicDisorderRegionsBatch.batch_subtracted_invalidated_gdr_records_vcf
    Array[File] batch_subtracted_invalidated_gdr_records_index = GenotypeGenomicDisorderRegionsBatch.batch_subtracted_invalidated_gdr_records_index
    Array[File] batch_gdr_subtracted_vcf = GenotypeGenomicDisorderRegionsBatch.batch_gdr_subtracted_vcf
    Array[File] batch_gdr_subtracted_index = GenotypeGenomicDisorderRegionsBatch.batch_gdr_subtracted_index

    Array[File] batch_subsetted_vcf = GenotypeGenomicDisorderRegionsBatch.batch_subsetted_vcf
  }
}

task RunRdTest {
  input{
    String output_prefix
    File rdtest_bed
    File rd_file
    File rd_index
    File median_file
    File depth_sepcutoff
    Int large_size_cutoff = 1000000
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 15,
                               disk_gb: ceil(40.0 + size(rd_file, "GiB") * 4),
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
      -v TRUE -g TRUE -p TRUE \
      -r ~{depth_sepcutoff} \
      -b intervals.bed \
      -c ~{rd_file} \
      -m ~{median_file} \
      -n ~{output_prefix} \
      -s ~{large_size_cutoff} -o ~{output_prefix}
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
                               mem_gb: 7.5,
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
                               mem_gb: 7.5,
                               disk_gb: ceil(50.0 + size(vcf, "GiB") * 3),
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
    mkdir tmp
    bcftools sort -T tmp/ ~{prefix}.new_revised_records.unsorted.vcf.gz -Oz -o ~{prefix}.new_revised_records.vcf.gz
    tabix ~{prefix}.new_revised_records.vcf.gz
    tabix ~{prefix}.original_revised_records.vcf.gz
    tabix ~{prefix}.subtracted.vcf.gz
  >>>
  output{
    File revised_records_vcf = "~{prefix}.new_revised_records.vcf.gz"
    File revised_records_index = "~{prefix}.new_revised_records.vcf.gz.tbi"
    File original_records_vcf = "~{prefix}.original_revised_records.vcf.gz"
    File original_records_index = "~{prefix}.original_revised_records.vcf.gz.tbi"
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
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bed = "~{prefix}.bed"
  }
  command <<<
    set -euo pipefail
    svtk vcf2bed ~{vcf} ~{prefix}.bed
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