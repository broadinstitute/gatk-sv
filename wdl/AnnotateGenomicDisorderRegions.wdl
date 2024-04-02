version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as util

workflow AnnotateGenomicDisorderRegions {
  input {
    String output_prefix
    Array[File] vcfs
    File genomic_disorder_regions_bed

    Float? overlap
    String? annotate_additional_args

    File? annotate_gdr_script

    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_annotate
  }

  scatter (i in range(length(vcfs))) {
    call AnnotateGenomicDisorderRegionsTask {
      input:
        prefix="~{output_prefix}.annotate_gdr",
        vcf=vcfs[i],
        region_bed=genomic_disorder_regions_bed,
        overlap=overlap,
        additional_args=annotate_additional_args,
        script=annotate_gdr_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_annotate
    }
  }

  output{
    # Cohort VCF outputs
    Array[File] gdr_annotated_vcf = AnnotateGenomicDisorderRegionsTask.out_vcf
    Array[File] gdr_annotated_vcf_index = AnnotateGenomicDisorderRegionsTask.out_vcf_index
  }
}

task AnnotateGenomicDisorderRegionsTask {
  input{
    String prefix
    String vcf
    File region_bed
    Float? overlap
    String? additional_args
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
    python ~{default="/opt/src/sv-pipeline/scripts/annotate_genomic_disorder_regions.py" script} \
      --region-bed ~{region_bed} \
      --vcf ~{vcf} --out annotate_test
      --input ~{region_bed} \
      --out ~{prefix} \
      ~{"--overlap " + overlap} \
      ~{additional_args}
  >>>
  output{
    File out_vcf = "~{prefix}.vcf.gz"
    File out_vcf_index = "~{prefix}.vcf.gz.tbi"
    File out_variant_rdtest_bed = "~{prefix}.padded_variants.rdtest.bed"
    File out_region_rdtest_bed = "~{prefix}.padded_regions.rdtest.bed"
    File out_manifest_tsv = "~{prefix}.manifest.tsv"
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
