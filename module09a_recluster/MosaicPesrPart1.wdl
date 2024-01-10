version 1.0

import "Structs.wdl"

workflow Mosaic{
  input{
    String name
    Array[File] pesr_vcfs
    File metrics
    File cutoffs
    File coverage_file
    File coverage_file_idx
    File fam_file
    File median_file
    String sv_pipeline_docker
  }
  call MergePesrVcfs {
    input:
      pesr_vcfs=pesr_vcfs,
      batch=name,
      sv_pipeline_docker=sv_pipeline_docker
  }
  call GetPotential{
    input:
      name=name,
      metrics=metrics,
      cutoffs=cutoffs,
      depth_vcf=MergePesrVcfs.merged_pesr_vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }
  output{
    File merged_pesr=MergePesrVcfs.merged_pesr_vcf
    File common_potential=GetPotential.common
  }

}
task MergePesrVcfs {
  input{
    Array[File] pesr_vcfs
    String batch
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    for VCF in ~{sep=" " pesr_vcfs}; do
      bcftools view --min-ac 1 $VCF |bgzip -c > temp.vcf.gz
      mv temp.vcf.gz $VCF
    done

    vcf-concat ~{sep=" " pesr_vcfs} \
      | vcf-sort -c \
      | bgzip -c > \
      ~{batch}.filtered_pesr_merged.vcf.gz
  >>>

  output {
    File merged_pesr_vcf = "~{batch}.filtered_pesr_merged.vcf.gz"
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

task GetPotential{
  input{
    String name
    File metrics
    File cutoffs
    File depth_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command<<<
    set -euo pipefail
    cut -f 1,2,3,4,36,37,9,15,27,4 ~{metrics} > ~{name}.rd.metrics
    awk '{if ($1!~"depth" && $4>5000) print}' ~{name}.rd.metrics|egrep "DUP|DEL"  > ~{name}.depth.metrics
    delmed=$(fgrep PESR ~{cutoffs}|fgrep Median|awk '{if($3=="RD_Median_Separation") print $2}')
    dupmed=$delmed
    delp=$(fgrep PESR ~{cutoffs}|fgrep RD_log_pval|awk '{if($3=="RD_log_pval") print $2}')
    dupp=$delp
    pe_p=$(fgrep PE_log_pval ~{cutoffs}|cut -f 2)
    sr_p=$(fgrep SR_sum_log_pval ~{cutoffs}|cut -f 2)
    pesr_p=$(fgrep PESR_log_pval ~{cutoffs}|cut -f 2)
    awk -v delp="$delp" -v delmed="$delmed" -v pe_p="$pe_p" -v sr_p="$sr_p" -v pesr_p="$pesr_p" '{if ($3=="DEL" && $8<delmed && $9>delp ) print}' ~{name}.depth.metrics > del.potentialmosaic.txt 
    awk -v dupp="$dupp" -v dupmed="$dupmed" -v pe_p="$pe_p" -v sr_p="$sr_p" -v pesr_p="$pesr_p" '{if ($3=="DUP" && $8<dupmed && $9>dupp ) print}' ~{name}.depth.metrics> dup.potentialmosaic.txt
    cat del.potentialmosaic.txt dup.potentialmosaic.txt |cut -f1 > potentialmosaic.txt
    tabix -f ~{depth_vcf}
    tabix -H ~{depth_vcf} > head.txt
    zcat ~{depth_vcf} |fgrep -w -f potentialmosaic.txt >body.txt
    cat head.txt body.txt |bgzip -c > test.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_calls.sh -x 1 test.vcf.gz test1.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_calls.sh -x 1 test1.vcf.gz test2.vcf.gz
    bash /opt/sv-pipeline/04_variant_resolution/scripts/stitch_fragmented_calls.sh -x 1 test2.vcf.gz test3.vcf.gz
    svtk vcf2bed test3.vcf.gz ~{name}.potentialmosaic.bed

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
    output{
        File common="~{name}.potentialmosaic.bed"
    }
}
