version 1.0
import "Structs.wdl"

workflow ProcessVCFs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxes
    String sv_pipeline_base_docker
  }

  scatter(i in range(length(vcfs))) {
    call CleanVCF { 
      input: 
        vcf = vcfs[i] ,
        vcf_idx = vcf_idxes[i],
        sv_pipeline_base_docker = sv_pipeline_base_docker
        }
  }

  output {
    Array[File] cleaned_vcfs = CleanVCF.pass_vcf
    Array[File] cleaned_idxes = CleanVCF.pass_idx
  }
}


task CleanVCF {
  input {
    File vcf
    File vcf_idx
    String sv_pipeline_base_docker
  }

  String prefix =  basename(vcf, ".vcf.gz")

  command <<<
    set -euo pipefail

    bcftools view -h ~{vcf} > ~{basename}.PASS_1.vcf
    bcftools view -H ~{vcf} | \
      awk '{ if ($NF=="0|0" || $NF=="0|1" || $NF=="1|0" || $NF=="1|1") print }' \
      >> ~{basename}.PASS_1.vcf
    bgzip -f ~{basename}.PASS_1.vcf
    tabix -p vcf ~{basename}.PASS_1.vcf.gz

    bcftools view -h ~{vcf} > ~{basename}.PASS_2.vcf
    bcftools view -H ~{vcf} | \
      awk '{ if ($NF!="0|0" && $NF!="0|1" && $NF!="1|0" && $NF!="1|1") print }' | \
      cut -f1-9 | \
      sed -e "s/\$/\t1|1/" \
      >> ~{basename}.PASS_2.vcf
    bgzip -f ~{basename}.PASS_2.vcf
    tabix -p vcf ~{basename}.PASS_2.vcf.gz


    bcftools concat \
      -Oz \
      -o ~{basename}.PASS.vcf.gz \
      ~{basename}.PASS_1.vcf.gz \
      ~{basename}.PASS_2.vcf.gz

    tabix -p vcf ~{basename}.PASS.vcf.gz
  >>>

  output {
    File pass_vcf = "~{basename}.PASS.vcf.gz"
    File pass_idx = "~{basename}.PASS.vcf.gz.tbi"
 }

  runtime {
    docker: "biocontainers/bcftools:v1.17-1-deb_cv1"
    cpu: 1
    memory: "4G"
    disks: "local-disk 20 HDD"
  }
}