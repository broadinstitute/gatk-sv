## Copyright Broad Institute, 2022
## 
##
## Consolidate boost scores per sample across all batches and write those scores
## directly into an input VCF
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"


workflow AnnotateVcfWithBoostScores {
  input {
    File vcf
    File vcf_idx
    Array[File] boost_score_tarballs
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_override_subset_vcf
    RuntimeAttr? runtime_override_annotate_vcf
  }

  # Scatter over tarballs of boost scores
  scatter ( boost_res in boost_score_tarballs ) {

    # Subset VCF to samples in tarball
    call SubsetVcf {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        boost_tarball=boost_res,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_vcf
    }

    # Annotate boost scores
    call AnnotateBoostScores {
      input:
        vcf=SubsetVcf.subsetted_vcf,
        vcf_idx=SubsetVcf.subsetted_vcf_idx,
        docker=boost_docker,
        runtime_attr_override=runtime_override_annotate_vcf
    }

  }

  # Column-wise merge of all annotated VCFs
  # TODO: implement this

  output {}
}


task SubsetVcf {
  input {
    File vcf
    File vcf_idx
    File boost_tarball
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String vcf_out_prefix = basename(boost_tarball, ".tar.gz")

  Float input_size = size(select_all([vcf, boost_tarball]), "GB")
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
      mem_gb: 4,
      disk_gb: ceil(base_disk_gb + (input_size * 10.0)),
      cpu_cores: 1,
      preemptible_tries: 3,
      max_retries: 1,
      boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
      memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
      disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
      cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
      preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
      maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
      docker: sv_base_mini_docker
      bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    # Decompress & reorganize tarball
    mkdir boost_results
    tar -xzvf ~{boost_tarball} --directory boost_results/
    find boost_results -name "*.boost_filtered" \
    | xargs -I {} mv {} boost_results/

    # Get list of sample IDs
    find boost_results/ -name "*.boost_filtered" \
    | xargs -I {} basename {} \
    | sed 's/\.boost_filtered/\t/g' | cut -f1 \
    | sort -Vk1,1 | uniq \
    > samples.list

    # Subset & index VCF 
    bcftools view \
      -S samples.list \
      --force-samples \
      -O z \
      -o ~{vcf_out_prefix}.subsetted.vcf.gz \
      ~{vcf}
    tabix -p vcf ~{vcf_out_prefix}.subsetted.vcf.gz
  >>>

  output {
    File subsetted_vcf = "~{vcf_out_prefix}.subsetted.vcf.gz"
    File subsetted_vcf_idx = "~{vcf_out_prefix}.subsetted.vcf.gz.tbi"
  }
}


task AnnotateBoostScores {
  input {
    File vcf
    File vcf_idx
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String vcf_out_prefix = basename(vcf, ".vcf.gz") + "boost_anno"

  Float input_size = size(vcf, "GB")
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
      mem_gb: 4,
      disk_gb: ceil(base_disk_gb + (input_size * 20.0)),
      cpu_cores: 1,
      preemptible_tries: 3,
      max_retries: 1,
      boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
      memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
      disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
      cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
      preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
      maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
      docker: sv_pipeline_docker
      bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    # TODO: implement this

  >>>


}