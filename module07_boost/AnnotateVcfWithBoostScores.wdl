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
    String sv_pipeline_base_docker
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
        boost_tarball=boost_res,
        sv_pipeline_base_docker=sv_pipeline_base_docker,
        runtime_attr_override=runtime_override_annotate_vcf
    }
  }

  # Column-wise merge of all annotated VCFs
  call MergeVcfs {
    input:
      vcfs=AnnotateBoostScores.annotated_vcf,
      vcf_idxs=AnnotateBoostScores.annotated_vcf_idx,
      out_prefix=basename(vcf, ".vcf.gz") + ".boost_annotated",
      sv_base_mini_docker=sv_base_mini_docker
  }

  output {
    File annotated_vcf = MergeVcfs.merged_vcf
    File annotated_vcf_idx = MergeVcfs.merged_vcf_idx
  }
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
    File boost_tarball
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  String out_prefix = basename(vcf, ".vcf.gz") + ".boost_anno"

  Float input_size = size([vcf, boost_tarball], "GB")
  Float base_disk_gb = 100.0
  Float compression_factor = 20.0
  RuntimeAttr runtime_default = object {
      mem_gb: 6,
      disk_gb: ceil(base_disk_gb + (input_size * compression_factor)),
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
      docker: sv_pipeline_base_docker
      bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    # Get list of variant IDs and samples in VCF
    bcftools query -f '%ID\n' ~{vcf} > vids.list
    bcftools query -l ~{vcf} > samples.list

    # Decompress Boost results and subset to variant IDs present in VCF
    mkdir boost_results
    mkdir boost_results/raw
    tar -xzvf ~{boost_tarball} --directory boost_results/raw/
    while read sample; do
      infile=$( find boost_results -name "$sample.boost_filtered" | sed -n '1p' )
      if ! [ -z $infile ]; then
        fgrep -wf vids.list $infile > boost_results/$sample.scores.tsv
        echo -e "$sample\tboost_results/$sample.scores.tsv" >> boost_anno_inputs.tsv
      fi
    done < samples.list

    # Add Boost scores to VCF
    /opt/sv-pipeline/07_filtering/add_boost_scores_to_vcf.py \
      --vcf ~{vcf} \
      --boost-tsv boost_anno_inputs.tsv \
      --outfile ~{out_prefix}.vcf.gz
    tabix -p vcf  ~{out_prefix}.vcf.gz
  >>>

  output {
    File annotated_vcf = "~{out_prefix}.vcf.gz"
    File annotated_vcf_idx = "~{out_prefix}.vcf.gz.tbi"
  }
}


task MergeVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String out_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcfs, "GB")
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
      mem_gb: 4,
      disk_gb: ceil(base_disk_gb + (input_size * 5.0)),
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

    bcftools merge \
      -m id \
      -O z \
      -o ~{out_prefix}.vcf.gz \
      ~{sep=" " vcfs}
    tabix -p vcf ~{out_prefix}.vcf.gz
  >>>

  output {
    File merged_vcf = "~{out_prefix}.vcf.gz"
    File merged_vcf_idx = "~{out_prefix}.vcf.gz.tbi"
  }
}
