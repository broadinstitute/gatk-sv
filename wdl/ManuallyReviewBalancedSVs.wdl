version 1.0

import "Structs.wdl"
import "ManuallyReviewBalancedSVsPerBatch.wdl" as batch_rev

workflow ManuallyReviewBalancedSVs {
  input {
    String prefix

    Array[File] cohort_vcfs  # cohort vcf or vcfs sharded by contig
    Array[File] batch_manta_tloc_vcfs
    Array[File] batch_pe_files

    Array[String] batches
    Array[File] samples_in_batches

    Int min_size

    File generate_pe_tabix_py_script # for development

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_select_ctx
    RuntimeAttr? runtime_attr_select_cpx
    RuntimeAttr? runtime_attr_select_inv
    RuntimeAttr? runtime_attr_subset_samples
    RuntimeAttr? runtime_attr_combine_tlocs
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_generate_script
    RuntimeAttr? runtime_attr_collect_pe
    RuntimeAttr? runtime_attr_concat_ctx
    RuntimeAttr? runtime_attr_concat_cpx
    RuntimeAttr? runtime_attr_concat_inv

  }

  # select svs cohort-wide and merge across contigs
  # then select samples in batch in separate step
  call SelectSVType as SelectCTX {
    input:
      vcfs = cohort_vcfs,
      svtype = "CTX",
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_select_ctx
  }
  call SelectSVType as SelectCPX {
    input:
      vcfs = cohort_vcfs,
      svtype = "CPX",
      min_size=min_size,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_select_cpx
  }
  call SelectSVType as SelectINV {
    input:
      vcfs = cohort_vcfs,
      svtype = "INV",
      min_size=min_size,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_select_inv
  }

  scatter (i in range(length(batches))) {
    call batch_rev.ManuallyReviewBalancedSVsPerBatch as ManuallyReviewCTXPerBatch {
      input:
        batch = batches[i],
        svtype = "CTX",
        cohort_vcf = SelectCTX.svtype_vcf,
        cohort_vcf_index = SelectCTX.svtype_vcf_index,
        batch_pe_file = batch_pe_files[i],
        batch_manta_tloc_vcf = batch_manta_tloc_vcfs[i],
        batch_samples = samples_in_batches[i],
        generate_pe_tabix_py_script=generate_pe_tabix_py_script,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_subset_samples=runtime_attr_subset_samples,
        runtime_attr_combine_tlocs=runtime_attr_combine_tlocs,
        runtime_attr_vcf2bed=runtime_attr_vcf2bed,
        runtime_attr_generate_script=runtime_attr_generate_script,
        runtime_attr_collect_pe=runtime_attr_collect_pe
    }
    call batch_rev.ManuallyReviewBalancedSVsPerBatch as ManuallyReviewCPXPerBatch {
      input:
        batch = batches[i],
        svtype = "CPX",
        cohort_vcf = SelectCPX.svtype_vcf,
        cohort_vcf_index = SelectCPX.svtype_vcf_index,
        batch_pe_file = batch_pe_files[i],
        batch_samples = samples_in_batches[i],
        generate_pe_tabix_py_script=generate_pe_tabix_py_script,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_subset_samples=runtime_attr_subset_samples,
        runtime_attr_combine_tlocs=runtime_attr_combine_tlocs,
        runtime_attr_vcf2bed=runtime_attr_vcf2bed,
        runtime_attr_generate_script=runtime_attr_generate_script,
        runtime_attr_collect_pe=runtime_attr_collect_pe
    }
    call batch_rev.ManuallyReviewBalancedSVsPerBatch as ManuallyReviewINVPerBatch {
      input:
        batch = batches[i],
        svtype = "INV",
        cohort_vcf = SelectINV.svtype_vcf,
        cohort_vcf_index = SelectINV.svtype_vcf_index,
        batch_pe_file = batch_pe_files[i],
        batch_samples = samples_in_batches[i],
        generate_pe_tabix_py_script=generate_pe_tabix_py_script,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_subset_samples=runtime_attr_subset_samples,
        runtime_attr_combine_tlocs=runtime_attr_combine_tlocs,
        runtime_attr_vcf2bed=runtime_attr_vcf2bed,
        runtime_attr_generate_script=runtime_attr_generate_script,
        runtime_attr_collect_pe=runtime_attr_collect_pe
    }
  }

  call ConcatEvidences as ConcatCTXEvidences {
    input:
      prefix = "~{prefix}.CTX",
      evidences = ManuallyReviewCTXPerBatch.batch_pe_evidence,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_ctx
  }

  call ConcatEvidences as ConcatCPXEvidences {
    input:
      prefix = "~{prefix}.CPX",
      evidences = ManuallyReviewCPXPerBatch.batch_pe_evidence,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_cpx
  }

  call ConcatEvidences as ConcatINVEvidences {
    input:
      prefix = "~{prefix}.INV",
      evidences = ManuallyReviewINVPerBatch.batch_pe_evidence,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_inv
  }

  output {
    File ctx_evidence = ConcatCTXEvidences.concat_evidence
    File cpx_evidence = ConcatCPXEvidences.concat_evidence
    File inv_evidence = ConcatINVEvidences.concat_evidence
  }
}


task SelectSVType {
  input {
    String prefix
    Array[File] vcfs
    String svtype
    Int? min_size
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String size_filter_cmd = if defined(min_size) then " && INFO/SVLEN>=~{min_size}" else ""

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 3 * size(vcfs, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    while read vcf; do
      bcftools view \
        -i "INFO/SVTYPE=='~{svtype}'~{size_filter_cmd}" \
        $vcf \
        | bcftools +fill-tags \
        -O z \
        -o subset."$(basename $file)" \
        -- -t AC,AN,AF
      echo subset."$(basename $file)" >> for_concat.txt
    done < ~{write_lines(vcfs)}

    bcftools concat --no-version --naive -O z -o ~{prefix}.~{svtype}.vcf.gz --file-list for_concat.txt
    tabix ~{prefix}.~{svtype}.vcf.gz
  >>>

  output {
    File svtype_vcf = "~{prefix}.~{svtype}.vcf.gz"
    File svtype_vcf_index = "~{prefix}.~{svtype}.vcf.gz.tbi"
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


task ConcatEvidences{
  input {
    Array[File] evidences
    String prefix
    String sv_base_mini_docker
   RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    # get only one header
    zcat ~{evidences[0]} | awk 'NR==1' | bgzip -c > ~{prefix}.PE_review.txt.gz
    while read FILE; do
      zcat $FILE | sed '1d'
    done < ~{write_lines(evidences)} \
      | bgzip -c \
      >> ~{prefix}.PE_review.txt.gz

  >>>

  output {
    File concat_evidence = "~{prefix}.PE_review.txt.gz"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

