version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "MainVcfQc.wdl" as qc

workflow FilterGenotypes {

  input {
    File vcf
    String? output_prefix
    File ploidy_table
    File sl_cutoff_table
    File? optimized_sl_cutoff_table
    Float no_call_rate_cutoff = 0.05 # Set to 1 to disable NCR filtering
    String? sl_filter_args # Explicitly set SL arguments - see apply_sl_filter.py

    Int filter_vcf_records_per_shard = 20000
    String header_drop_fields = "FILTER/LOW_QUALITY,FORMAT/TRUTH_CN_EQUAL,FORMAT/GT_FILTER,FORMAT/CONC_ST,INFO/STATUS,INFO/TRUTH_AC,INFO/TRUTH_AN,INFO/TRUTH_AF,INFO/TRUTH_VID,INFO/CNV_CONCORDANCE,INFO/GENOTYPE_CONCORDANCE,INFO/HET_PPV,INFO/HET_SENSITIVITY,INFO/HOMVAR_PPV,INFO/HOMVAR_SENSITIVITY,INFO/MINSL,INFO/NON_REF_GENOTYPE_CONCORDANCE,INFO/SL_MAX,INFO/SL_MEAN,INFO/VAR_PPV,INFO/VAR_SENSITIVITY,INFO/VAR_SPECIFICITY"

    # For MainVcfQc
    File primary_contigs_fai
    File? ped_file
    Array[Array[String]]? site_level_comparison_datasets    # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to directory with one BED per population]
    Array[Array[String]]? sample_level_comparison_datasets  # Array of two-element arrays, one per dataset, each of format [prefix, gs:// path to per-sample tarballs]
    File? sample_renaming_tsv # File with mapping to rename per-sample benchmark sample IDs for compatibility with cohort
    Boolean run_qc = true
    String qc_bcftools_preprocessing_options = "-i 'FILTER=\"PASS\" || FILTER=\"MULTIALLELIC\"'"
    Int qc_sv_per_shard = 2500
    Int qc_samples_per_shard = 600
    RuntimeAttr? runtime_override_plot_qc_per_family
    RuntimeAttr? runtime_override_per_sample_benchmark_plot
    
    String sv_base_mini_docker
    String sv_pipeline_docker
  }

  String output_prefix_ = if defined(output_prefix) then select_first([output_prefix]) else basename(vcf, ".vcf.gz")

  call tasks_cohort.ScatterVcf as ScatterForFilter {
    input:
      vcf=vcf,
      records_per_shard=filter_vcf_records_per_shard,
      prefix="~{output_prefix_}.filter_genotypes_scatter",
      sv_pipeline_docker=sv_pipeline_docker
  }

  scatter ( i in range(length(ScatterForFilter.shards)) ) {
    call FilterVcf {
      input:
        vcf=ScatterForFilter.shards[i],
        ploidy_table=ploidy_table,
        ncr_threshold=no_call_rate_cutoff,
        sl_cutoff_table=select_first([optimized_sl_cutoff_table, sl_cutoff_table]),
        args=sl_filter_args,
        output_prefix="~{output_prefix_}.filter_genotypes.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=FilterVcf.out,
      naive=true,
      outfile_prefix="~{output_prefix_}.filter_genotypes",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call SanitizeHeader {
    input:
      vcf = ConcatVcfs.concat_vcf,
      vcf_index = ConcatVcfs.concat_vcf_idx,
      drop_fields = header_drop_fields,
      prefix = "~{output_prefix_}.filter_genotypes.sanitized",
      sv_pipeline_docker = sv_pipeline_docker
  }

  if (run_qc) {
    call qc.MainVcfQc {
      input:
        vcfs=[SanitizeHeader.out],
        prefix="~{output_prefix_}.filter_genotypes",
        ped_file=ped_file,
        bcftools_preprocessing_options=qc_bcftools_preprocessing_options,
        primary_contigs_fai=primary_contigs_fai,
        site_level_comparison_datasets=site_level_comparison_datasets,
        sample_level_comparison_datasets=sample_level_comparison_datasets,
        sample_renaming_tsv=sample_renaming_tsv,
        sv_per_shard=qc_sv_per_shard,
        runtime_override_per_sample_benchmark_plot=runtime_override_per_sample_benchmark_plot,
        runtime_override_plot_qc_per_family=runtime_override_plot_qc_per_family,
        samples_per_shard=qc_samples_per_shard,
        sv_pipeline_qc_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  output {
    File filtered_vcf = SanitizeHeader.out
    File filtered_vcf_index = SanitizeHeader.out_index
    File? main_vcf_qc_tarball = MainVcfQc.sv_vcf_qc_output
  }
}

task FilterVcf {
  input {
    File vcf
    File ploidy_table
    Float ncr_threshold
    File? sl_cutoff_table
    File? script
    String? args
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(100 + size(vcf, "GB") * 2 + size(sl_cutoff_table, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Convert format
    python ~{default="/opt/sv-pipeline/scripts/apply_sl_filter.py" script} \
      --vcf ~{vcf} \
      --out tmp.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      --ncr-threshold ~{ncr_threshold} \
      ~{if defined(sl_cutoff_table) then "--sl-cutoff-table " + sl_cutoff_table else ""} \
      ~{args}

    bcftools +fill-tags tmp.vcf.gz -- -t AC,AN,AF \
      | bcftools view --no-update -i 'SVTYPE=="CNV" || AC>0' -Oz -o ~{output_prefix}.vcf.gz
    tabix ~{output_prefix}.vcf.gz
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

task SanitizeHeader {
  input {
    File vcf
    File vcf_index
    String prefix
    String drop_fields
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + size(vcf, "GiB") * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    bcftools view --no-version -h ~{vcf} > header.vcf

    grep -v -e ^"##bcftools" header.vcf \
      -e ^"##source=depth" \
      -e ^"##source=cleanvcf" \
      -e ^"##ALT=<ID=UNR" \
      | sed 's/Split read genotype quality/Split-read genotype quality/g' \
      | sed 's/##ALT=<ID=BND,Description="Translocation">/##ALT=<ID=BND,Description="Breakend">/g' > newheader.vcf

    bcftools reheader -h newheader.vcf ~{vcf} \
      | bcftools annotate -x ~{drop_fields} \
        --no-version \
        -O z \
        -o ~{prefix}.vcf.gz

    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}
