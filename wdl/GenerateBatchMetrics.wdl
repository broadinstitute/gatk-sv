version 1.0

import "RDTest.wdl" as rdt
import "TasksMakeCohortVcf.wdl" as taskscohort
import "Utils.wdl" as util
import "TestUtils.wdl" as tu

workflow GenerateBatchMetrics {
  input {
    String batch

    File depth_vcf
    File? melt_vcf
    File? scramble_vcf
    File? wham_vcf
    File? manta_vcf

    File pe_file
    File sr_file
    File baf_file
    File rd_file

    File median_file
    File mean_coverage_file
    File ploidy_table

    Int records_per_shard_agg
    Int records_per_shard_rdtest

    String? additional_gatk_args_agg

    File? svtk_to_gatk_script

    String chr_x
    String chr_y

    Float? java_mem_fraction

    File rmsk
    File segdups
    File ped_file
    File autosome_contigs
    File allosome_contigs
    File reference_dict

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    File? primary_contigs_list  # required if run_module_metrics = true

    String gatk_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_base_docker
    String linux_docker

    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_subset_ped
    RuntimeAttr? runtime_attr_sample_list
    RuntimeAttr? runtime_attr_aggregate_tests
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_agg
    RuntimeAttr? runtime_attr_split_rd_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
    RuntimeAttr? runtime_attr_get_male_only
    RuntimeAttr? runtime_attr_metrics_file_metrics
    RuntimeAttr? runtime_attr_annotate_overlap
  }

  String prefix = "~{batch}.batch_metrics"

  call util.GetSampleIdsFromVcf {
    input:
      vcf = depth_vcf,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  call util.SubsetPedFile {
    input:
      ped_file = ped_file,
      sample_list = GetSampleIdsFromVcf.out_file,
      subset_name = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_ped
  }

  call GetSampleLists {
    input:
      ped_file = SubsetPedFile.ped_subset_file,
      samples_list = GetSampleIdsFromVcf.out_file,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_sample_list
  }

  Array[File] vcfs_ = select_all([depth_vcf, manta_vcf, melt_vcf, scramble_vcf, wham_vcf])
  scatter (i in range(length(vcfs_))) {
    File vcfs_index_ = vcfs_[i] + ".tbi"
  }

  call taskscohort.ConcatVcfs as ConcatInputVcfs {
    input:
      vcfs=vcfs_,
      vcfs_idx=vcfs_index_,
      allow_overlaps=true,
      outfile_prefix="~{prefix}.concat_input_vcfs",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_vcfs
  }

  call taskscohort.ScatterVcf {
    input:
      vcf=ConcatInputVcfs.concat_vcf,
      records_per_shard = records_per_shard_agg,
      prefix = "~{prefix}.scatter_vcf",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_scatter_vcf
  }

  scatter ( i in range(length(ScatterVcf.shards)) ) {
    call FormatVcfForGatk {
      input:
        vcf=ScatterVcf.shards[i],
        ploidy_table=ploidy_table,
        output_prefix="~{prefix}.format.shard_~{i}",
        script=svtk_to_gatk_script,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_format
    }
    call SVRegionOverlap {
      input:
        vcf = FormatVcfForGatk.out,
        vcf_index = FormatVcfForGatk.out_index,
        reference_dict = reference_dict,
        output_prefix = "~{prefix}.region_overlap.shard_~{i}",
        region_files = [segdups, rmsk],
        region_file_indexes = [segdups + ".tbi", rmsk + ".tbi"],
        region_names = ["SEGDUP", "RMSK"],
        java_mem_fraction=java_mem_fraction,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_annotate_overlap
    }
    call AggregateSVEvidence {
      input:
        vcf = SVRegionOverlap.out,
        vcf_index = SVRegionOverlap.out_index,
        output_prefix = "~{prefix}.aggregate.shard_~{i}",
        mean_coverage_file = mean_coverage_file,
        ploidy_table=ploidy_table,
        pe_file = pe_file,
        pe_file_index = if defined(pe_file) then select_first([pe_file]) + ".tbi" else pe_file,
        sr_file = sr_file,
        sr_file_index = if defined(sr_file) then select_first([sr_file]) + ".tbi" else sr_file,
        baf_file = baf_file,
        baf_file_index = if defined(baf_file) then select_first([baf_file]) + ".tbi" else baf_file,
        chr_x = chr_x,
        chr_y = chr_y,
        additional_args=additional_gatk_args_agg,
        java_mem_fraction = java_mem_fraction,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_agg
    }
  }

  call taskscohort.ConcatVcfs as ConcatOutputVcfs {
    input:
      vcfs=AggregateSVEvidence.out,
      vcfs_idx=AggregateSVEvidence.out_index,
      naive=true,
      outfile_prefix="~{prefix}.concat",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_vcfs
  }

  call GetMaleOnlyVariantIDs {
    input:
      vcf = ConcatInputVcfs.concat_vcf,
      female_samples = GetSampleLists.female_list,
      male_samples = GetSampleLists.male_list,
      contig = select_first([chr_x, "chrX"]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_get_male_only
  }

  call rdt.RDTest {
    input:
      vcf = ConcatInputVcfs.concat_vcf,
      prefix = "~{prefix}.rdtest",
      coveragefile = select_first([rd_file]),
      medianfile = median_file,
      ped_file = ped_file,
      autosome_contigs = autosome_contigs,
      split_size = records_per_shard_rdtest,
      flags = "",
      allosome_contigs = allosome_contigs,
      ref_dict = reference_dict,
      samples = GetSampleIdsFromVcf.out_file,
      male_samples = GetSampleLists.male_list,
      female_samples = GetSampleLists.female_list,
      male_only_variant_ids = GetMaleOnlyVariantIDs.male_only_variant_ids,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
      linux_docker = linux_docker,
      runtime_attr_rdtest = runtime_attr_rdtest,
      runtime_attr_split_rd_vcf = runtime_attr_split_rd_vcf,
      runtime_attr_merge_allo = runtime_attr_merge_allo,
      runtime_attr_merge_stats = runtime_attr_merge_stats
  }

  call AggregateTests {
    input:
      vcf = ConcatOutputVcfs.concat_vcf,
      vcf_index = ConcatOutputVcfs.concat_vcf_idx,
      prefix = "~{prefix}.aggregate_tests",
      rdtest = RDTest.rdtest,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_aggregate_tests
  }

  Boolean run_module_metrics_ = if defined(run_module_metrics) then select_first([run_module_metrics]) else true
  if (run_module_metrics_) {
    call tu.MetricsFileMetrics {
      input:
        metrics_file = AggregateTests.out,
        contig_list = select_first([primary_contigs_list]),
        common = false,
        prefix = "GenerateBatchMetrics.~{batch}",
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_metrics_file_metrics
    }
  }

  output {
    File metrics = AggregateTests.out
    File? metrics_file_batchmetrics = MetricsFileMetrics.out
  }
}

task GetSampleLists {
  input {
    File ped_file
    File samples_list
    String sv_base_docker
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
    File male_list = "male.list"
    File female_list = "female.list"
  }
  command <<<

    set -eu
    awk -v sex=1 '($5==sex) {print $2}' ~{ped_file} > ped_males.list
    awk -v sex=2 '($5==sex) {print $2}' ~{ped_file} > ped_females.list

    python3 <<CODE
    with open("ped_males.list",'r') as ped_m, open("ped_females.list",'r') as ped_f:
      male_samples = set([x.strip() for x in ped_m.readlines() if x.strip()])
      female_samples = set([x.strip() for x in ped_f.readlines() if x.strip()])
      with open("male.list", 'w') as samples_m, open("female.list",'w') as samples_f, open("~{samples_list}",'r') as samples:
        for line in samples:
          if line.strip():
            if (line.strip() in male_samples):
              samples_m.write(line)
            if (line.strip() in female_samples):
              samples_f.write(line)
    CODE
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task AggregateCallers {
  input {
    String batch
    Array[File] input_metrics
    String sv_pipeline_base_docker
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

  String output_file = "${batch}.metrics"

  output {
    File metrics = "~{output_file}"
  }
  command <<<

    set -eu
    python3 <<CODE
    import pandas as pd
    metrics = ["~{sep='", "' input_metrics}"]
    dfs=[]
    for df in metrics:
      dfs.append(pd.read_table(df))
    df = pd.concat(dfs)
    df.to_csv("~{output_file}", index=False, sep='\t')
    CODE
        
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task AggregateSVEvidence {
  input {
    File vcf
    File vcf_index
    String output_prefix

    File mean_coverage_file
    File ploidy_table
    File? pe_file
    File? pe_file_index
    File? sr_file
    File? sr_file_index
    File? baf_file
    File? baf_file_index

    String chr_x
    String chr_y

    String? additional_args

    Float? java_mem_fraction
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }


  parameter_meta {
    pe_file: {
               localization_optional: true
             }
    sr_file: {
               localization_optional: true
             }
    baf_file: {
                localization_optional: true
              }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2.5 + size(sr_file, "GB")),
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

    function getJavaMem() {
      # get JVM memory in MiB by getting total memory from /proc/meminfo
      # and multiplying by java_mem_fraction
      cat /proc/meminfo \
        | awk -v MEM_FIELD="$1" '{
          f[substr($1, 1, length($1)-1)] = $2
        } END {
          printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
        }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)
    echo "JVM memory: $JVM_MAX_MEM"

    gatk --java-options "-Xmx${JVM_MAX_MEM}" AggregateSVEvidence \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      --sample-coverage ~{mean_coverage_file} \
      --ploidy-table ~{ploidy_table} \
      --x-chromosome-name ~{chr_x} \
      --y-chromosome-name ~{chr_y} \
      ~{"--discordant-pairs-file " + pe_file} \
      ~{"--split-reads-file " + sr_file} \
      ~{"--baf-file " + baf_file} \
      ~{additional_args}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GetMaleOnlyVariantIDs {
  input {
    File vcf
    File female_samples
    File male_samples
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File male_only_variant_ids = "male_only_variant_ids.txt"
  }
  command <<<
    set -euxo pipefail
    bcftools view -t ~{contig} -S ~{male_samples} ~{vcf} | bcftools view --min-ac 1 | bcftools query -f '%ID\n' > variant_ids_in_males.txt
    bcftools view -t ~{contig} -S ~{female_samples} ~{vcf} | bcftools view --min-ac 1 | bcftools query -f '%ID\n' > variant_ids_in_females.txt
    awk 'NR==FNR{a[$0];next} !($0 in a)' variant_ids_in_females.txt variant_ids_in_males.txt > male_only_variant_ids.txt
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

task FormatVcfForGatk {
  input {
    File vcf
    File ploidy_table
    File? script
    String? remove_infos
    String? remove_formats
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2.0),
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
    python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" script} \
      --vcf ~{vcf} \
      --out ~{output_prefix}.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      ~{"--remove-infos " + remove_infos} \
      ~{"--remove-formats " + remove_formats} \
      --fix-end
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

task SVRegionOverlap {
  input {
    File vcf
    File vcf_index
    File reference_dict
    String output_prefix
    Array[File] region_files
    Array[File] region_file_indexes
    Array[String] region_names

    String? region_set_rule
    String? region_merging_rule
    Int? region_padding

    Boolean? suppress_overlap_fraction
    Boolean? suppress_endpoint_counts

    Float? java_mem_fraction

    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2.0),
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

    function getJavaMem() {
      # get JVM memory in MiB by getting total memory from /proc/meminfo
      # and multiplying by java_mem_fraction
      cat /proc/meminfo \
        | awk -v MEM_FIELD="$1" '{
          f[substr($1, 1, length($1)-1)] = $2
        } END {
          printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
        }'
      }
    JVM_MAX_MEM=$(getJavaMem MemTotal)
    echo "JVM memory: $JVM_MAX_MEM"

    gatk --java-options "-Xmx${JVM_MAX_MEM}" SVRegionOverlap \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      --sequence-dictionary ~{reference_dict} \
      --region-file ~{sep=" --region-file " region_files} \
      --region-name ~{sep=" --region-name " region_names} \
      ~{"--region-set-rule " + region_set_rule} \
      ~{"--region-merging-rule " + region_merging_rule} \
      ~{"--region-padding " + region_padding} \
      --suppress-overlap-fraction ~{default="false" suppress_overlap_fraction} \
      --suppress-endpoint-counts ~{default="false" suppress_endpoint_counts}

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task AggregateTests {
  input {
    File vcf
    File vcf_index
    String prefix
    File? rdtest
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(50 + size(vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{prefix}.metrics.tsv"
  }
  command <<<
    /opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py \
      -v ~{vcf} \
      ~{"-r " + rdtest} \
      ~{prefix}.metrics.tsv
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