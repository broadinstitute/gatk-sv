version 1.0

import "GATKSVGenotype.wdl" as svg
import "GATKSVDepth.wdl" as gatksv_depth

workflow GATKSVJoinSamples {
  input {
    String batch
    Array[String] samples

    File sv_vcf
    File sr_file
    File pe_file
    Array[File] counts
    Array[File] ploidy_calls
    File sample_mean_depth_file
    File sample_median_count_file

    # Filtering options
    File? inclusion_intervals_depth_only
    File? exclusion_intervals_depth_only
    Int min_size_depth_only = 5000
    Float min_overlap_fraction_depth_only = 0.5
    Boolean require_breakend_overlap_depth_only = false

    File? inclusion_intervals_non_depth_only
    File? exclusion_intervals_non_depth_only
    Int min_size_non_depth_only = 50
    Float min_overlap_fraction_non_depth_only = 0
    Boolean require_breakend_overlap_non_depth_only = true

    # Clustering options
    Float defragment_min_sample_set_fraction_overlap = 0.9
    Float defragment_padding_fraction = 0.5

    String cluster1_algorithm = "SINGLE_LINKAGE"
    String cluster1_breakpoint_summary_strategy = "MEDIAN_START_MEDIAN_END"
    Float cluster1_depth_overlap_fraction = 0.9
    Float cluster1_mixed_overlap_fraction = 0.9
    Float cluster1_pesr_overlap_fraction = 0.9
    Int cluster1_depth_breakend_window = 50
    Int cluster1_mixed_breakend_window = 50
    Int cluster1_pesr_breakend_window = 50

    String cluster2_algorithm = "MAX_CLIQUE"
    String cluster2_breakpoint_summary_strategy = "MEDIAN_START_MEDIAN_END"
    Float cluster2_depth_overlap_fraction = 0.5
    Float cluster2_mixed_overlap_fraction = 0.5
    Float cluster2_pesr_overlap_fraction = 0.5
    Int cluster2_depth_breakend_window = 1000
    Int cluster2_mixed_breakend_window = 500
    Int cluster2_pesr_breakend_window = 500

    # Evidence aggregation options
    Int pesr_aggregation_records_per_shard = 1000
    Int pe_inner_window = 100
    Int pe_outer_window = 500
    Int sr_window = 200

    # Condense read counts
    Int small_cnv_condense_num_bins = 2
    Int small_cnv_condense_bin_size = 200
    Int large_cnv_condense_num_bins = 20
    Int large_cnv_condense_bin_size = 2000

    Int num_intervals_per_scatter = 1000

    Int large_cnv_padding = 1
    Int small_cnv_padding = 1
    Int depth_train_max_iter = 2000
    Int depth_predictive_samples = 100
    Int depth_predictive_iter = 10
    Int depth_discrete_samples = 1000

    Float? depth_mu_eps
    Float? depth_alpha_ref
    Float? depth_alpha_non_ref
    Float? depth_var_phi

    String depth_train_device = 'cpu'
    String depth_infer_device = 'cpu'

    File contig_list
    File ref_fasta_dict
    File ref_fasta_fai
    File ref_fasta
    File genome_file

    String linux_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String sv_pipeline_base_docker
    String gatk_docker
    String condense_counts_docker

    RuntimeAttr? runtime_attr_merge
    RuntimeAttr? runtime_attr_filter_contig
    RuntimeAttr? runtime_attr_filter_non_depth
    RuntimeAttr? runtime_attr_filter_depth_1
    RuntimeAttr? runtime_attr_filter_depth_2
    RuntimeAttr? runtime_attr_shard_clustered
    RuntimeAttr? runtime_attr_aggregate_pe
    RuntimeAttr? runtime_attr_aggregate_sr
    RuntimeAttr? runtime_attr_tar_ploidy_calls
    RuntimeAttr? runtime_attr_concat
    RuntimeAttr? runtime_attr_small_intervals
    RuntimeAttr? runtime_attr_intersect_intervals
    RuntimeAttr? runtime_attr_counts_to_intervals
    RuntimeAttr? runtime_attr_defrag
    RuntimeAttr? runtime_attr_cluster_1
    RuntimeAttr? runtime_attr_cluster_2
    RuntimeAttr? runtime_attr_posteriors
    RuntimeAttr? runtime_attr_condense_counts
    RuntimeAttr? runtime_attr_override_make_bincov
    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_train
    RuntimeAttr? runtime_attr_infer
  }

  File sr_index_ = sr_file + ".tbi"
  File pe_index_ = pe_file + ".tbi"

  Array[String] contigs = read_lines(contig_list)

  scatter (contig in contigs) {
    # Get records on this contig
    call SelectVariants as FilterContig {
      input:
        vcf = sv_vcf,
        vcf_index = sv_vcf + ".tbi",
        output_name = "~{batch}.~{contig}.filter_contig",
        contig = contig,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_filter_contig
    }
    # Filter non-depth-only records by size and location
    call FilterSVs as FilterNonDepth {
      input:
        vcf = FilterContig.out,
        vcf_index = FilterContig.out_index,
        output_name = "~{batch}.~{contig}.filter_non_depth",
        non_depth_only = true,
        min_size = min_size_non_depth_only,
        min_overlap_fraction = min_overlap_fraction_non_depth_only,
        require_breakend_overlap = require_breakend_overlap_non_depth_only,
        inclusion_intervals = inclusion_intervals_non_depth_only,
        exclusion_intervals = exclusion_intervals_non_depth_only,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_filter_non_depth
    }
    # Get depth-only records
    call FilterSVs as FilterDepth1 {
      input:
        vcf = FilterContig.out,
        vcf_index = FilterContig.out_index,
        output_name = "~{batch}.~{contig}.filter_depth_1",
        depth_only = true,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_filter_depth_1
    }
    # Defragment depth-only records
    call ClusterVariants as Defragment {
      input:
        vcf = FilterDepth1.out,
        vcf_index = FilterDepth1.out_index,
        output_name = "~{batch}.~{contig}.defragment",
        vid_prefix = "SV_" + contig + "_",
        contig = contig,
        algorithm = "DEFRAGMENT_CNV",
        min_sample_set_fraction_overlap = defragment_min_sample_set_fraction_overlap,
        defrag_padding_fraction = defragment_padding_fraction,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_defrag
    }
    # Filter defragmented depth-only records by size and location
    call FilterSVs as FilterDepth2 {
      input:
        vcf = Defragment.out,
        vcf_index = Defragment.out_index,
        output_name = "~{batch}.~{contig}.filter_depth_2",
        inclusion_intervals = inclusion_intervals_depth_only,
        exclusion_intervals = exclusion_intervals_depth_only,
        min_size = min_size_depth_only,
        min_overlap_fraction = min_overlap_fraction_depth_only,
        require_breakend_overlap = require_breakend_overlap_depth_only,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_filter_depth_2
    }
    # Merge depth and non-depth records
    call svg.ConcatVcfs as ConcatDepthAndNonDepth {
      input:
        vcfs = [FilterNonDepth.out, FilterDepth2.out],
        vcfs_idx = [FilterNonDepth.out_index, FilterDepth2.out_index],
        outfile_prefix = "~{batch}.~{contig}.concat_depth_and_non_depth",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_concat
    }
    # First cluster step (default: strict single-linkage)
    call ClusterVariants as Cluster1 {
      input:
        vcf = ConcatDepthAndNonDepth.out,
        vcf_index = ConcatDepthAndNonDepth.out_index,
        vid_prefix = "SV_" + contig + "_",
        algorithm = cluster1_algorithm,
        convert_inv = true,
        breakpoint_summary_strategy = cluster1_breakpoint_summary_strategy,
        depth_overlap_fraction = cluster1_depth_overlap_fraction,
        mixed_overlap_fraction = cluster1_mixed_overlap_fraction,
        pesr_overlap_fraction = cluster1_pesr_overlap_fraction,
        depth_breakend_window = cluster1_depth_breakend_window,
        mixed_breakend_window = cluster1_mixed_breakend_window,
        pesr_breakend_window = cluster1_pesr_breakend_window,
        output_name = "~{batch}.~{contig}.cluster_1",
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_cluster_1
    }
    # Second cluster step (default: less strict max-clique)
    call ClusterVariants as Cluster2 {
      input:
        vcf = Cluster1.out,
        vcf_index = Cluster1.out_index,
        vid_prefix = "SV_" + contig + "_",
        algorithm = cluster2_algorithm,
        breakpoint_summary_strategy = cluster2_breakpoint_summary_strategy,
        depth_overlap_fraction = cluster2_depth_overlap_fraction,
        mixed_overlap_fraction = cluster2_mixed_overlap_fraction,
        pesr_overlap_fraction = cluster2_pesr_overlap_fraction,
        depth_breakend_window = cluster2_depth_breakend_window,
        mixed_breakend_window = cluster2_mixed_breakend_window,
        pesr_breakend_window = cluster2_pesr_breakend_window,
        output_name = "~{batch}.~{contig}.cluster_2",
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_cluster_2
    }
  }
  # Concatenate contig shards
  call svg.ConcatVcfs as ConcatClustered {
    input:
      vcfs = Cluster2.out,
      vcfs_idx = Cluster2.out_index,
      outfile_prefix = "~{batch}.clustered",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  call svg.ShardVcf as ShardClustered {
    input:
      vcf = ConcatClustered.out,
      vcf_index = ConcatClustered.out_index,
      records_per_shard = pesr_aggregation_records_per_shard,
      basename = "~{batch}.clustered",
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_shard_clustered
  }
  scatter (i in range(length(ShardClustered.out))) {
    # Collect PE evidence
    call AggregatePESREvidence as AggregatePE {
      input:
        vcf = ShardClustered.out[i],
        output_name = "~{batch}.shard_~{i}.aggregate_pe",
        pe_file = pe_file,
        pe_index = pe_index_,
        pe_inner_window = pe_inner_window,
        pe_outer_window = pe_outer_window,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_aggregate_pe
    }
    # Collect SR evidence
    call AggregatePESREvidence as AggregateSR {
      input:
        vcf = AggregatePE.out,
        output_name = "~{batch}.shard_~{i}.aggregate_sr",
        sr_file = sr_file,
        sr_index = sr_index_,
        sample_mean_depth_file = sample_mean_depth_file,
        sr_window = sr_window,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_aggregate_sr
    }
  }
  # Concatenate PESR aggreagtion shards
  call svg.ConcatVcfs as ConcatAggregated {
    input:
      vcfs = AggregateSR.out,
      vcfs_idx = AggregateSR.out_index,
      outfile_prefix = "~{batch}.aggregate_pesr",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  # Combine ploidy calls into a single file
  call TarFiles as TarPloidyCalls {
    input:
      files = ploidy_calls,
      name = "~{batch}.ploidy_calls",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_tar_ploidy_calls
  }
  # Depth model for large CNVs
  call gatksv_depth.GATKSVDepth as DepthLarge {
    input:
      batch = batch,
      cnv_size_name = "large_cnv",
      vcf = ConcatAggregated.out,
      samples = samples,
      counts = counts,
      sample_median_count_file = sample_median_count_file,
      ploidy_calls_tar = TarPloidyCalls.out,
      condense_num_bins = large_cnv_condense_num_bins,
      condense_bin_size = large_cnv_condense_bin_size,
      include_depth_only = true,
      cnv_size_conditional = ">= 5000",
      cnv_padding = large_cnv_padding,
      num_intervals_per_scatter = num_intervals_per_scatter,
      alpha_ref = depth_alpha_ref,
      train_max_iter = depth_train_max_iter,
      train_device = depth_train_device,
      predictive_samples = depth_predictive_samples,
      predictive_iter = depth_predictive_iter,
      discrete_samples = depth_discrete_samples,
      infer_device = depth_infer_device,
      ref_fasta_dict = ref_fasta_dict,
      ref_fasta_fai = ref_fasta_fai,
      genome_file = genome_file,
      linux_docker = linux_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      gatk_docker = gatk_docker,
      condense_counts_docker = condense_counts_docker,
      runtime_attr_small_intervals = runtime_attr_small_intervals,
      runtime_attr_override_make_bincov = runtime_attr_override_make_bincov,
      runtime_attr_intersect_intervals = runtime_attr_intersect_intervals,
      runtime_attr_counts_to_intervals = runtime_attr_counts_to_intervals,
      runtime_attr_condense_counts = runtime_attr_condense_counts,
      runtime_attr_scatter = runtime_attr_scatter,
      runtime_attr_train = runtime_attr_train,
      runtime_attr_infer = runtime_attr_infer,
      runtime_attr_concat = runtime_attr_concat
  }
  # Depth model for small CNVs
  call gatksv_depth.GATKSVDepth as DepthSmall {
    input:
      batch = batch,
      cnv_size_name = "small_cnv",
      vcf = ConcatAggregated.out,
      samples = samples,
      counts = counts,
      sample_median_count_file = sample_median_count_file,
      ploidy_calls_tar = TarPloidyCalls.out,
      condense_num_bins = small_cnv_condense_num_bins,
      condense_bin_size = small_cnv_condense_bin_size,
      include_depth_only = false,
      cnv_size_conditional = "< 5000",
      cnv_padding = small_cnv_padding,
      num_intervals_per_scatter = num_intervals_per_scatter,
      mu_eps = depth_mu_eps,
      alpha_ref = depth_alpha_ref,
      alpha_non_ref = depth_alpha_non_ref,
      var_phi = depth_var_phi,
      train_max_iter = depth_train_max_iter,
      train_device = depth_train_device,
      predictive_samples = depth_predictive_samples,
      predictive_iter = depth_predictive_iter,
      discrete_samples = depth_discrete_samples,
      infer_device = depth_infer_device,
      ref_fasta_dict = ref_fasta_dict,
      ref_fasta_fai = ref_fasta_fai,
      genome_file = genome_file,
      linux_docker = linux_docker,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      gatk_docker = gatk_docker,
      condense_counts_docker = condense_counts_docker,
      runtime_attr_small_intervals = runtime_attr_small_intervals,
      runtime_attr_override_make_bincov = runtime_attr_override_make_bincov,
      runtime_attr_intersect_intervals = runtime_attr_intersect_intervals,
      runtime_attr_counts_to_intervals = runtime_attr_counts_to_intervals,
      runtime_attr_condense_counts = runtime_attr_condense_counts,
      runtime_attr_scatter = runtime_attr_scatter,
      runtime_attr_train = runtime_attr_train,
      runtime_attr_infer = runtime_attr_infer,
      runtime_attr_concat = runtime_attr_concat
  }
  # Aggregate depth evidence
  call AggregateDepth {
    input:
      vcf = ConcatAggregated.out,
      vcf_index = ConcatAggregated.out_index,
      output_name = "~{batch}.aggregate_depth",
      depth_posterior_vcfs = [DepthLarge.out, DepthSmall.out],
      depth_posterior_vcfs_indexes = [DepthLarge.out_index, DepthSmall.out_index],
      ploidy_calls_tar = TarPloidyCalls.out,
      ref_fasta_dict = ref_fasta_dict,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_posteriors
  }

  output {
    File joined_vcf = AggregateDepth.out
    File joined_vcf_index = AggregateDepth.out_index

    File large_cnv_depth_vcf = DepthLarge.out
    File large_cnv_depth_vcf_index = DepthLarge.out_index
    File large_cnv_depth_counts = DepthLarge.depth_file
    File large_cnv_depth_counts_index = DepthLarge.depth_file_index

    File small_cnv_depth_vcf = DepthSmall.out
    File small_cnv_depth_vcf_index = DepthSmall.out_index
    File small_cnv_depth_counts = DepthSmall.depth_file
    File small_cnv_depth_counts_index = DepthSmall.depth_file_index
  }
}

task AggregateDepth {
  input {
    File vcf
    File vcf_index
    String output_name
    File ploidy_calls_tar
    Array[File] depth_posterior_vcfs
    Array[File] depth_posterior_vcfs_indexes
    File ploidy_calls_tar
    File ref_fasta_dict
    String gatk_path = "/gatk/gatk"
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 35,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_name}.vcf.gz"
    File out_index = "~{output_name}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail

    # Extract ploidy call tarballs
    mkdir ploidy-calls
    tar xzf ~{ploidy_calls_tar} -C ploidy-calls
    cd ploidy-calls/
    for file in *.tar.gz; do
      name=$(basename $file .tar.gz)
      mkdir $name
      tar xzf $file -C $name/
    done
    cd ../
    ls ploidy-calls/*/contig_ploidy.tsv > ploidy_files.list

    # Create arguments file
    echo "--cnv-intervals-vcf ~{sep=" --cnv-intervals-vcf " depth_posterior_vcfs}" > args.txt
    while read line; do
      echo "--ploidy-calls-file $line" >> args.txt
    done < ploidy_files.list

    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" SVAggregateDepth \
      --arguments_file args.txt \
      --variant ~{vcf} \
      --output ~{output_name}.vcf.gz \
      --sequence-dictionary ~{ref_fasta_dict} \
      --genotype-depth-calls

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ClusterVariants {
  input {
    File vcf
    File vcf_index
    String? contig
    String output_name

    String algorithm
    String vid_prefix
    Boolean convert_inv = false
    Boolean enable_cnv = false

    Float? min_sample_set_fraction_overlap
    Float? defrag_padding_fraction

    String? breakpoint_summary_strategy

    Float? depth_overlap_fraction
    Float? mixed_overlap_fraction
    Float? pesr_overlap_fraction
    Int? depth_breakend_window
    Int? mixed_breakend_window
    Int? pesr_breakend_window

    String gatk_path = "/gatk/gatk"
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcf: {
      localization_optional: true
    }
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

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_name}.vcf.gz"
    File out_index = "~{output_name}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail
    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" SVCluster \
      --disable-sequence-dictionary-validation \
      -V ~{vcf} \
      -O ~{output_name}.vcf.gz \
      --algorithm ~{algorithm} \
      --variant-prefix ~{vid_prefix} \
      ~{"-L " + contig} \
      ~{"--breakpoint-summary-strategy " + breakpoint_summary_strategy} \
      ~{"--min-sample-set-fraction-overlap " + min_sample_set_fraction_overlap} \
      ~{"--defrag-padding-fraction " + defrag_padding_fraction} \
      ~{"--depth-overlap-fraction " + depth_overlap_fraction} \
      ~{"--mixed-overlap-fraction " + mixed_overlap_fraction} \
      ~{"--pesr-overlap-fraction " + pesr_overlap_fraction} \
      ~{"--depth-breakend-window " + depth_breakend_window} \
      ~{"--mixed-breakend-window " + mixed_breakend_window} \
      ~{"--pesr-breakend-window " + pesr_breakend_window} \
      ~{if convert_inv then "--convert-inv-to-bnd" else ""} \
      ~{if enable_cnv then "--enable-cnv" else ""}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task AggregatePESREvidence {
  input {
    File vcf
    String output_name

    File? sr_file
    File? sr_index
    File? sample_mean_depth_file

    File? pe_file
    File? pe_index

    Int? pe_inner_window
    Int? pe_outer_window
    Int? sr_window

    String gatk_path = "/gatk/gatk"
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    sr_file: {
               localization_optional: true
             }
    pe_file: {
               localization_optional: true
             }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 15,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_name}.vcf.gz"
    File out_index = "~{output_name}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail
    tabix ~{vcf}
    ~{gatk_path} --java-options "-Xmx~{java_mem_mb}m" AggregatePairedEndAndSplitReadEvidence \
      -V ~{vcf} \
      -O ~{output_name}.vcf.gz \
      ~{"--split-reads-file " + sr_file} \
      ~{"--discordant-pairs-file " + pe_file} \
      ~{"--sample-coverage " + sample_mean_depth_file} \
      ~{"--pe-inner-window " + pe_inner_window} \
      ~{"--pe-outer-window " + pe_outer_window} \
      ~{"--sr-window " + sr_window}
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task FilterSVs {
  input {
    File vcf
    File vcf_index

    String output_name
    String? contig

    Boolean depth_only = false
    Boolean non_depth_only = false
    Boolean require_breakend_overlap = false
    File? inclusion_intervals
    File? exclusion_intervals
    Int? min_size
    Float? min_overlap_fraction

    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.0,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_name}.vcf.gz"
    File out_index = "~{output_name}.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    gatk --java-options -Xmx~{java_mem_mb}M SVSelectVariants \
      -V ~{vcf} \
      -O ~{output_name}.vcf.gz \
      ~{if depth_only then "--depth-only" else ""} \
      ~{if non_depth_only then "--non-depth-only" else ""} \
      ~{"-L " + inclusion_intervals} \
      ~{"-XL " + exclusion_intervals} \
      ~{"--min-size " + min_size} \
      ~{"--min-overlap-fraction " + min_overlap_fraction} \
      ~{if require_breakend_overlap then "--require-breakend-overlap" else ""}

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

task SelectVariants {
  input {
    File vcf
    File vcf_index
    String output_name
    String? contig
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_name}.vcf.gz"
    File out_index = "~{output_name}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail
    gatk --java-options -Xmx~{java_mem_mb}M SelectVariants \
      -V ~{vcf} \
      -O ~{output_name}.vcf.gz \
      ~{"-L " + contig}
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

task TarFiles {
  input {
    Array[File] files
    String name
    String linux_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "tmp/~{name}.tar.gz"
  }
  command <<<
    mkdir tmp
    cd tmp
    while read -r file; do
      mv $file .
    done < ~{write_lines(files)}
    tar czf ~{name}.tar.gz *
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}