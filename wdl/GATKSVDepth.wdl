version 1.0

import "CollectCoverage.wdl" as cov
import "GATKSVGenotype.wdl" as svg
import "GermlineCNVTasks.wdl" as gcnv_tasks
import "MakeBincovMatrix.wdl" as mbm

workflow GATKSVDepth {
  input {
    String batch
    String cnv_size_name
    File vcf
    Array[String] samples
    Array[File] counts
    File sample_median_count_file
    File ploidy_calls_tar

    # Condense read counts
    Int condense_num_bins
    Int condense_bin_size

    # CNVIntervals
    Boolean include_depth_only
    String cnv_size_conditional
    Int cnv_padding

    # ScatterIntervals
    Int num_intervals_per_scatter

    # Training
    Float? mu_eps
    Float? alpha_ref
    Float? alpha_non_ref
    Float? var_phi
    Int train_max_iter
    String train_device

    # Inference
    Int predictive_samples
    Int predictive_iter
    Int discrete_samples
    String infer_device

    # Reference
    File ref_fasta_dict
    File ref_fasta_fai
    File genome_file

    # Dockers
    String linux_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_docker
    String gatk_docker
    String condense_counts_docker

    # Runtime attributes
    RuntimeAttr? runtime_attr_small_intervals
    RuntimeAttr? runtime_attr_override_make_bincov
    RuntimeAttr? runtime_attr_intersect_intervals
    RuntimeAttr? runtime_attr_counts_to_intervals
    RuntimeAttr? runtime_attr_condense_counts
    RuntimeAttr? runtime_attr_scatter
    RuntimeAttr? runtime_attr_train
    RuntimeAttr? runtime_attr_infer
    RuntimeAttr? runtime_attr_concat
  }

  File vcf_index = vcf + ".tbi"

  if (condense_num_bins > 1) {
    scatter (i in range(length(samples))) {
      call cov.CondenseReadCounts {
        input:
          counts = counts[i],
          sample = samples[i],
          num_bins = condense_num_bins,
          expected_bin_size = condense_bin_size,
          condense_counts_docker = condense_counts_docker,
          runtime_attr_override = runtime_attr_condense_counts
      }
    }
  }

  call CNVIntervals {
    input:
      vcf = vcf,
      vcf_index = vcf_index,
      ref_fasta_fai = ref_fasta_fai,
      size_conditional = cnv_size_conditional,
      include_depth_only = include_depth_only,
      padding = cnv_padding,
      prefix = "~{batch}.~{cnv_size_name}.intervals",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_small_intervals
  }

  scatter (i in range(length(samples))) {
    call IntersectCountsWithIntervals {
      input:
        counts = select_first([CondenseReadCounts.out, counts])[i],
        interval_list = CNVIntervals.out,
        genome_file = genome_file,
        output_name = "~{batch}.~{cnv_size_name}.counts.tsv",
        gzip = true,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_intersect_intervals
    }
  }

  call cov.CountsToIntervals {
    input:
      counts = IntersectCountsWithIntervals.out[0],
      output_name = "~{batch}.~{cnv_size_name}.intersected_intervals",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_counts_to_intervals
  }

  call mbm.MakeBincovMatrix {
    input:
      samples = samples,
      count_files = IntersectCountsWithIntervals.out,
      batch = "~{batch}.~{cnv_size_name}",
      sv_base_mini_docker = sv_base_mini_docker,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_override_make_bincov
  }

  call gcnv_tasks.ScatterIntervals {
    input:
      interval_list = CountsToIntervals.out,
      num_intervals_per_scatter = num_intervals_per_scatter,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_scatter
  }

  scatter (i in range(length(ScatterIntervals.scattered_interval_lists))) {
    String model_name = "~{batch}.depth.~{cnv_size_name}.shard_~{i}"
    call SVTrainDepth {
      input:
        depth_file = MakeBincovMatrix.merged_bincov,
        intervals = ScatterIntervals.scattered_interval_lists[i],
        ploidy_calls_tar = ploidy_calls_tar,
        sample_median_count_file = sample_median_count_file,
        ref_fasta_dict = ref_fasta_dict,
        model_name = model_name,
        mu_eps = mu_eps,
        alpha_ref = alpha_ref,
        alpha_non_ref = alpha_non_ref,
        var_phi = var_phi,
        max_iter = train_max_iter,
        device = train_device,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_train
    }
    call SVInferDepth {
      input:
        model_tar = SVTrainDepth.out,
        ref_fasta_dict = ref_fasta_dict,
        model_name = model_name,
        output_name = model_name,
        predictive_samples = predictive_samples,
        predictive_iter = predictive_iter,
        discrete_samples = discrete_samples,
        gatk_docker = gatk_docker,
        device = infer_device,
        runtime_attr_override = runtime_attr_infer
    }
  }

  call svg.ConcatVcfs {
    input:
      vcfs = SVInferDepth.out,
      vcfs_idx = SVInferDepth.out_index,
      outfile_prefix = "~{batch}.depth.~{cnv_size_name}",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat
  }

  output {
    File out = ConcatVcfs.out
    File out_index = ConcatVcfs.out_index
    File depth_file = MakeBincovMatrix.merged_bincov
    File depth_file_index = MakeBincovMatrix.merged_bincov_idx
  }
}

task CNVIntervals {
  input {
    File vcf
    File vcf_index
    File ref_fasta_fai
    Boolean include_depth_only
    String size_conditional
    String prefix
    Int padding
    String sv_pipeline_docker
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

  output {
    File out = "~{prefix}.intervals.bed"
  }
  command <<<
    set -euo pipefail
    python > intervals.bed <<EOF
import sys
from pysam import VariantFile
vcf = VariantFile('~{vcf}')
types = set(['DEL', 'DUP', 'BND'])
for record in vcf.fetch():
  if record.info['SVTYPE'] in types \
  ~{if !include_depth_only then "and record.info['ALGORITHMS'] != 'depth'" else ""} \
  and record.chrom == record.info['CHR2'] \
  and record.info['STRANDS'][0] != record.info['STRANDS'][1] \
  and record.stop - record.pos ~{size_conditional}:
    fields = [record.chrom, str(record.pos), str(record.stop)]
    print('\t'.join(fields))
EOF

    bedtools slop -b ~{padding} -i intervals.bed -g ~{ref_fasta_fai} \
      | bedtools merge \
      > ~{prefix}.intervals.bed

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

task SVTrainDepth {
  input {
    File depth_file
    File intervals
    File ploidy_calls_tar
    File sample_median_count_file
    File ref_fasta_dict

    Float? mu_eps
    Float? alpha_ref
    Float? alpha_non_ref
    Float? var_phi

    String model_name
    String gatk_docker
    String device
    Int? max_iter
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    depth_file: {
                  localization_optional: true
                }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 15,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{model_name}.model_files.tar.gz"
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
    while read line; do
    echo "--ploidy-calls-file $line" >> args.txt
    done < ploidy_files.list

    mkdir svmodel
    gatk --java-options -Xmx~{java_mem_mb}M SVTrainDepth \
      -L ~{intervals} \
      --depth-file ~{depth_file} \
      --coverage-file ~{sample_median_count_file} \
      --output-name ~{model_name} \
      --output-dir svmodel \
      --arguments_file args.txt \
      --sequence-dictionary ~{ref_fasta_dict} \
      --jit \
      ~{"--max-iter " + max_iter} \
      ~{"--mu-eps " + mu_eps} \
      ~{"--alpha-ref " + alpha_ref} \
      ~{"--alpha-non-ref " + alpha_non_ref} \
      ~{"--var_phi " + var_phi}

    tar czf ~{model_name}.model_files.tar.gz svmodel/*
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

task SVInferDepth {
  input {
    File model_tar
    File ref_fasta_dict
    Int predictive_samples
    Int predictive_iter
    Int discrete_samples
    String model_name
    String output_name
    String gatk_docker
    String device
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 2.5,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 0
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File out = "~{output_name}.vcf.gz"
    File out_index = "~{output_name}.vcf.gz.tbi"
  }
  command <<<

    set -eo pipefail
    mkdir svmodel
    tar xzf ~{model_tar} svmodel/

    gatk --java-options -Xmx~{java_mem_mb}M SVInferDepth \
      --output ~{output_name}.vcf.gz \
      --predictive-samples ~{predictive_samples} \
      --predictive-iter ~{predictive_iter} \
      --discrete-samples ~{discrete_samples} \
      --model-name ~{model_name} \
      --model-dir svmodel \
      --sequence-dictionary ~{ref_fasta_dict} \
      --jit

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

task IntersectCountsWithIntervals {
  input {
    File counts
    File interval_list
    File genome_file
    String output_name
    Boolean gzip = true
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 0.9,
                               disk_gb: 10,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = if gzip then output_name + ".gz" else output_name
  }
  command <<<

    set -euo pipefail
    zgrep -B9999999999 -m1 -v "^@" ~{counts} > ~{output_name}
    zgrep -v "^@" ~{counts} \
      | tail -n +2 \
      | bedtools intersect -wa -sorted -g ~{genome_file} -u -a stdin -b ~{interval_list} \
      >> ~{output_name}
    if ~{gzip}; then
      bgzip ~{output_name}
    fi

  >>>
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
