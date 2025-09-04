version 1.0

import "Structs.wdl"
import "RecalibrateGq.wdl" as recalibrate_gq
import "TasksMakeCohortVcf.wdl" as tasks_cohort

workflow TrainGenotypeFilteringModel {

  input {
    File vcf
    String? output_prefix
    File? truth_json
    File gq_recalibrator_model_file

    Array[String] recalibrate_gq_args = []
    Array[File] genome_tracks = []
    Float fmax_beta = 0.4
    Int optimize_vcf_records_per_shard = 50000

    String linux_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
  }

  String output_prefix_ = if defined(output_prefix) then select_first([output_prefix]) else basename(vcf, ".vcf.gz")

  call recalibrate_gq.RecalibrateGq {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      gq_recalibrator_model_file=gq_recalibrator_model_file,
      recalibrate_gq_args=recalibrate_gq_args,
      genome_tracks=genome_tracks,
      gatk_docker=gatk_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker
  }

  if (defined(truth_json)) {
		call tasks_cohort.ScatterVcf as ScatterForOptimization {
			input:
				vcf=RecalibrateGq.filtered_vcf,
				records_per_shard=optimize_vcf_records_per_shard,
				prefix="~{output_prefix_}.filter_genotypes_scatter",
				sv_pipeline_docker=sv_pipeline_docker
		}

		scatter ( i in range(length(ScatterForOptimization.shards)) ) {
			call MakeVcfTable {
				input:
					vcf=ScatterForOptimization.shards[i],
					truth_json=select_first([truth_json]),
					output_prefix="~{output_prefix_}.vcf_table_shard_~{i}",
					sv_pipeline_docker=sv_pipeline_docker
			}
		}

		call tasks_cohort.ConcatHeaderedTextFiles {
			input:
				text_files=MakeVcfTable.out,
				gzipped=true,
				output_filename="~{output_prefix_}.vcf_table.tsv.gz",
				linux_docker=linux_docker
		}

		call OptimizeCutoffs {
			input:
				table=ConcatHeaderedTextFiles.out,
				fmax_beta=fmax_beta,
				output_prefix="~{output_prefix_}.sl_optimization",
				sv_pipeline_docker=sv_pipeline_docker
		}
	}

  output {
		File unfiltered_recalibrated_vcf = RecalibrateGq.filtered_vcf
		File unfiltered_recalibrated_vcf_index = RecalibrateGq.filtered_vcf_index

		File? vcf_optimization_table = ConcatHeaderedTextFiles.out
    File? sl_cutoff_table = OptimizeCutoffs.sl_cutoff_table
    File? sl_cutoff_qc_tarball = OptimizeCutoffs.out
  }
}

task MakeVcfTable {
  input {
    File vcf
    File truth_json
    File? script
    String? args
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2 + size(truth_json, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.tsv.gz"
  }
  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/scripts/make_sl_table.py" script} \
      --vcf ~{vcf} \
      --truth-json ~{truth_json} \
      --out ~{output_prefix}.tsv.gz \
      ~{args}
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

task OptimizeCutoffs {
  input {
    File table
    File? script
    Float fmax_beta
    String? args
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(100 + size(table, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.tar.gz"
    File sl_cutoff_table = "~{output_prefix}/~{output_prefix}.sl_cutoff_table.tsv"
  }
  command <<<
    set -euo pipefail

    mkdir ~{output_prefix}
    python ~{default="/opt/sv-pipeline/scripts/optimize_sl_cutoffs.py" script} \
      --table ~{table} \
      --out-dir ~{output_prefix} \
      --out-name ~{output_prefix} \
      --beta ~{fmax_beta} \
      ~{args}

    tar czf ~{output_prefix}.tar.gz -C ~{output_prefix}/ .
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