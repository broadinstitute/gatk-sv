version development

import "TasksRunGqRecalibrator.wdl" as Utils
import "TrainGqRecalibrator.wdl" as TrainGqRecalibrator

workflow RunGqRecalibrator {
    input {
        File vcf
        File? vcf_index
        Array[File] genome_tracts
        File gq_recalibrator_model_file
        Array[String] recalibrate_gq_args = []
        Boolean fix_vcf = true
        Boolean shard_vcf = false
        Int records_per_shard = 5000
        Float? recalibrator_mem_gb_overhead = 1.5
        String gatk_docker
        String module03_docker
        String sv_base_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_override_shard_vcf
        RuntimeAttr? runtime_override_fix_vcf
        RuntimeAttr? runtime_override_concat_vcf
    }

    if(shard_vcf) {
        call Utils.SplitVcf as ShardVcf {
            input:
                vcf=vcf,
                vcf_idx=vcf_index,
                prefix=basename(vcf, ".vcf.gz") + "_sharded",
                min_vars_per_shard=records_per_shard,
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_override_shard_vcf
        }
    }

    Array[File] vcf_shards_ = select_first([ShardVcf.vcf_shards, [vcf]])

    scatter ( shard in vcf_shards_ ) {
        if(fix_vcf) {
            call TrainGqRecalibrator.FixVcf as FixVcf {
                input:
                    vcf=shard,
                    module03_docker=module03_docker,
                    runtime_attr_override=runtime_override_fix_vcf
            }
        }

        if(!fix_vcf) {
            call TrainGqRecalibrator.IndexVcf as IndexVcf {
                input:
                    vcf=shard,
                    sv_base_docker=sv_base_docker
            }
        }

        File vcf_ = select_first([FixVcf.fixed_vcf, shard])
        File vcf_index_ = select_first([FixVcf.fixed_vcf_index, IndexVcf.index_file])

        call ApplyGqRecalibratorFilter {
            input:
                vcf=vcf_,
                vcf_index=vcf_index_,
                genome_tracts=genome_tracts,
                gq_recalibrator_model_file=gq_recalibrator_model_file,
                recalibrate_gq_args = recalibrate_gq_args,
                gatk_docker=gatk_docker
        }
    }

    call Utils.ConcatVcfs as ConcatShards {
        input:
            vcfs=ApplyGqRecalibratorFilter.filtered_vcf,
            vcfs_idx=ApplyGqRecalibratorFilter.filtered_vcf_index,
            outfile_prefix=basename(vcf, ".vcf.gz") + ".recalibrated_GQ",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_concat_vcf
    }

    output {
        File filtered_vcf = ConcatShards.concat_vcf
        File filtered_vcf_index = ConcatShards.concat_vcf_idx
    }
}


task ApplyGqRecalibratorFilter {
    input {
        File vcf
        File vcf_index
        Array[File] genome_tracts
        File gq_recalibrator_model_file
        Array[String] recalibrate_gq_args = []
        String gatk_docker
        Float mem_gb_java = 9.0
        Float mem_gb_overhead = 1.5
    }

    String args_str = if length(recalibrate_gq_args) > 0 then sep(" ", recalibrate_gq_args) else ""

    Int base_disk_gb = 30
    Int disk_gb = round(base_disk_gb + size([vcf, vcf_index, gq_recalibrator_model_file], "GiB") + size(genome_tracts, "GiB"))
    Float mem_gb = mem_gb_java + mem_gb_overhead
    String filtered_vcf_name = sub(sub(basename(vcf), ".gz", ""), ".vcf", "_gq_recalibrated.vcf.gz")

    runtime {
        docker: gatk_docker
        cpu: 1
        preemptible: 3
        max_retries: 0
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail

        ln -s ~{vcf} .
        ln -s ~{vcf_index} .

        mem_kb_java_actual=$(grep -m1 MemTotal /proc/meminfo \
                             | awk '{printf("%.0f\n", $2 - ~{mem_gb_overhead} * 1048576)}')

        gatk --java-options "-Xmx${mem_kb_java_actual}K" XGBoostMinGqVariantFilter \
          --mode "Filter" \
          --variant ~{basename(vcf)} \
          --genome-tract ~{sep=" --genome-tract " genome_tracts} \
          --model-file ~{gq_recalibrator_model_file} \
          --output ~{filtered_vcf_name} \
          ~{args_str}
    >>>

    output {
        File filtered_vcf = filtered_vcf_name
        File filtered_vcf_index = filtered_vcf_name + ".tbi"
    }
}
