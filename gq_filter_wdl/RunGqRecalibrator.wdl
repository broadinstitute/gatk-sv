version development

import "TrainGqRecalibrator.wdl" as TrainGqRecalibrator

workflow RunGqRecalibrator {
    input {
        File vcf
        File? vcf_index
        Array[File] genome_tracts
        File gq_recalibrator_model_file
        Array[String] recalibrate_gq_args = []
        Boolean fix_vcf = true
        String gatk_docker
        String module03_docker
        String sv_base_docker
    }

    if(fix_vcf) {
        call TrainGqRecalibrator.FixVcf as FixVcf {
            input:
                vcf=vcf,
                module03_docker=module03_docker
        }
    }

    if(!defined(vcf_index) && !fix_vcf) {
        call TrainGqRecalibrator.IndexVcf as IndexVcf {
            input:
                vcf=vcf,
                sv_base_docker=sv_base_docker
        }
    }

    File vcf_ = select_first([FixVcf.fixed_vcf, vcf])
    File vcf_index_ = select_first([FixVcf.fixed_vcf_index, vcf_index, IndexVcf.index_file])

    call ApplyGqRecalibratorFilter {
        input:
            vcf=vcf_,
            vcf_index=vcf_index_,
            genome_tracts=genome_tracts,
            gq_recalibrator_model_file=gq_recalibrator_model_file,
            recalibrate_gq_args = recalibrate_gq_args,
            gatk_docker=gatk_docker
    }

    output {
        File filtered_vcf = ApplyGqRecalibratorFilter.filtered_vcf
        File filtered_vcf_index = ApplyGqRecalibratorFilter.filtered_vcf_index
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

    Int disk_gb = round(1000 + size([vcf, vcf_index, gq_recalibrator_model_file], "GiB") + size(genome_tracts, "GiB"))
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
