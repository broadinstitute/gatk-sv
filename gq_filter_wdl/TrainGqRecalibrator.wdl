version development

import "GetTruthOverlap.wdl" as GetTruthOverlap

workflow TrainGqRecalibrator {
    input {
        File train_vcf
        File? train_vcf_index
        Array[File] truth_vcfs
        Array[File]? truth_vcf_indices
        Map[String, File]? vapor_files
        File ped_file
        Array[File] genome_tracts
        File? optimal_overlap_cutoffs
        File? gq_recalibrator_model_file
        Array[String] train_args = []
        Array[String] get_truth_overlap_args = []
        Boolean fix_vcf = true
        String module03_docker
        String gatk_docker
        String sv_base_docker
    }

    if(fix_vcf) {
        call FixVcf {
            input:
                vcf=train_vcf,
                module03_docker=module03_docker,
                index_output_vcf=true
        }
    }

    if(!defined(train_vcf_index) && !fix_vcf) {
        call IndexVcf {
            input:
                vcf=train_vcf,
                sv_base_docker=sv_base_docker
        }
    }

    File train_vcf_ = select_first([FixVcf.fixed_vcf, train_vcf])
    File train_vcf_index_ = select_first([FixVcf.fixed_vcf_index, train_vcf_index, IndexVcf.index_file])

    call GetTruthOverlap.GetTruthOverlap as GetTruthOverlap {
        input:
            test_vcfs=[train_vcf_],
            test_vcf_indices=[train_vcf_index_],
            truth_vcfs=truth_vcfs,
            truth_vcf_indices=truth_vcf_indices,
            vapor_files=vapor_files,
            ped_files=[ped_file],
            optimal_overlap_cutoffs=optimal_overlap_cutoffs,
            get_truth_overlap_args=get_truth_overlap_args,
            module03_docker=module03_docker,
            sv_base_docker=sv_base_docker
    }

    call GetTruthOverlap.GetVcfSize as GetVcfSize {
        input:
            vcf=train_vcf_,
            vcf_index=train_vcf_index_,
            sv_base_docker=sv_base_docker
    }

    call TrainGqRecalibratorFilter {
        input:
            train_vcf=train_vcf_,
            train_vcf_index=train_vcf_index_,
            ped_file=ped_file,
            truth_file=GetTruthOverlap.truth_overlap_info,
            genome_tracts=genome_tracts,
            gq_recalibrator_model_file=gq_recalibrator_model_file,
            train_args=train_args,
            gatk_docker=gatk_docker,
            num_entries=GetVcfSize.num_entries
    }

    output {
        File clean_vcf=train_vcf_
        File clean_vcf_index=train_vcf_index_
        File truth_overlap_info = GetTruthOverlap.truth_overlap_info
        File output_optimal_overlap_cutoffs = GetTruthOverlap.output_optimal_overlap_cutoffs
        File output_gq_recalibrator_model_file = TrainGqRecalibratorFilter.output_gq_recalibrator_model_file
    }
}


task FixVcf {
    input {
        File vcf
        String module03_docker
        Boolean index_output_vcf = true
    }

    String fixed_vcf_name = sub(sub(basename(vcf), ".gz$", ""), ".vcf$", "_fixed.vcf.gz")
    String index_file_name = fixed_vcf_name + ".tbi"

    Float uncompress_scale = 10
    Int disk_gb = 1000 + round((2 + uncompress_scale) * size(vcf, "GiB"))
    Int mem_gb_overhead = 2
    Float mem_scale = 2.0
    Float mem_gb = mem_gb_overhead + mem_scale * size(vcf, "GiB")

    runtime {
        docker: module03_docker
        cpu: 1
        preemptible: 3
        max_retries: 0
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail

        module03 fix_vcf ~{vcf} ~{fixed_vcf_name} --index-output-vcf ~{index_output_vcf}
    >>>

    output {
        File fixed_vcf = fixed_vcf_name
        File? fixed_vcf_index = index_file_name
    }
}


task IndexVcf {
    input {
        File vcf
        String sv_base_docker
    }
    parameter_meta {
        vcf: {
          localization_optional: true
        }
    }

    String index_file_name = basename(vcf) + ".tbi"

    Int disk_gb = 10
    Float mem_gb_overhead = 1.0
    Float mem_gb_java = 2.0
    Float mem_gb = mem_gb_java + mem_gb_overhead

    runtime {
        docker: sv_base_docker
        cpu: 1
        preemptible: 3
        max_retries: 0
        memory: mem_gb + " GiB"
        disks: "local-disk " + disk_gb + " HDD"
    }

    command <<<
        set -eu -o pipefail
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

        bcftools index --tbi --threads 2 ~{vcf} -o ~{index_file_name}
    >>>

    output {
        File index_file = index_file_name
    }
}


task TrainGqRecalibratorFilter {
    input {
        File train_vcf
        File train_vcf_index
        File ped_file
        File truth_file
        Array[File] genome_tracts
        File? gq_recalibrator_model_file # can be passed to do extra rounds of training on existing model
        Array[String] train_args = []
        String gatk_docker
        Int? num_entries
        Float mem_scale_vcf_size = 25.2
        Float mem_scale_num_entries = "3.6e-7"
        Float mem_gb_overhead = 1.5
    }

    Int disk_gb = round(1000 + size([train_vcf, train_vcf_index, ped_file, truth_file], "GiB") +
                        size(genome_tracts, "GiB") + size(gq_recalibrator_model_file, "GiB"))
    Float mem_gb_java = if defined(num_entries)
        then 3.0 + mem_scale_num_entries * select_first([num_entries])
        else 3.0 + mem_scale_vcf_size * size(train_vcf, "GiB")
    Float mem_gb = min(mem_gb_java + mem_gb_overhead, 624.0)
    String model_file_name = if defined(gq_recalibrator_model_file)
        then basename(select_first([gq_recalibrator_model_file]))
        else "gq_recalibrator.model"
    String args_str = if length(train_args) > 0 then sep(" ", train_args) else ""

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
        ln -s ~{train_vcf} .
        ln -s ~{train_vcf_index} .
        if ~{defined(gq_recalibrator_model_file)}; then
            cp ~{gq_recalibrator_model_file} .
        fi

        mem_kb_java_actual=$(grep -m1 MemTotal /proc/meminfo \
                             | awk '{printf("%.0f\n", $2 - ~{mem_gb_overhead} * 1048576)}')

        NUM_PHYSICAL_CPUS=$(grep ^"core id" /proc/cpuinfo | sort -u | wc -l)

        gatk --java-options "-Xmx${mem_kb_java_actual}K" XGBoostMinGqVariantFilter \
          --mode "Train" \
          --variant ./$(basename ~{train_vcf}) \
          --pedigree ~{ped_file} \
          --truth-file ~{truth_file} \
          --genome-tract ~{sep=" --genome-tract " genome_tracts} \
          --model-file ~{model_file_name} \
          --n-threads-xgboost $NUM_PHYSICAL_CPUS \
          ~{args_str}

    >>>

    output {
        File output_gq_recalibrator_model_file = model_file_name
    }
}
