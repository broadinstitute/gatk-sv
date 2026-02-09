version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as util

workflow ShardedAnnotateVcf {
    input {
        File vcf
        File vcf_idx
        String contig
        String prefix

        File lps_tsv
        File sample_pop_assignments
        File ped_file
        File par_bed
        Array[String]? strip_info_fields

        Int records_per_shard
        String sv_pipeline_docker
        String sv_base_mini_docker
        String gatk_docker

        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_strip_info_fields
        RuntimeAttr? runtime_attr_compute_AFs
    }

    call MiniTasks.ScatterVcf {
        input:
            vcf = vcf,
            vcf_index = vcf_idx,
            records_per_shard = records_per_shard,
            contig = contig,
            prefix = prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_scatter_vcf
    }

    scatter (i in range(length(ScatterVcf.shards))) {
        if (defined(strip_info_fields)) {
            call StripInfoFields {
                input:
                    vcf = ScatterVcf.shards[i],
                    vcf_index = ScatterVcf.shards_index[i],
                    info_fields = select_first([strip_info_fields]),
                    prefix = "~{prefix}.~{i}.stripped",
                    docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_strip_info_fields
            }
        }
        
        call ComputeAFs {
            input:
                vcf = select_first([StripInfoFields.stripped_vcf, ScatterVcf.shards[i]]),
                lps_tsv = lps_tsv,
                sample_pop_assignments = sample_pop_assignments,
                ped_file = ped_file,
                par_bed = par_bed,
                prefix = "~{prefix}.~{i}.AF",
                docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_compute_AFs
        }
    }

    output {
        Array[File] sharded_annotated_vcf = ComputeAFs.af_vcf
        Array[File] sharded_annotated_vcf_idx = ComputeAFs.af_vcf_idx
    }
}

task StripInfoFields {
    input {
        File vcf
        File vcf_index
        String prefix
        Array[String] info_fields
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x INFO/~{sep=",INFO/" info_fields} \
            -Oz -o ~{prefix}.vcf.gz \
            ~{vcf}
      
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File stripped_vcf = "~{prefix}.vcf.gz"
        File stripped_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }
    
    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ComputeAFs {
    input {
        File vcf
        File sample_pop_assignments
        File ped_file
        File par_bed
        File lps_tsv
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py "~{vcf}" stdout \
            ~{"-p " + sample_pop_assignments} \
            ~{"-f " + ped_file} \
            ~{"--par " + par_bed} \
            ~{"-l " + lps_tsv} \
            | bgzip -c > "~{prefix}.vcf.gz"

        tabix -p vcf "~{prefix}.vcf.gz"
    >>>
    
    output {
        File af_vcf = "~{prefix}.vcf.gz"
        File af_vcf_idx = "~{prefix}.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(vcf, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 3,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
