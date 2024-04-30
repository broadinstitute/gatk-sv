version 1.0
    
import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks


workflow GetBatchedVcf {

    input {
        File vcf_file
        String prefix
        File samples
        Int records_per_shard
        String variant_interpretation_docker
        String sv_pipeline_updates_docker
        RuntimeAttr? runtime_attr_batch_vcf
        RuntimeAttr? runtime_override_shard_vcf
        RuntimeAttr? runtime_attr_merge
    }

    call MiniTasks.ScatterVcf as SplitVcf {
        input:
            vcf=vcf_file,
            prefix=prefix,
            records_per_shard=records_per_shard,
            sv_pipeline_docker=sv_pipeline_updates_docker,
            runtime_attr_override=runtime_override_shard_vcf
    }

    scatter (shard in SplitVcf.shards) {
        call BatchVcf {
            input:
                vcf_file = shard,
                samples = samples,
                variant_interpretation_docker=variant_interpretation_docker,
                runtime_attr_override = runtime_attr_batch_vcf
        }
    }

    call Merge {
        input:
            vcf_files = BatchVcf.subset_vcf,
            variant_interpretation_docker=variant_interpretation_docker,
            runtime_attr_override = runtime_attr_merge
    }

    output {
        File split_vcf = Merge.merged_vcf
    }

}

task BatchVcf {

    input {
        File vcf_file
        File samples
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(samples, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(10 + input_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File subset_vcf = "filtered.vcf.gz"
    }

    command {
        set -exuo pipefail
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools view -I -S ~{samples} ~{vcf_file} -O z -o filtered.vcf.gz
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}

task Merge {

    input {
        Array[File] vcf_files
        String variant_interpretation_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf_files, "GB")
    Float base_mem_gb = 3.75

    RuntimeAttr default_attr = object {
        mem_gb: base_mem_gb,
        disk_gb: ceil(10 + input_size * 1.5),
        cpu_cores: 1,
        preemptible_tries: 2,
        max_retries: 1,
        boot_disk_gb: 8
    }
    
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File merged_vcf = "Merged.vcf.gz"
    }

    command {
        set -exuo pipefail

        for vcf in ~{sep=' ' vcf_files}
        do
            export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            tabix -p vcf $vcf
        done

        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        bcftools concat ~{sep=' ' vcf_files} -O z -o Merged.vcf.gz
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
        docker: variant_interpretation_docker
    }
}
