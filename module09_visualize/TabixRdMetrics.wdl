#this script is developed to fix false postiive inversions in gnomad v3 vcf that are >1Mb in size but have PE_GT=0
#INVs with no passing samples have filter column revised from "PASS" to "UNRESOLVED"

#developing workdir on erisone: /data/talkowski/xuefang/data/gnomad_V3/module08/step9_sm_depyh_only_dup_fix

version 1.0

import "Structs.wdl"

workflow TabixRdMetrics {

    input {
        Array[File] bed_metrics
        Array[File] bed_metrics_idxes
        Array[File] region_bed_list
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_extract_region_metrics
    }

    scatter (i in range(length(bed_metrics))){

        call ExtractRegionMetrics{
            input:
                bed_metric = bed_metrics[i],
                bed_metric_idx = bed_metrics_idxes[i],
                bed = region_bed_list[i],
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_extract_region_metrics
        }
    }
    output{
        Array[File] targeted_metrics = ExtractRegionMetrics.targeted_metric
        Array[File] targeted_metrics_idxes = ExtractRegionMetrics.targeted_metric_idx
    }
}



task ExtractRegionMetrics{
    input{
        File bed_metric
        File bed_metric_idx
        File bed

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 15, 
        disk_gb: ceil(40.0 + size(bed_metric, "GiB") * 4),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(bed_metric,".txt.gz")

    command<<<
        set -euo pipefail
        sort -k1,1 -k2,2n ~{bed} | bedtools merge -i - > sorted.bed
        tabix -h -R sorted.bed ~{bed_metric} | bgzip -c > ~{prefix}.targeted.txt.gz
        tabix -b 2 -e 2 ~{prefix}.targeted.txt.gz
    >>>

    output{
        File targeted_metric = "~{prefix}.targeted.txt.gz"
        File targeted_metric_idx = "~{prefix}.targeted.txt.gz.tbi"
    }

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


