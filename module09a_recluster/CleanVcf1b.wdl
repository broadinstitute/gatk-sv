version 1.0

import "Structs.wdl"
import "CleanVcf5.wdl" as CleanVcf5
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CleanVcf1b {
    input {
        File intermediate_vcf
        String prefix
        Int records_per_shard

        String sv_pipeline_docker
        String sv_base_mini_docker
        String sv_pipeline_updates_docker

        RuntimeAttr? runtime_attr_override_subset_large_cnvs
        RuntimeAttr? runtime_attr_override_sort_bed
        RuntimeAttr? runtime_attr_override_intersect_bed
        RuntimeAttr? runtime_attr_override_build_dict
        RuntimeAttr? runtime_attr_override_scatter
        RuntimeAttr? runtime_attr_override_filter_vcf
        RuntimeAttr? runtime_override_concat_vcfs
        RuntimeAttr? runtime_override_cat_multi_cnvs
    }

    call SubsetLargeCNVs {
        input:
            vcf=intermediate_vcf,
            prefix="~{prefix}.subset_large_cnvs",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override_subset_large_cnvs
    }

    call Vcf2Bed {
        input:
            vcf=SubsetLargeCNVs.out,
            prefix="~{prefix}.subset_large_cnvs",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override_subset_large_cnvs
    }

    call SortBed {
        input:
            bed=Vcf2Bed.out,
            prefix="~{prefix}.subset_large_cnvs.sorted",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_override_sort_bed
    }

    call BedtoolsIntersect {
        input:
            bed=SortBed.out,
            prefix="~{prefix}.bedtools_intersect",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_override_intersect_bed
    }

    call BuildGenoNormalReviseDictionary {
        input:
            filtered_vcf=SubsetLargeCNVs.out,
            intersect_bed=BedtoolsIntersect.out,
            prefix="~{prefix}.geno_normal_revise",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override_build_dict
    }

    call MiniTasks.ScatterVcf {
        input:
            vcf=intermediate_vcf,
            records_per_shard=records_per_shard,
            prefix="~{prefix}.scatter_vcf",
            sv_pipeline_docker=sv_pipeline_updates_docker,
            runtime_attr_override=runtime_attr_override_scatter
    }

    scatter ( i in range(length(ScatterVcf.shards)) ) {
        call FilterVcf {
            input:
                intermediate_vcf=ScatterVcf.shards[i],
                dictionary_json_gz=BuildGenoNormalReviseDictionary.out,
                prefix="~{prefix}.filter_vcf.shard_~{i}",
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_override=runtime_attr_override_filter_vcf
        }
    }

    call MiniTasks.ConcatVcfs as ConcatCleanVcf1bShards {
        input:
            vcfs=FilterVcf.out,
            naive=true,
            sort_vcf_list=true,
            outfile_prefix="~{prefix}.concat_vcfs",
            sv_base_mini_docker=sv_pipeline_updates_docker,
            runtime_attr_override=runtime_override_concat_vcfs
    }

    call MiniTasks.CatUncompressedFiles as ConcatMultiCnvs  {
        input:
            shards=FilterVcf.multi_cnvs,
            outfile_name="~{prefix}.multi.cnvs.txt",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_cat_multi_cnvs
    }

    output {
        File normal = ConcatCleanVcf1bShards.concat_vcf
        File multi = ConcatMultiCnvs.outfile
    }
}

task SubsetLargeCNVs {
    input {
        File vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 2.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        bcftools view --no-version \
            -i '(INFO/SVTYPE=="DEL" || INFO/SVTYPE=="DUP") && INFO/SVLEN>=5000' \
            ~{vcf} \
            | bgzip \
            > ~{prefix}.vcf.gz
    >>>
    output {
        File out = "~{prefix}.vcf.gz"
    }
}

task Vcf2Bed {
    input {
        File vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 2.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        svtk vcf2bed --no-header ~{vcf} stdout \
            | awk -F'\t' -v OFS='\t' '{if ($6=="") $6="blanksample";print $0}' \
            | gzip -1 \
            > ~{prefix}.bed.gz
    >>>
    output {
        File out = "~{prefix}.bed.gz"
    }
}

task SortBed {
    input {
        File bed
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 10.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        mkdir tmp
        zcat ~{bed} \
            | sort -T tmp -k1,1 -k2,2n \
            | gzip -1 \
            > ~{prefix}.bed.gz
    >>>
    output {
        File out = "~{prefix}.bed.gz"
    }
}

task BedtoolsIntersect {
    input {
        File bed
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(bed, "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 10.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        bedtools intersect -sorted -wa -wb -a <(zcat ~{bed}) -b <(zcat ~{bed}) \
            | awk -F'\t' -v OFS='\t' '$4!=$10 && $5!=$11' \
            | gzip -1 \
            > ~{prefix}.bed.gz
    >>>
    output {
        File out = "~{prefix}.bed.gz"
    }
}

task BuildGenoNormalReviseDictionary {
    input {
        File filtered_vcf
        File intersect_bed
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([filtered_vcf, intersect_bed], "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 2.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        python /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1b_build_dict.py ~{filtered_vcf} ~{intersect_bed} \
            | gzip -1 \
            > ~{prefix}.json.gz
    >>>
    output {
        File out = "~{prefix}.json.gz"
    }
}

task FilterVcf {
    input {
        File intermediate_vcf
        File dictionary_json_gz
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size([intermediate_vcf, dictionary_json_gz], "GB")
    RuntimeAttr runtime_default = object {
                                      mem_gb: 3.75,
                                      disk_gb: ceil(10.0 + input_size * 2.0),
                                      cpu_cores: 1,
                                      preemptible_tries: 3,
                                      max_retries: 1,
                                      boot_disk_gb: 10
                                  }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        python /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part1b_filter.py ~{dictionary_json_gz} ~{intermediate_vcf} \
            | bgzip \
            > ~{prefix}.vcf.gz
        mv multi.cnvs.txt ~{prefix}.multi.cnvs.txt
    >>>
    output {
        File out = "~{prefix}.vcf.gz"
        File multi_cnvs = "~{prefix}.multi.cnvs.txt"
    }
}
