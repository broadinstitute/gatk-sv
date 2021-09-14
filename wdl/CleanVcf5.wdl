version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow CleanVcf5 {
    input {
        File normal_revise_vcf
        File revise_vcf_lines
        File ped_file
        File sex_chr_revise
        File multi_ids
        File? outlier_samples_list

        File? shard_script
        File? make_clean_gq_script
        File? find_redundant_sites_script

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_override
    }

    call ScatterVcf as ScatterVcfPart5 {
        input:
            vcf=normal_revise_vcf,
            records_per_shard = 5000,
            shard_script=shard_script,
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override
    }

    scatter ( clean_gq_shard in ScatterVcfPart5.shards ) {
        call CleanVcf5MakeCleanGQ {
            input:
                revise_vcf_lines=revise_vcf_lines,
                normal_revise_vcf=clean_gq_shard,
                ped_file=ped_file,
                sex_chr_revise=sex_chr_revise,
                multi_ids=multi_ids,
                outlier_samples_list=outlier_samples_list,
                make_clean_gq_script=make_clean_gq_script,
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_override=runtime_attr_override
        }
    }

    call MiniTasks.ConcatVcfs as GatherMultiallelicShards {
        input:
            vcfs = CleanVcf5MakeCleanGQ.multiallelic_vcf,
            outfile_prefix="multiallelic",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_override
    }

    call CleanVcf5FindRedundantMultiallelics {
        input:
            multiallelic_vcf=GatherMultiallelicShards.concat_vcf,
            find_redundant_sites_script=find_redundant_sites_script,
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override
    }

    call CleanVcf5Polish {
        input:
            clean_gq_vcfs=CleanVcf5MakeCleanGQ.clean_gq_vcf,
            no_sample_lists=CleanVcf5MakeCleanGQ.no_sample_list,
            redundant_multiallelics_list=CleanVcf5FindRedundantMultiallelics.redundant_multiallelics_list,
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override
    }

    output {
        File polished=CleanVcf5Polish.polished
    }
}

task ScatterVcf {
    input {
        File vcf
        File? shard_script
        String prefix = "shard"
        Int? records_per_shard = 5000
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 10.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
        set -eu -o pipefail
        python3 ~{default="/opt/sv-pipeline/04-variant-resolution/scripts/shard_pysam.py" shard_script} ~{vcf} ~{prefix} ~{records_per_shard}
    >>>
    output {
        Array[File] shards = glob("shard.*.vcf.gz")
    }
}

task CleanVcf5MakeCleanGQ {
    input {
        File revise_vcf_lines
        File normal_revise_vcf
        File ped_file
        File sex_chr_revise
        File multi_ids
        File? outlier_samples_list
        File? make_clean_gq_script
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    # generally assume working memory is ~3 * inputs
    Float input_size = size(
                       select_all([revise_vcf_lines, normal_revise_vcf, ped_file, sex_chr_revise, multi_ids, outlier_samples_list]),
                       "GB")
    Float base_disk_gb = 10.0
    Float input_mem_scale = 3.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
        set -eu -o pipefail

        ~{if defined(outlier_samples_list) then "ln ~{outlier_samples_list} outliers.txt" else "touch outliers.txt"}

        # put the revise lines into a normal VCF format
        bcftools view -h ~{normal_revise_vcf} > header.txt
        cat header.txt <(zcat ~{revise_vcf_lines} | grep . | tr " " "\t") | bgzip -c > revise.vcf.lines.vcf.gz

        python3 ~{default="/opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part5_new_update_records.py" make_clean_gq_script} \
            revise.vcf.lines.vcf.gz \
            ~{normal_revise_vcf} \
            ~{ped_file} \
            ~{sex_chr_revise} \
            ~{multi_ids} \
            outliers.txt

        tabix cleanGQ.vcf.gz
    >>>

    output {
        File clean_gq_vcf="cleanGQ.vcf.gz"
        File clean_gq_vcf_idx="cleanGQ.vcf.gz.tbi"
        File multiallelic_vcf="multiallelic.vcf.gz"
        File no_sample_list = "no_called_samples.list"
    }
}

task CleanVcf5FindRedundantMultiallelics {
    input {
        File multiallelic_vcf
        File? find_redundant_sites_script
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    # generally assume working memory is ~3 * inputs
    Float input_size = size(multiallelic_vcf, "GB")
    Float base_disk_gb = 10.0
    Float input_mem_scale = 4.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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

        python3 ~{default="/opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part5_find_redundant_multiallelics.py" find_redundant_sites_script} \
            ~{multiallelic_vcf} \
            redundant_multiallelics.list

    >>>

    output {
        File redundant_multiallelics_list="redundant_multiallelics.list"
    }
}


task CleanVcf5Polish {
    input {
        Array[File] clean_gq_vcfs
        Array[File] no_sample_lists
        File redundant_multiallelics_list
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    # generally assume working memory is ~3 * inputs
    Float input_size = size(clean_gq_vcfs, "GB")
    Float base_disk_gb = 10.0
    Float input_mem_scale = 4.0
    Float input_disk_scale = 5.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
                                      cpu_cores: 4,
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

        cat ~{redundant_multiallelics_list} ~{sep=" " no_sample_lists} > ids_to_remove.list
        /usr/local/bin/bcftools concat --no-version --output-type u --file-list ~{write_lines(clean_gq_vcfs)} | \
            /usr/local/bin/bcftools view --no-version \
                --exclude 'ID=@ids_to_remove.list' \
                --output-type z -o polished.need_reheader.vcf.gz --threads 3

        # do the last bit of header cleanup
        bcftools view -h polished.need_reheader.vcf.gz | awk 'NR == 1' > new_header.vcf
        ls -l new_header.vcf
        bcftools view -h polished.need_reheader.vcf.gz \
            | awk 'NR > 1' \
            | egrep -v "CIPOS|CIEND|RMSSTD|EVENT|INFO=<ID=UNRESOLVED,|source|varGQ|bcftools|ALT=<ID=UNR|INFO=<ID=MULTIALLELIC," \
            | sort -k1,1 >> new_header.vcf
        ls -l new_header.vcf
        bcftools reheader polished.need_reheader.vcf.gz -h new_header.vcf -o polished.vcf.gz
        ls -l polished.vcf.gz
    >>>

    output {
        File polished="polished.vcf.gz"
    }
}
