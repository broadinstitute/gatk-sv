version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks

workflow CleanVcf5 {
    input {
        File normal_revise_vcf
        File revise_vcf_lines
        File ped_file
        File sex_chr_revise
        File multi_ids
        File? outlier_samples_list

        String prefix
        String contig
        Int records_per_shard

        File? make_clean_gq_script
        File? find_redundant_sites_script

        String sv_base_mini_docker
        String sv_pipeline_docker

        Int? threads_per_task
        RuntimeAttr? runtime_attr_override_scatter
        RuntimeAttr? runtime_attr_override_make_cleangq
        RuntimeAttr? runtime_attr_override_find_redundant_multiallelics
        RuntimeAttr? runtime_attr_override_polish
    }

    call tasks.ScatterVcf {
        input:
            vcf=normal_revise_vcf,
            records_per_shard = records_per_shard,
            prefix = "~{prefix}.scatter_vcf",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override_scatter
    }

    scatter ( i in range(length(ScatterVcf.shards)) ) {
        call MakeCleanGQ {
            input:
                revise_vcf_lines=revise_vcf_lines,
                normal_revise_vcf=ScatterVcf.shards[i],
                ped_file=ped_file,
                sex_chr_revise=sex_chr_revise,
                multi_ids=multi_ids,
                outlier_samples_list=outlier_samples_list,
                make_clean_gq_script=make_clean_gq_script,
                prefix="~{prefix}.make_clean_gq.shard_~{i}",
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_override=runtime_attr_override_make_cleangq
        }
    }

    call FindRedundantMultiallelics {
        input:
            multiallelic_vcfs=MakeCleanGQ.multiallelic_vcf,
            find_redundant_sites_script=find_redundant_sites_script,
            prefix="~{prefix}.find_redundant_multiallelics",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override_find_redundant_multiallelics
    }

    call Polish {
        input:
            clean_gq_vcfs=MakeCleanGQ.clean_gq_vcf,
            no_sample_lists=MakeCleanGQ.no_sample_list,
            redundant_multiallelics_list=FindRedundantMultiallelics.redundant_multiallelics_list,
            prefix="~{prefix}.polish",
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_attr_override_polish
    }

    output {
        File polished=Polish.polished
    }
}

task MakeCleanGQ {
    input {
        File revise_vcf_lines
        File normal_revise_vcf
        File ped_file
        File sex_chr_revise
        File multi_ids
        File? outlier_samples_list
        File? make_clean_gq_script
        String prefix
        Int? threads = 2
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    # generally assume working memory is ~3 * inputs
    Float input_size = size(
                       select_all([revise_vcf_lines, normal_revise_vcf, ped_file, sex_chr_revise, multi_ids, outlier_samples_list]),
                       "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
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

        python3 ~{default="/opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part5_update_records.py" make_clean_gq_script} \
            --threads_per_file ~{threads} \
            revise.vcf.lines.vcf.gz \
            ~{normal_revise_vcf} \
            ~{ped_file} \
            ~{sex_chr_revise} \
            ~{multi_ids} \
            outliers.txt \
            ~{prefix}

        bcftools view -G -O z ~{prefix}.multiallelic.vcf.gz > ~{prefix}.multiallelic.sites.vcf.gz
        tabix ~{prefix}.cleanGQ.vcf.gz
    >>>

    output {
        File clean_gq_vcf=prefix + ".cleanGQ.vcf.gz"
        File clean_gq_vcf_idx=prefix + ".cleanGQ.vcf.gz.tbi"
        File multiallelic_vcf=prefix + ".multiallelic.sites.vcf.gz"
        File no_sample_list = prefix + ".no_called_samples.list"
    }
}

task FindRedundantMultiallelics {
    input {
        Array[File] multiallelic_vcfs
        File? find_redundant_sites_script
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    # generally assume working memory is ~3 * inputs
    Float input_size = size(multiallelic_vcfs, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
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
        VCFS="~{write_lines(multiallelic_vcfs)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        bcftools concat --no-version --output-type z --file-list vcfs_sorted.list --output multiallelic.vcf.gz

        python3 ~{default="/opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part5_find_redundant_multiallelics.py" find_redundant_sites_script} \
            multiallelic.vcf.gz \
            ~{prefix}.list

    >>>

    output {
        File redundant_multiallelics_list="~{prefix}.list"
    }
}


task Polish {
    input {
        Array[File] clean_gq_vcfs
        Array[File] no_sample_lists
        File redundant_multiallelics_list
        String prefix
        String sv_pipeline_docker
        Int threads = 2
        RuntimeAttr? runtime_attr_override
    }

    # generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
    # generally assume working memory is ~3 * inputs
    Float input_size = size(clean_gq_vcfs, "GB")
    Float base_disk_gb = 10.0

    RuntimeAttr runtime_default = object {
                                      mem_gb: 16,
                                      disk_gb: ceil(base_disk_gb + input_size * 5.0),
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

        VCFS="~{write_lines(clean_gq_vcfs)}"
        cat $VCFS | awk -F '/' '{print $NF"\t"$0}' | sort -k1,1V | awk '{print $2}' > vcfs_sorted.list
        cat ~{redundant_multiallelics_list} ~{sep=" " no_sample_lists} > ids_to_remove.list
        bcftools concat --no-version --output-type u --file-list vcfs_sorted.list | \
            bcftools view --no-version \
                --exclude 'ID=@ids_to_remove.list' \
                --output-type z -o polished.need_reheader.vcf.gz --threads ~{threads}

        # do the last bit of header cleanup
        bcftools view -h polished.need_reheader.vcf.gz > original_header.vcf
        cat original_header.vcf | fgrep '##fileformat' > new_header.vcf
        cat original_header.vcf \
            | egrep -v "CIPOS|CIEND|RMSSTD|EVENT|INFO=<ID=UNRESOLVED,|source|varGQ|bcftools|ALT=<ID=UNR|INFO=<ID=MULTIALLELIC|GATKCommandLine|#CHROM|##contig|##fileformat" \
            | sort >> new_header.vcf
        # Don't sort contigs lexicographically, which would result in incorrect chr1, chr10, chr11, ... ordering
        cat original_header.vcf | fgrep '##contig' >> new_header.vcf
        cat original_header.vcf | fgrep '#CHROM' >> new_header.vcf
        bcftools reheader polished.need_reheader.vcf.gz -h new_header.vcf -o ~{prefix}.vcf.gz
    >>>

    output {
        File polished="~{prefix}.vcf.gz"
    }
}
