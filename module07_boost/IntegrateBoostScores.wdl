##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

## Copyright Broad Institute, 2020
## 
## This WDL pipeline implements Duphold 
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

version 1.0

import "Structs.wdl"

workflow IntegrateBoostResultsAcrossBatches{
    input{
        Array[File] boost_models
        File boost_cutoff_table
        File contig_list
        String prefix
        String sv_benchmark_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_boost_integration_across_batches
        RuntimeAttr? runtime_attr_split_boost_stat_by_contig
        RuntimeAttr? runtime_attr_concat_stat
    }

    Array[String] contigs = transpose(read_tsv(contig_list))[0]

    scatter (boost_model in boost_models){
        call IntegrateBoostScoreAcrossBatches{
            input:
                boost_model = boost_model,
                boost_cutoff_table = boost_cutoff_table,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_boost_integration_across_batches
        }
    }

    scatter (contig in contigs){
        call SplitBoostStatByContig{
            input:
                boost_stats = IntegrateBoostScoreAcrossBatches.boost_stat
                contig = contig,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_split_boost_stat_by_contig
        }
    }

    call ConcatStats{
        input:
            shard_stats = SplitBoostStatByContig.per_contig_stat,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_stat
    }
}

Task IntegrateBoostScoreAcrossBatches{
    input{
        File boost_model
        File boost_cutoff_table
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File boost_stat = "~{boostname}.stat"
    }

    String boostname = basename(BoostModel, ".tar.gz")
    
    command <<<
        set -Eeuo pipefail

        mkdir boost_models/
        tar zxvf ~{BoostModel} -C boost_models/
        python /src/Integrate_Boost_Score_across_Samples.py \
            --cff_table ~{BoostCutoffTable} \
            --path "boost_models/~{boostname}/" \
            --output "~{boostname}.stat"
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

Task SplitBoostStatByContig{
    input{
        Array[File] boost_stats
        String contig
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File per_contig_stat = "~{boost_stat_name}.~{contig}.stat.gz"
    }

    String boost_stat_name = basename(BoostStat,".stat")

    command <<<
        set -Eeuo pipefail

        while read BoostStat; do
            sed -e 's/\./\t/g' ~{BoostStat} | awk '{if ($2=="~{contig}") print}' > ~{contig}.~{BoostStat}
        done < ~{write_lines(boost_stats)}

        Rscript Integrate_Boost_Score_across_Batches.R \
            --prefix ~{boost_stat_name} \
            --chr ~{contig}

        bgzip ~{boost_stat_name}.~{contig}.stat
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ConcatStats {
    input {
        Array[File] shard_stats
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
    # be held in memory or disk while working, potentially in a form that takes up more space)
    Float input_size = size(shard_stats, "GB")
    Float compression_factor = 5.0
    Float base_disk_gb = 5.0
    Float base_mem_gb = 2.0
    RuntimeAttr runtime_default = object {
        mem_gb: base_mem_gb + compression_factor * input_size,
        disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
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
        set -eu
        head -1 shard_stats[0] > ~{prefix}.stat

        set -o pipefail
        while read SPLIT; do
            zcat $SPLIT | tail -n+2 >> ~{prefix}.stat 
        done < ~{write_lines(shard_stats)} 

        bgzip ~{prefix}.stat 
    >>>

    output {
        File merged_stat = "~{prefix}.stat.gz"
    }
}


