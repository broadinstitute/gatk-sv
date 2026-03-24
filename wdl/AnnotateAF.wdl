version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow AnnotateAF {
    input {
        File vcf
        File vcf_index
        Array[String] contigs
        String prefix

        File sample_pop_assignments
        File ped_file
        File par_bed
        File lps_tsv
        Array[String]? strip_info_fields

        Int records_per_shard
        String gatk_sv_lr_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_subset_tsv
        RuntimeAttr? runtime_attr_split_ref_bed
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_strip_info_fields
        RuntimeAttr? runtime_attr_compute_AFs
        RuntimeAttr? runtime_attr_concat
    }

    scatter (contig in contigs) {
        call SubsetLpsTsvToContig {
            input:
                tsv = lps_tsv,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                docker = gatk_sv_lr_docker,
                runtime_attr_override = runtime_attr_subset_tsv
        }

        call sharded_annotate_vcf.ShardedAnnotateVcf {
            input:
                vcf = vcf,
                vcf_idx = vcf_index,
                contig = contig,
                prefix = "~{prefix}.~{contig}",
                sample_pop_assignments = sample_pop_assignments,
                ped_file = ped_file,
                par_bed = par_bed,
                lps_tsv = SubsetLpsTsvToContig.subset_tsv,
                strip_info_fields = strip_info_fields,
                records_per_shard = records_per_shard,
                sv_pipeline_docker = gatk_sv_lr_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
                runtime_attr_strip_info_fields = runtime_attr_strip_info_fields,
                runtime_attr_compute_AFs = runtime_attr_compute_AFs
        }
    }

    call MiniTasks.ConcatVcfs {
        input:
            vcfs = flatten(ShardedAnnotateVcf.sharded_annotated_vcf),
            vcfs_idx = flatten(ShardedAnnotateVcf.sharded_annotated_vcf_idx),
            allow_overlaps = true,
            outfile_prefix = "~{prefix}.annotated",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat
    }

    output {
        File af_annotated_vcf = ConcatVcfs.concat_vcf
        File af_annotated_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}

task SubsetLpsTsvToContig {
    input {
        File tsv
        String contig
        String prefix
        String docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -eou pipefail

        zcat -f ~{tsv} \
            | awk -F'\t' -v contig="~{contig}" '
                BEGIN {
                    if (contig ~ /^chr/) contig = substr(contig, 4)
                }

                {
                    split($1, arr, ",")
                    split(arr[1], parts, "-")
                    chr = parts[1]
                    if (chr ~ /^chr/) chr = substr(chr, 4)
                    if (chr ~ /^trid/) header = $0
                    else if (chr == contig) data[++n] = $0
                }

                END {
                    print header
                    for (i = 1; i <= n; i++) print data[i]
                }
            ' \
            > ~{prefix}.tsv
    >>>

    output {
        File subset_tsv = "~{prefix}.tsv"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 4,
        disk_gb: 2 * ceil(size(tsv, "GB")) + 5,
        boot_disk_gb: 10,
        preemptible_tries: 2,
        max_retries: 0
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
