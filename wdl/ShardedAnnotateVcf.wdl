version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as util
import "AnnotateFunctionalConsequences.wdl" as func
import "AnnotateExternalAFPerShard.wdl" as eaf

workflow ShardedAnnotateVcf {
    input {
        File vcf
        File vcf_idx
        String prefix
        String contig

        File? protein_coding_gtf
        File? noncoding_bed
        Int? promoter_window
        Int? max_breakend_as_cnv_length
        String? svannotate_additional_args
        Array[String]? strip_info_fields

        Boolean annotate_external_af
        Boolean annotate_internal_af
        Boolean annotate_functional_consequences

        File? sample_pop_assignments # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
        File? sample_keep_list
        File? ped_file                # Used for M/F AF calculations
        File? par_bed
        File? allosomes_list
        Int   sv_per_shard

        File? ref_bed              # File with external allele frequencies
        String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
        Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

        String sv_pipeline_docker
        String sv_base_mini_docker
        String gatk_docker

        RuntimeAttr? runtime_attr_split_ref_bed
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_subset_vcf_by_samples_list
        RuntimeAttr? runtime_attr_strip_info_fields
        RuntimeAttr? runtime_attr_svannotate
        RuntimeAttr? runtime_attr_compute_AFs
        RuntimeAttr? runtime_attr_modify_vcf
        RuntimeAttr? runtime_attr_split_query_vcf
        RuntimeAttr? runtime_attr_bedtools_closest
        RuntimeAttr? runtime_attr_select_matched_svs
    }

    if (annotate_external_af && defined(ref_bed)) {
        call eaf.SplitRefBed {
            input:
            bed = select_first([ref_bed]),
            contig = contig,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_split_ref_bed
        }
    }

    call MiniTasks.ScatterVcf {
        input:
            vcf = vcf,
            vcf_index = vcf_idx,
            prefix = prefix,
            records_per_shard = sv_per_shard,
            contig = contig,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_scatter_vcf
    }

    scatter (i in range(length(ScatterVcf.shards))) {
        String shard_prefix = "~{prefix}.~{contig}.~{i}"

        if (defined(sample_keep_list)) {
            call util.SubsetVcfBySamplesList {
                input:
                    vcf = ScatterVcf.shards[i],
                    list_of_samples = select_first([sample_keep_list]),
                    sv_base_mini_docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_subset_vcf_by_samples_list
            }
        }

        File vcf_initial = select_first([SubsetVcfBySamplesList.vcf_subset, ScatterVcf.shards[i]])
        File vcf_idx_initial = select_first([SubsetVcfBySamplesList.vcf_subset_index, ScatterVcf.shards_index[i]])
        if (defined(strip_info_fields)) {
            call StripInfoFields {
                input:
                    vcf = vcf_initial,
                    vcf_index = vcf_idx_initial,
                    prefix = shard_prefix,
                    info_fields = select_first([strip_info_fields]),
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_strip_info_fields
            }
        }

        File vcf_for_chain = select_first([StripInfoFields.out_vcf, vcf_initial])
        File vcf_idx_for_chain = select_first([StripInfoFields.out_vcf_idx, vcf_idx_initial])
        if (annotate_functional_consequences) {
            call func.AnnotateFunctionalConsequences {
                input:
                    vcf = vcf_for_chain,
                    vcf_index = vcf_idx_for_chain,
                    prefix = shard_prefix,
                    protein_coding_gtf = protein_coding_gtf,
                    noncoding_bed = noncoding_bed,
                    promoter_window = promoter_window,
                    max_breakend_as_cnv_length = max_breakend_as_cnv_length,
                    additional_args = svannotate_additional_args,
                    gatk_docker = gatk_docker,
                    runtime_attr_svannotate = runtime_attr_svannotate
            }
        }

        File vcf_after_func = select_first([AnnotateFunctionalConsequences.annotated_vcf, vcf_for_chain])
        File vcf_idx_after_func = select_first([AnnotateFunctionalConsequences.annotated_vcf_index, vcf_idx_for_chain])
        if (annotate_internal_af) {
            call ComputeAFs {
                input:
                    vcf = vcf_after_func,
                    prefix = shard_prefix,
                    sample_pop_assignments = sample_pop_assignments,
                    ped_file = ped_file,
                    par_bed = par_bed,
                    allosomes_list = allosomes_list,
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_compute_AFs
            }
        }

        File vcf_after_af = select_first([ComputeAFs.af_vcf, vcf_after_func])
        File vcf_idx_after_af = select_first([ComputeAFs.af_vcf_idx, vcf_idx_after_func])
        if (annotate_external_af && defined(ref_bed)) {
            call eaf.AnnotateExternalAFPerShard {
                input:
                    vcf = vcf_after_af,
                    vcf_idx = vcf_idx_after_af,
                    split_ref_bed_del = select_first([SplitRefBed.del]),
                    split_ref_bed_dup = select_first([SplitRefBed.dup]),
                    split_ref_bed_ins = select_first([SplitRefBed.ins]),
                    split_ref_bed_inv = select_first([SplitRefBed.inv]),
                    split_ref_bed_bnd = select_first([SplitRefBed.bnd]),
                    population = select_first([population]),
                    ref_prefix = select_first([ref_prefix]),
                    prefix = "~{prefix}.~{contig}.~{i}",
                    sv_base_mini_docker = sv_base_mini_docker,
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_modify_vcf = runtime_attr_modify_vcf,
                    runtime_attr_split_query_vcf = runtime_attr_split_query_vcf,
                    runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
                    runtime_attr_select_matched_svs = runtime_attr_select_matched_svs
            }
        }
        File final_vcf = select_first([AnnotateExternalAFPerShard.annotated_vcf, vcf_after_af])
        File final_vcf_idx = select_first([AnnotateExternalAFPerShard.annotated_vcf_tbi, vcf_idx_after_af])
    }

    output {
        Array[File] sharded_annotated_vcf = final_vcf
        Array[File] sharded_annotated_vcf_idx = final_vcf_idx
    }
}

task StripInfoFields {
    input {
        File vcf
        File vcf_index
        String prefix
        Array[String] info_fields
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        bcftools annotate \
            -x INFO/~{sep=",INFO/" info_fields} \
            -Oz -o ~{prefix}.stripped.vcf.gz \
            ~{vcf}
      
        tabix -p vcf ~{prefix}.stripped.vcf.gz
    >>>

    output {
        File out_vcf = "~{prefix}.stripped.vcf.gz"
        File out_vcf_idx = "~{prefix}.stripped.vcf.gz.tbi"
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
  }

task ComputeAFs {
    input {
        File vcf
        String prefix
        File? sample_pop_assignments
        File? ped_file
        File? par_bed
        File? allosomes_list
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euo pipefail

        /opt/sv-pipeline/05_annotation/scripts/compute_AFs.py "~{vcf}" stdout \
            ~{"-p " + sample_pop_assignments} \
            ~{"-f " + ped_file} \
            ~{"--par " + par_bed} \
            ~{"--allosomes-list " + allosomes_list} \
            | bgzip -c \
            > "~{prefix}.wAFs.vcf.gz"

        tabix -p vcf "~{prefix}.wAFs.vcf.gz"
    >>>
    
    output {
        File af_vcf = "~{prefix}.wAFs.vcf.gz"
        File af_vcf_idx = "~{prefix}.wAFs.vcf.gz.tbi"
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
