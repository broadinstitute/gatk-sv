version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "AnnotateGenomicContext.wdl" as AnnotateGenomicContext
import "TasksBenchmark.wdl" as mini_tasks

# Workflow to annotate vcf file with genomic context
workflow ReClusterCleanVcfSVs {
    input {
        File vcf
        File vcf_index
        File Repeat_Masks
        File Simple_Repeats
        File Segmental_Duplicates
        File? annotated_SV_origin

        String prefix
        String vid_prefix
        String svtype
        String Genomic_Context
        Int min_size
        Int max_size

        Int num_vids
        Int num_samples
        Int dist
        Int svsize
        Float frac
        Float sample_overlap

        String sv_base_mini_docker
        String sv_benchmark_docker
        String sv_pipeline_docker

        # overrides for MiniTasks
        RuntimeAttr? runtime_override_vcf_to_bed
        RuntimeAttr? runtime_attr_override_anno_gc
        RuntimeAttr? runtime_attr_override_inte_gc
        RuntimeAttr? runtime_override_extract_SV_sites
        RuntimeAttr? runtime_attr_override_svtk_vcfcluster
        RuntimeAttr? runtime_attr_override_extract_svid
        RuntimeAttr? runtime_attr_override_integrate_vcfs
        RuntimeAttr? runtime_attr_override_sort_reclustered_vcf
    }

    if(!defined(annotated_SV_origin)){
        call AnnotateGenomicContext.AnnotateSVsWithGenomicContext as AnnotateSVsWithGenomicContext {
            input:
                vcf = vcf,
                vcf_index = vcf_index,
                Repeat_Masks = Repeat_Masks,
                Simple_Repeats = Simple_Repeats,
                Segmental_Duplicates = Segmental_Duplicates,

                sv_base_mini_docker  =sv_base_mini_docker,
                sv_benchmark_docker = sv_benchmark_docker,
                sv_pipeline_docker = sv_pipeline_docker,

                runtime_override_extract_SV_sites = runtime_override_extract_SV_sites,
                runtime_override_vcf_to_bed = runtime_override_vcf_to_bed,
                runtime_attr_override_anno_gc = runtime_attr_override_anno_gc,
                runtime_attr_override_inte_gc = runtime_attr_override_inte_gc
        }
    }

    File annotated_SVs = select_first([annotated_SV_origin, AnnotateSVsWithGenomicContext.annotated_SVs])

    call Extract_SVID_by_GC {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            SVID_GC = annotated_SVs,
            Genomic_Context = Genomic_Context,
            svtype = svtype,
            min_size = min_size,
            max_size = max_size,

            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_extract_svid
    }

    call SvtkVcfCluster {
        input:
            vcf = Extract_SVID_by_GC.out_vcf,
            vcf_idx = Extract_SVID_by_GC.out_idx,
            prefix = prefix,
            vid_prefix = vid_prefix,
            num_vids = num_vids,
            num_samples = num_samples,
            dist = dist,
            frac = frac,
            sample_overlap = sample_overlap,
            svsize = svsize,
            svtype = svtype, 
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_override_svtk_vcfcluster
    }

    output{
        File reclustered_SV = SvtkVcfCluster.vcf_out
        File reclustered_SV_idx = SvtkVcfCluster.vcf_out_idx
    }
}



task SvtkVcfCluster {
    input {
        File vcf
        File vcf_idx
        String prefix
        String vid_prefix
        Int num_vids
        Int num_samples
        Int dist
        Float frac
        Float sample_overlap
        File? exclude_list
        File? exclude_list_idx
        Int svsize
        String svtype
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float default_mem_gb = 20 + (200.0 * (num_vids / 19000.0))
    RuntimeAttr runtime_default = object {
        mem_gb: default_mem_gb,
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 50.0),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
        }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        set -euo pipefail
        ~{if defined(exclude_list) && !defined(exclude_list_idx) then "tabix -p bed ~{exclude_list}" else ""}

        #Run clustering
        svtk vcfcluster <(echo "~{vcf}") - \
                -d ~{dist} \
                -f ~{frac} \
                ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
                -z ~{svsize} \
                -p ~{vid_prefix} \
                -t ~{svtype} \
                -o ~{sample_overlap} \
                --preserve-ids \
                --preserve-genotypes \
                --preserve-header \
                | bcftools sort -O z -o ~{prefix}.vcf.gz - 
        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_out = "~{prefix}.vcf.gz"
        File vcf_out_idx = "~{prefix}.vcf.gz.tbi"
    }
}


task Extract_SVID_by_GC{
    input {
        File? SVID_GC
        File vcf
        File vcf_index
        String Genomic_Context
        String svtype
        Int min_size
        Int max_size
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    Float vcf_size = size(vcf, "GiB")
    Int vm_disk_size = ceil(vcf_size * 2)

    RuntimeAttr runtime_default = object {
        mem_gb: 1,
        disk_gb: vm_disk_size,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_benchmark_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, ".vcf.gz")

    command <<<
        set -euxo pipefail

        python /src/extract_SVs_by_svtype_genomiccontect.py \
        ~{vcf} \
        ~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz \
        ~{SVID_GC} \
        ~{svtype} \
        ~{Genomic_Context} \
        ~{min_size} \
        ~{max_size}

        tabix -p vcf ~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz 

    >>>

    output {
        File out_vcf = "~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz"
        File out_idx = "~{prefix}.~{svtype}.~{Genomic_Context}.vcf.gz.tbi"
    } 
}


