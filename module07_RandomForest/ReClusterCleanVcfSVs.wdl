version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "AnnotateGenomicContext.wdl" as AnnotateGenomicContext
import "ShardSVsByGenomicContext.wdl" as ShardSVsByGenomicContext

# Workflow to annotate vcf file with genomic context
workflow ReClusterCleanVcfSVs {
    input {
        File vcf
        File vcf_index
        File Repeat_Masks
        File Simple_Repeats
        File Segmental_Duplicates

        String prefix
        String vid_prefix
        String svtype
        String Genomic_Context

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

    call ShardSVsByGenomicContext.ShardSVsByGenomicContext as ShardSVsByGenomicContext {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            SVID_GC = AnnotateSVsWithGenomicContext.annotated_SVs,
            Genomic_Context = Genomic_Context,
            svtype = svtype,

            sv_base_mini_docker  =sv_base_mini_docker,
            sv_benchmark_docker = sv_benchmark_docker,

            runtime_attr_override_extract_svid = runtime_attr_override_extract_svid
    }

    call SvtkVcfCluster {
        input:
            vcf = ShardSVsByGenomicContext.sharded_vcf,
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

    call IntegrateReClusterdVcfs{
        input:
            vcf_all = vcf,
            vcf_all_idx = vcf_index,
            vcf_recluster = SvtkVcfCluster.vcf_out,
            vcf_recluster_idx = SvtkVcfCluster.vcf_idx,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_override_integrate_vcfs
    }

    call SortReClusterdVcfs{
        input:
            vcf_1 = IntegrateReClusterdVcfs.reclustered_Part1,
            vcf_2 = IntegrateReClusterdVcfs.reclustered_Part2,
            vcf_1_idx = IntegrateReClusterdVcfs.reclustered_Part1_idx,
            vcf_2_idx = IntegrateReClusterdVcfs.reclustered_Part2_idx,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_override_sort_reclustered_vcf
    }

    output{
        File SVID_GC = AnnotateSVsWithGenomicContext.annotated_SVs
        File resolve_vcf = SvtkVcfCluster.vcf_out
        File svid_annotation = IntegrateReClusterdVcfs.SVID_anno
        File out_vcf = SortReClusterdVcfs.sorted_vcf
        File out_vcf_idx = SortReClusterdVcfs.sorted_vcf_idx
    }
}



task SvtkVcfCluster {
    input {
        File vcf
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

    Float default_mem_gb = 3.75 + (120.0 * (num_vids / 19000.0) * (num_samples / 140000.0))
    RuntimeAttr runtime_default = object {
        mem_gb: default_mem_gb,
        disk_gb: ceil(10.0 + size(vcf, "GiB") * 5.0),
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
        File vcf_idx = "~{prefix}.vcf.gz.tbi"
    }
}

task IntegrateReClusterdVcfs{
    input{
        File vcf_all
        File vcf_all_idx
        File vcf_recluster
        File vcf_recluster_idx

        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 7.5,
        disk_gb: ceil(10.0 + size(vcf_all, "GiB")*4),
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

    String prefix = basename(vcf_recluster, ".vcf.gz")
    String prefix_all = basename(vcf_all, ".vcf.gz")
    
    command <<<
        svtk vcf2bed -i MEMBERS ~{vcf_recluster} ~{prefix}.bed
        cut -f4,7 ~{prefix}.bed > ~{prefix}.SVID_anno
        
        python3 <<CODE

        import os
        import pysam

        clustered_vs_origin_SVID = {}
        SVID_list = []
        fin=open("~{prefix}.SVID_anno")
        for line in fin:
            pin=line.strip().split()
            if not pin[0]=='name' and len(pin[1].split(','))>1:
                SVID_list.append(pin[0])
                for i in pin[1].split(','):
                    if not i in clustered_vs_origin_SVID.keys():
                        clustered_vs_origin_SVID[i] = pin[0]
        fin.close()

        fin = pysam.VariantFile("~{vcf_all}")
        fin2 = pysam.VariantFile("~{vcf_recluster}")
        fo = pysam.VariantFile("~{prefix_all}.ReClustered_unsort.vcf.gz", 'w', header = fin2.header)
        fo2 = pysam.VariantFile("~{prefix}.ReClustered_unsort.vcf.gz", 'w', header = fin2.header)
        for record in fin:
            if not record.id in clustered_vs_origin_SVID.keys():
                fo.write(record)
        fin.close()

        for record in fin2:
            if record.id in SVID_list:
                print(record.id)
                fo2.write(record)
        fin2.close()
        fo.close()
        fo2.close()
        print("done")
        CODE

        tabix -p vcf "~{prefix_all}.ReClustered_unsort.vcf.gz"
        tabix -p vcf "~{prefix}.ReClustered_unsort.vcf.gz"
    >>>

    output{
        File SVID_anno = "~{prefix}.SVID_anno"
        File reclustered_Part1 = "~{prefix_all}.ReClustered_unsort.vcf.gz"
        File reclustered_Part2 = "~{prefix}.ReClustered_unsort.vcf.gz"
        File reclustered_Part1_idx = "~{prefix_all}.ReClustered_unsort.vcf.gz.tbi"
        File reclustered_Part2_idx = "~{prefix}.ReClustered_unsort.vcf.gz.tbi"
    }
}

task SortReClusterdVcfs{
    input{
        File vcf_1
        File vcf_2
        File vcf_1_idx
        File vcf_2_idx

        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 7.5,
        disk_gb: ceil(10.0 + size(vcf_1, "GiB")*2 + size(vcf_2, "GiB")*2),
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

    String prefix_all = basename(vcf_1, ".ReClustered_unsort.vcf.gz")
    
    command <<<
        echo "~{vcf_1}" >> vcf_list
        echo "~{vcf_2}" >> vcf_list
        bcftools concat --allow-overlaps --output-type z --file-list vcf_list --output ~{prefix_all}.ReClustered.vcf.gz
        tabix -p vcf ~{prefix_all}.ReClustered.vcf.gz

    >>>

    output{
        File sorted_vcf = "~{prefix_all}.ReClustered.vcf.gz"
        File sorted_vcf_idx = "~{prefix_all}.ReClustered.vcf.gz.tbi"
    }
}

