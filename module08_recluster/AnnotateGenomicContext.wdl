version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"

# Workflow to annotate vcf file with genomic context
workflow AnnotateSVsWithGenomicContext {
    input {
        File vcf
        File vcf_index
        File Repeat_Masks
        File Simple_Repeats
        File Segmental_Duplicates

        String sv_base_mini_docker
        String sv_benchmark_docker
        String sv_pipeline_docker

        # overrides for MiniTasks
        RuntimeAttr? runtime_override_extract_SV_sites
        RuntimeAttr? runtime_override_vcf_to_bed
        RuntimeAttr? runtime_attr_override_anno_gc
        RuntimeAttr? runtime_attr_override_inte_gc
    }

    call ExtractSitesFromVcf {
        input:
            vcf = vcf,
            vcf_index = vcf_index,
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_override_extract_SV_sites
    }

    call Vcf2Bed{
        input:
            vcf = ExtractSitesFromVcf.out,
            vcf_index = ExtractSitesFromVcf.out_idx, 
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override=runtime_override_vcf_to_bed
    }

    call AnnotateGenomicContext{
        input:
            bed_gz = Vcf2Bed.out,
            simp_rep = Simple_Repeats,
            seg_dup = Segmental_Duplicates,
            rep_mask = Repeat_Masks,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_anno_gc
    }

    call IntegrateGenomicContext{
        input:
            bed_gz = Vcf2Bed.out,
            le_bp_vs_sr = AnnotateGenomicContext.le_bp_vs_sr,
            le_bp_vs_sd = AnnotateGenomicContext.le_bp_vs_sd,
            le_bp_vs_rm = AnnotateGenomicContext.le_bp_vs_rm,
            ri_bp_vs_sr = AnnotateGenomicContext.ri_bp_vs_sr,
            ri_bp_vs_sd = AnnotateGenomicContext.ri_bp_vs_sd,
            ri_bp_vs_rm = AnnotateGenomicContext.ri_bp_vs_rm,
            lg_cnv_vs_sr = AnnotateGenomicContext.lg_cnv_vs_sr,
            lg_cnv_vs_sd = AnnotateGenomicContext.lg_cnv_vs_sd,
            lg_cnv_vs_rm = AnnotateGenomicContext.lg_cnv_vs_rm,
            sv_benchmark_docker = sv_benchmark_docker,
            runtime_attr_override = runtime_attr_override_inte_gc
    }

    output{
        File annotated_SVs = IntegrateGenomicContext.anno
    }
}

task ExtractSitesFromVcf {
    input {
        File vcf
        File vcf_index
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    Float vcf_size = size(vcf, "GiB")
    Int vm_disk_size = ceil(vcf_size * 1.2)

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, ".vcf.gz")
    command <<<
        set -euxo pipefail

        zcat ~{vcf} | cut -f1-10 |bgzip > ~{prefix}.sites.vcf.gz
        tabix ~{prefix}.sites.vcf.gz
    >>>

    output {
        File out = "~{prefix}.sites.vcf.gz"
        File out_idx = "~{prefix}.sites.vcf.gz.tbi"
    }
}

task Vcf2Bed{
    input {
        File vcf
        File vcf_index
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float vcf_size = size(vcf, "GiB")
    Int vm_disk_size = ceil(vcf_size * 1.2)

    RuntimeAttr runtime_default = object {
        mem_gb: 2,
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
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, ".vcf.gz")
    command <<<
        set -euxo pipefail

        svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} ~{prefix}.bed
        cut -f1-4,7,8 ~{prefix}.bed | bgzip > ~{prefix}.bed.gz
    >>>

    output {
        File out = "~{prefix}.bed.gz"
    }        
}

task AnnotateGenomicContext{
    input {
        File bed_gz
        File simp_rep
        File seg_dup
        File rep_mask
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr runtime_default = object {
        mem_gb: 7.5,
        disk_gb: 15,
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
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(bed_gz, ".bed.gz")
    command <<<
        set -euxo pipefail

        zcat ~{bed_gz} | awk '{print $1,$2,$2,$4,$5}'    | sed -e 's/ /\t/g' > ~{prefix}.le_bp
        zcat ~{bed_gz} | awk '{print $1,$3,$3,$4,$5}'    | sed -e 's/ /\t/g' > ~{prefix}.ri_bp
        zcat ~{bed_gz} | awk '{if ($5=="DEL" || $5=="DUP" || $5=="CNV" ) print}' | awk '{if ($3-$2>5000) print}' | cut -f1-5 >    ~{prefix}.lg_cnv

        bedtools coverage -a ~{prefix}.le_bp -b ~{simp_rep} | awk '{if ($9>0) print}'> ~{prefix}.le_bp.vs.SR
        bedtools coverage -a ~{prefix}.le_bp -b ~{seg_dup}    | awk '{if ($9>0) print}'> ~{prefix}.le_bp.vs.SD
        bedtools coverage -a ~{prefix}.le_bp -b ~{rep_mask}    | awk '{if ($9>0) print}'> ~{prefix}.le_bp.vs.RM

        bedtools coverage -a ~{prefix}.ri_bp -b ~{simp_rep} | awk '{if ($9>0) print}'> ~{prefix}.ri_bp.vs.SR
        bedtools coverage -a ~{prefix}.ri_bp -b ~{seg_dup}    | awk '{if ($9>0) print}'> ~{prefix}.ri_bp.vs.SD
        bedtools coverage -a ~{prefix}.ri_bp -b ~{rep_mask}    | awk '{if ($9>0) print}'> ~{prefix}.ri_bp.vs.RM

        bedtools coverage -a ~{prefix}.lg_cnv -b ~{simp_rep}    > ~{prefix}.lg_cnv.vs.SR
        bedtools coverage -a ~{prefix}.lg_cnv -b ~{seg_dup}    > ~{prefix}.lg_cnv.vs.SD
        bedtools coverage -a ~{prefix}.lg_cnv -b ~{rep_mask}    > ~{prefix}.lg_cnv.vs.RM

    >>>

    output {
        File le_bp_vs_sr = "~{prefix}.le_bp.vs.SR"
        File le_bp_vs_sd = "~{prefix}.le_bp.vs.SD"
        File le_bp_vs_rm = "~{prefix}.le_bp.vs.RM"
        File ri_bp_vs_sr = "~{prefix}.ri_bp.vs.SR"
        File ri_bp_vs_sd = "~{prefix}.ri_bp.vs.SD"
        File ri_bp_vs_rm = "~{prefix}.ri_bp.vs.RM"
        File lg_cnv_vs_sr = "~{prefix}.lg_cnv.vs.SR"
        File lg_cnv_vs_sd = "~{prefix}.lg_cnv.vs.SD"
        File lg_cnv_vs_rm = "~{prefix}.lg_cnv.vs.RM"
    }        
}

task IntegrateGenomicContext{
    input {
        File bed_gz
        File le_bp_vs_sr
        File le_bp_vs_sd
        File le_bp_vs_rm
        File ri_bp_vs_sr
        File ri_bp_vs_sd
        File ri_bp_vs_rm
        File lg_cnv_vs_sr
        File lg_cnv_vs_sd
        File lg_cnv_vs_rm

        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr runtime_default = object {
        mem_gb: 1,
        disk_gb: 10,
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

    String prefix = basename(bed_gz, ".bed.gz")
    command <<<
        set -euxo pipefail

        Rscript /src/integrate_Genomic_Content_annotations.R \
        --bed ~{bed_gz} \
        --out ~{prefix}.GC \
        --le_bp_vs_sr ~{le_bp_vs_sr} \
        --le_bp_vs_sd ~{le_bp_vs_sd} \
        --le_bp_vs_rm ~{le_bp_vs_rm} \
        --ri_bp_vs_sr ~{ri_bp_vs_sr} \
        --ri_bp_vs_sd ~{ri_bp_vs_sd} \
        --ri_bp_vs_rm ~{ri_bp_vs_rm} \
        --lg_cnv_vs_sr ~{lg_cnv_vs_sr} \
        --lg_cnv_vs_sd ~{lg_cnv_vs_sd} \
        --lg_cnv_vs_rm ~{lg_cnv_vs_rm}

    >>>

    output {
        File anno = "~{prefix}.GC"
    }        
}










