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
import "AnnotateGenomicContext.wdl" as annotate_genomic_context
import "TasksBenchmark.wdl" as mini_tasks
workflow SplitPerSiteVCF{
    input{
        Array[File] cleanVcfs
        File simp_rep
        File seg_dup
        File rep_mask
        File? ped_file

        String prefix
        String sv_pipeline_docker
        String sv_base_mini_docker
        String sv_benchmark_docker

        RuntimeAttr? runtime_attr_override_vcf2bed
        RuntimeAttr? runtime_attr_override_anno_gc
        RuntimeAttr? runtime_attr_override_inte_gc
        RuntimeAttr? runtime_attr_override_concat_beds
        RuntimeAttr? runtime_attr_override_inte_features
    }

    scatter(vcf in cleanVcfs){
        call vcf2bed{
            input:
                vcf = vcf,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_vcf2bed
        }

        call AnnotateGenomicContext{
            input:
                bed_gz = vcf2bed.bed_gz,
                simp_rep = simp_rep,
                seg_dup = seg_dup,
                rep_mask = rep_mask,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override  = runtime_attr_override_anno_gc
        }

        call IntegrateGenomicContext{
            input:
                bed_gz = vcf2bed.bed_gz,
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

        call IntegrateSiteFeatures{
            input:
                bed_gz = vcf2bed.bed_gz,
                gc_anno = IntegrateGenomicContext.anno,
                sample_list = vcf2bed.sample_list,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_inte_features
        }
    }

    call mini_tasks.ConcatBeds as ConcatAnnos{
        input:
            prefix = "~{prefix}.anno",
            shard_bed_files = IntegrateSiteFeatures.anno,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_concat_beds
    }

    call mini_tasks.ConcatBeds as ConcatBeds{
        input:
            prefix = "~{prefix}.vcf2bed",
            shard_bed_files = vcf2bed.bed_gz,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_concat_beds
    }

    output{
        Array[File] vcf_to_bed = vcf2bed.bed_gz
        File anno_bed = ConcatAnnos.merged_bed_file
    }}


task vcf2bed{
    input{
        File vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(2.0 +  size(vcf, "GB")),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,'.vcf.gz')

    command<<<
        svtk vcf2bed -i SVTYPE -i SVLEN -i ALGORITHMS -i EVIDENCE  ~{vcf} ~{prefix}.bed
        bgzip ~{prefix}.bed
        zcat ~{vcf} | grep -m1 CHROM | cut -f10- | sed -e 's/\t/\n/g' > ~{prefix}.sample
    >>>

    output{
        File bed_gz = "~{prefix}.bed.gz"
        File sample_list = "~{prefix}.sample"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
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
        mem_gb: 5,
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

task IntegrateSiteFeatures{
    input {
        File bed_gz
        File gc_anno
        File sample_list
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 5,
        disk_gb: 20,
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

        Rscript /src/integrate_GenomicAnno_Freq_annotations.R \
        --bed ~{bed_gz} \
        --gc_anno ~{gc_anno} \
        --sample_list ~{sample_list} \
        --out ~{prefix}.anno 

        bgzip ~{prefix}.anno 

    >>>

    output {
        File anno = "~{prefix}.anno.gz"
    }
}




