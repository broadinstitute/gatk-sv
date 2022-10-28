##########################################################################################

## Copyright Broad Institute, 2022
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
import "CollectPEMetricsForCPX.wdl" as collect_pe_metrics_for_cpx

workflow ManualReview{
    input{ 
        File SVID_to_Remove #large CNVs, mCNVs that failed manual review, or redundant large CNVs that overlaps with Seg Dup and mCNVs
        
        File LINE1_Ref
        File HERVK_Ref

        Array[File] vcf_files #list of vcf to be manually revised
        Array[File] vcf_idxes

        Array[String] batch_name_list
        Array[File] PE_metrics
        Array[File] PE_metrics_idxes

        String prefix
        Int n_per_split
        File sample_PE_metrics


        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_override_vcf2bed 
        RuntimeAttr? runtime_attr_override_vcf2bed_sm
        RuntimeAttr? runtime_attr_override_extract_cpx_ctx
        RuntimeAttr? runtime_attr_override_extract_bnd_del
        RuntimeAttr? runtime_attr_override_bnd_vs_mei
        RuntimeAttr? runtime_attr_override_collect_pe
        RuntimeAttr? runtime_attr_override_split_script
        RuntimeAttr? runtime_attr_override_calcu_pe_stat
        RuntimeAttr? runtime_attr_override_concat_evidence
        RuntimeAttr? runtime_attr_override_generate_cpx_review_script
    }

    scatter(i in range(length(vcf_files))){

        call Vcf2Bed as Vcf2Bed{
            input:
                vcf = vcf_files[i],
                vcf_idx = vcf_idxes[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_vcf2bed
        }

        call SplitCpxCtx{
            input:
                bed = Vcf2Bed.bed,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_extract_cpx_ctx
        }

        call GenerateCpxReviewScript{
            input:
                bed = SplitCpxCtx.cpx_ctx_bed,
                sample_PE_metrics = sample_PE_metrics,
                sv_benchmark_docker = sv_benchmark_docker,
                runtime_attr_override = runtime_attr_override_generate_cpx_review_script
        }

        call collect_pe_metrics_for_cpx.CollectPEMetricsForCPX{
            input:
                batch_name_list = batch_name_list,
                PE_metrics = PE_metrics,
                PE_metrics_idxes = PE_metrics_idxes,
                PE_collect_script = GenerateCpxReviewScript.pe_evidence_collection_script,
                prefix = prefix,
                n_per_split = n_per_split,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override_collect_pe = runtime_attr_override_collect_pe,
                runtime_attr_override_split_script = runtime_attr_override_split_script,
                runtime_attr_override_calcu_pe_stat = runtime_attr_override_calcu_pe_stat,
                runtime_attr_override_concat_evidence = runtime_attr_override_concat_evidence
        }

        call SplitBndDel{
            input:
                vcf = vcf_files[i],
                vcf_idx = vcf_idxes[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_extract_bnd_del
        }

        call Vcf2Bed as vcf_to_bed_bnd_del{
            input:
                vcf = SplitBndDel.bnd_del_vcf,
                vcf_idx = SplitBndDel.bnd_del_vcf_idx,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_vcf2bed_sm
        }
        call BNDvsMEI{
            input:
                bed = vcf_to_bed_bnd_del.bed,
                LINE1 = LINE1_Ref,
                HERVK = HERVK_Ref,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_bnd_vs_mei
        }
    }
}
   
task Vcf2Bed{
    input{
        File vcf
        File vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(vcf, ".vcf.gz")
    command<<<
        set -euo pipefail
        svtk vcf2bed -i ALL --include-filters ~{vcf} ~{prefix}.bed
        bgzip ~{prefix}.bed
    >>>

    output{
        File bed = "~{prefix}.bed.gz"
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

task SplitCpxCtx{
    input{
        File bed
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(bed, ".bed.gz")

    command<<<
        set -eu
        zcat ~{bed} | head -1 > ~{prefix}.cpx_ctx.bed

        set -euo pipefail
        zcat ~{bed} | grep CPX | grep -v UNRESOLVE >> ~{prefix}.cpx_ctx.bed
        zcat ~{bed} | grep CTX >> ~{prefix}.cpx_ctx.bed
        zcat ~{bed} | grep INS | grep INV >> ~{prefix}.cpx_ctx.bed
        bgzip ~{prefix}.cpx_ctx.bed
    >>>

    output{
        File cpx_ctx_bed = "~{prefix}.cpx_ctx.bed.gz"
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

task GenerateCpxReviewScript{
    input{
        File bed
        File sample_PE_metrics
        String sv_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(bed, ".bed.gz")

    command<<<
        set -euo pipefail
    
        python /src/reformat_CPX_bed_and_generate_script.py \
        -i ~{bed}\
        -s ~{sample_PE_metrics} \
        -p CPX_CTX_disINS.PASS.PE_evidences \
        -c collect_PE_evidences.CPX_CTX_disINS.PASS.sh \
        -r ~{prefix}.svelter

    >>>

    output{
        File pe_evidence_collection_script = "collect_PE_evidences.CPX_CTX_disINS.PASS.sh"
        File svelter = "~{prefix}.svelter"
    }

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

task SplitBndDel{
    input{
        File vcf
        File vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(vcf, ".vcf.gz")
    command<<<
        set -euo pipefail

        python3 <<CODE

        import os
        import pysam
        fin=pysam.VariantFile("~{vcf}")
        fo=pysam.VariantFile("~{prefix}.bnd_del.vcf.gz", 'w', header = fin.header)
        for record in fin:
            if record.info['SVTYPE'] in ['BND']:
                if record.info['STRANDS']=="+=" and record.info['SVLEN']!=-1:
                    fo.write(record)
        fin.close()
        fo.close()
        CODE

        tabix -p vcf ~{prefix}.bnd_del.vcf.gz

    >>>

    output{
        File bnd_del_vcf = "~{prefix}.bnd_del.vcf.gz"
        File bnd_del_vcf_idx = "~{prefix}.bnd_del.vcf.gz.tbi"
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

task BNDvsMEI{
    input{
        File bed
        File LINE1
        File HERVK
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    command<<<
        set -euo pipefail

        bedtools coverage -wo -a ~{bed} -b ~{LINE1} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_LINE1/' > manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv 
        bedtools coverage -wo -a ~{bed} -b ~{HERVK} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_HERVK/' >> manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv 
    >>>

    output{
        File mei_del_SVID = "manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv"
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





