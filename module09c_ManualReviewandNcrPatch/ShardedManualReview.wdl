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
import "HailMerge.wdl" as HailMerge
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "RemoveDuplicateEvents.wdl" as remove_duplicate_events
import "CollectPEMetricsForCPX.wdl" as collect_pe_metrics_for_cpx
import "CollectLgCnvSupportForCPX.wdl" as collect_lg_cnv_supp_for_cpx
import "ReviseVcfWithManualResults.wdl" as revise_vcf_with_manual_results

workflow ShardedManualReview{
    input{ 
        File SVID_to_Remove #large CNVs, mCNVs that failed manual review, or redundant large CNVs that overlaps with Seg Dup and mCNVs
        File SVID_to_Remove_DEL_bump #SVID list of DELs <5Kb that are enriched for outlier samples; the outlier samples are defined as samples that carry more than average (top 5%) DELs between 400bp and 1Kb, AC>1 and AF<1%
        File CTX_manual
        #File CPX_manual
        #File duplicated_SVID_manual
        
        File LINE1_Ref
        File HERVK_Ref

        File vcf
        File vcf_idx
        String contig #use chromosome name if the input vcf is per-contig; else, put "whole_genome"

        Array[String] batch_name_list
        Array[File] PE_metrics
        Array[File] PE_metrics_idxes
        Array[File] Depth_DEL_beds
        Array[File] Depth_DUP_beds

        String prefix
        Int n_per_split
        File? raw_SVs
        File sample_PE_metrics
        File sample_depth_calls

        Boolean use_hail = false
        Boolean run_fix_ends = false
        Boolean clean_del_bump = false
        Boolean run_remove_duplicates = false
        String? gcs_project

        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_vcf2bed 
        RuntimeAttr? runtime_attr_vcf2bed_sm
        RuntimeAttr? runtime_attr_bnd_vs_mei
        RuntimeAttr? runtime_attr_collect_pe
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_split_script
        RuntimeAttr? runtime_attr_calcu_pe_stat
        RuntimeAttr? runtime_attr_extract_cpx_ctx
        RuntimeAttr? runtime_attr_extract_bnd_del
        RuntimeAttr? runtime_attr_concat_evidence
        RuntimeAttr? runtime_attr_split_raw_SVs_per_chr
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_preconcat_sharded_cluster
        RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
        RuntimeAttr? runtime_attr_fix_header_sharded_cluster
        RuntimeAttr? runtime_attr_generate_cpx_review_script
        RuntimeAttr? runtime_attr_remove_duplicate_events_task
        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_extract_cpx_lg_cnv_by_batch
        RuntimeAttr? runtime_attr_seek_depth_supp_for_cpx
        RuntimeAttr? runtime_attr_concat_bed_Step1
        RuntimeAttr? runtime_attr_concat_bed_Step2
        RuntimeAttr? runtime_attr_add_raw_SVs
        RuntimeAttr? runtime_attr_fix_bad_ends
        RuntimeAttr? runtime_attr_override_clean_del_bump
        RuntimeAttr? runtime_attr_calcu_cpx_evidences
    }


    call MiniTasks.ScatterVcf{
        input:
          vcf = vcf,
          prefix = prefix,
          records_per_shard = n_per_split,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_scatter_vcf
      }

    scatter (i in range(length(ScatterVcf.shards))) {

        if (run_remove_duplicates){
            call remove_duplicate_events.RemoveDuplicateEventsTask as RemoveDuplicateEvents{
                input:
                    vcf = ScatterVcf.shards[i],
                    vcf_index = ScatterVcf.shards_idx[i],
                    prefix = "~{prefix}.~{i}",
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_remove_duplicate_events_task
            }
        }

        call Vcf2Bed as Vcf2Bed{
            input:
                vcf = select_first([RemoveDuplicateEvents.deduplicated_vcf,ScatterVcf.shards[i]]),
                vcf_idx = select_first([RemoveDuplicateEvents.deduplicated_vcf_index,ScatterVcf.shards_idx[i]]),
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_vcf2bed
        }

        call SplitCpxCtx{
            input:
                bed = Vcf2Bed.bed,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_extract_cpx_ctx
        }

        call collect_lg_cnv_supp_for_cpx.CollectLgCnvSupportForCPX{
            input:
                cpx_ctx_bed = SplitCpxCtx.cpx_ctx_bed,
                prefix = "~{prefix}.~{i}",

                sample_depth_calls = sample_depth_calls,
                batch_name_list = batch_name_list,
                Depth_DEL_beds = Depth_DEL_beds,
                Depth_DUP_beds = Depth_DUP_beds,

                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
    
                runtime_attr_concat_bed_Step1 = runtime_attr_concat_bed_Step1,
                runtime_attr_concat_bed_Step2 = runtime_attr_concat_bed_Step2,
                runtime_attr_seek_depth_supp_for_cpx = runtime_attr_seek_depth_supp_for_cpx,
                runtime_attr_extract_cpx_lg_cnv_by_batch = runtime_attr_extract_cpx_lg_cnv_by_batch,
                runtime_attr_generate_cnv_segments_from_cpx = runtime_attr_generate_cnv_segments_from_cpx
        }

        call GenerateCpxReviewScript{
            input:
                bed = SplitCpxCtx.cpx_ctx_bed,
                sample_PE_metrics = sample_PE_metrics,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_generate_cpx_review_script
        }

        call collect_pe_metrics_for_cpx.CollectPEMetricsForCPX{
            input:
                batch_name_list = batch_name_list,
                PE_metrics = PE_metrics,
                PE_metrics_idxes = PE_metrics_idxes,
                PE_collect_script = GenerateCpxReviewScript.pe_evidence_collection_script,
                prefix = "~{prefix}.~{i}",
                n_per_split = n_per_split,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_collect_pe = runtime_attr_collect_pe,
                runtime_attr_split_script = runtime_attr_split_script,
                runtime_attr_calcu_pe_stat = runtime_attr_calcu_pe_stat,
                runtime_attr_concat_evidence = runtime_attr_concat_evidence
        }

        call CalculateCpxEvidences{
            input:
                PE_collect_script = GenerateCpxReviewScript.pe_evidence_collection_script,
                PE_supp = CollectPEMetricsForCPX.evi_stat,
                depth_supp = CollectLgCnvSupportForCPX.lg_cnv_depth_supp,
                prefix = "~{prefix}.~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_calcu_cpx_evidences
        }

        call CalculateCtxEvidences{
            input:
                PE_collect_script = GenerateCpxReviewScript.pe_evidence_collection_script,
                PE_supp = CollectPEMetricsForCPX.evi_stat,
                depth_supp = CollectLgCnvSupportForCPX.lg_cnv_depth_supp,
                prefix = "~{prefix}.~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_calcu_cpx_evidences
        }

        #### compare BNDs vs. MEIs to rescue mobile element deletions
        call SplitBndDel{
            input:
                vcf = select_first([RemoveDuplicateEvents.deduplicated_vcf,ScatterVcf.shards[i]]),
                vcf_idx = select_first([RemoveDuplicateEvents.deduplicated_vcf_index,ScatterVcf.shards_idx[i]]),
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_extract_bnd_del
        }

        call Vcf2Bed as vcf_to_bed_bnd_del{
            input:
                vcf = SplitBndDel.bnd_del_vcf,
                vcf_idx = SplitBndDel.bnd_del_vcf_idx,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_vcf2bed_sm
        }

        call BNDvsMEI{
            input:
                bed = vcf_to_bed_bnd_del.bed,
                LINE1 = LINE1_Ref,
                HERVK = HERVK_Ref,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_bnd_vs_mei
        }

        call revise_vcf_with_manual_results.ReviseVcfWithManualResults as ReviseVcfWithManualResults_wo_raw{
            input:
                vcf_file = select_first([RemoveDuplicateEvents.deduplicated_vcf,ScatterVcf.shards[i]]),
                vcf_index = select_first([RemoveDuplicateEvents.deduplicated_vcf_index,ScatterVcf.shards_idx[i]]),
                SVID_to_Remove = SVID_to_Remove,
                MEI_DEL_Rescue = BNDvsMEI.mei_del_SVID,
                CPX_manual = CalculateCpxEvidences.manual_revise_CPX_results,
                CTX_manual = CTX_manual,
                prefix = "~{prefix}.~{i}",
                contig = contig,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_pipeline_hail_docker = sv_pipeline_hail_docker
        }        

        File manual_revised_vcf = ReviseVcfWithManualResults_wo_raw.revised_vcf
        File manual_revised_vcf_idx = ReviseVcfWithManualResults_wo_raw.revised_vcf_idx


        if (clean_del_bump){
            call CleanDelBump{
                input:
                    vcf = manual_revised_vcf,
                    vcf_idx = manual_revised_vcf_idx,
                    SVID_to_Remove = SVID_to_Remove_DEL_bump,
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_override_clean_del_bump
            }
        }

        File sharded_annotated_vcf = select_first([CleanDelBump.cleaned_vcf, ReviseVcfWithManualResults_wo_raw.revised_vcf])
        File sharded_annotated_vcf_idx = select_first([CleanDelBump.cleaned_vcf_idx, ReviseVcfWithManualResults_wo_raw.revised_vcf_idx])
    }

    if (length(sharded_annotated_vcf) == 0) {
        call MiniTasks.GetVcfHeaderWithMembersInfoLine as GetVcfHeader_annotated {
            input:
                vcf_gz=vcf,
                prefix="~{prefix}.manual_review",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_vcf_header_with_members_info_line
        }
    }

    if (length(sharded_annotated_vcf) > 0) {
        if (use_hail) {
            call HailMerge.HailMerge as ConcatVcfsHail_annotated {
                input:
                    vcfs=sharded_annotated_vcf,
                    prefix="~{prefix}.manual_review",
                    gcs_project=gcs_project,
                    sv_base_mini_docker=sv_base_mini_docker,
                    sv_pipeline_docker=sv_pipeline_docker,
                    sv_pipeline_hail_docker=sv_pipeline_hail_docker,
                    runtime_attr_preconcat=runtime_attr_preconcat_sharded_cluster,
                    runtime_attr_hail_merge=runtime_attr_hail_merge_sharded_cluster,
                    runtime_attr_fix_header=runtime_attr_fix_header_sharded_cluster
            }
        }

        if (!use_hail) {
            call MiniTasks.ConcatVcfs as ConcatVcfs_annotated {
                input:
                    vcfs=sharded_annotated_vcf,
                    vcfs_idx=sharded_annotated_vcf_idx,
                    allow_overlaps=true,
                    outfile_prefix="~{prefix}.manual_review",
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_concat_sharded_cluster
            }
        }
    }

    call MiniTasks.ConcatFiles as concat_cpx_evidences{
        input:
            shard_bed_files = CalculateCpxEvidences.manual_revise_CPX_results,
            prefix = "~{prefix}.CPX_evidence" ,
            sv_base_mini_docker = sv_base_mini_docker

    }

    call MiniTasks.ConcatFiles as concat_ctx_evidences{
        input:
            shard_bed_files = CalculateCtxEvidences.manual_revise_CTX_results,
            prefix = "~{prefix}.CTX_evidence" ,
            sv_base_mini_docker = sv_base_mini_docker

    }


    File reviewed_vcf = select_first([GetVcfHeader_annotated.out, ConcatVcfs_annotated.concat_vcf, ConcatVcfsHail_annotated.merged_vcf])
    File reviewed_vcf_idx = select_first([GetVcfHeader_annotated.out_idx, ConcatVcfs_annotated.concat_vcf_idx, ConcatVcfsHail_annotated.merged_vcf_index])



    if (defined(raw_SVs)){

        if (contig!="whole_genome"){
            call SplitRawSVsPerChr{
                input:
                    raw_SVs = raw_SVs,
                    contig = contig,
                    sv_base_mini_docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_split_raw_SVs_per_chr
            }
        }

        File split_raw_SVs = select_first([SplitRawSVsPerChr.raw_SV_per_chr, raw_SVs])

        call revise_vcf_with_manual_results.AddRawSVs{
            input:
                prefix = prefix,
                batch_name = contig,
                vcf_file = reviewed_vcf,
                raw_SVs = split_raw_SVs,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_add_raw_SVs
        }
    }

    output{
        File revised_output_vcf = select_first([AddRawSVs.vcf_with_raw_SVs, reviewed_vcf])
        File revised_output_vcf_idx = select_first([AddRawSVs.vcf_idx_with_raw_SVs, reviewed_vcf_idx])
        File cpx_evidences = concat_cpx_evidences.merged_file
        File ctx_evidences = concat_ctx_evidences.merged_file
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
        String sv_base_mini_docker
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
        zcat ~{bed} | head -1 > ~{prefix}.cpx_ctx.bed

        filterColumn=$(zcat ~{bed} | head -1 | tr "\t" "\n" | awk '$1=="FILTER" {print NR}')

        zcat ~{bed} | awk 'NR > 1' | { grep CPX || true; } | awk -v filter_column=${filterColumn} '$filter_column ~ /UNRESOLVED/' >> ~{prefix}.cpx_ctx.bed

        zcat ~{bed} | awk 'NR > 1' | { grep CTX || true; } >> ~{prefix}.cpx_ctx.bed

        zcat ~{bed} | awk 'NR > 1' | { grep INS || true; } | { grep INV || true; } >> ~{prefix}.cpx_ctx.bed

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
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GenerateCpxReviewScript{
    input{
        File bed
        File sample_PE_metrics
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


    String prefix  = basename(bed, ".bed.gz")

    command<<<
        set -euo pipefail
    
        python /opt/sv-pipeline/scripts/manual_review/reformat_CPX_bed_and_generate_script.py \
        -i ~{bed} \
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
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GenerateCnvSegmentFromCpx{
    input{
        File bed
        File sample_depth_calls
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


    String prefix  = basename(bed, ".bed.gz")

    command<<<
        set -euo pipefail

        python <<CODE

        import os
        import sys
        def split_cpx_interval(info):
            #eg of info: DUP_chrY:3125606-3125667
            out = [info.split('_')[0], info.split('_')[1].split(':')[0]] +  [int(i) for i in info.split('_')[1].split(':')[1].split('-')]
            return out[1:4]+[out[0]]

        def sample_depth_calls_readin(sample_depth):
            out = {}
            fin=open(sample_depth)
            for line in fin:
                pin=line.strip().split()
                out[pin[0]] = pin[1:]
            fin.close()
            return out

        def readin_cpx_cnv(bed_input,sample_depth_hash):
            CPX_CNV = []
            f_bed = os.popen(r'''zcat %s'''%(bed_input))
            for line in f_bed:
                pin=line.strip().split('\t')
                if pin[0][0]=='#':
                    pos_CPX_TYPE = pin.index('CPX_TYPE')
                    pos_CPX_INTERVALS = pin.index('CPX_INTERVALS')
                    pos_SOURCE = pin.index('SOURCE')
                else:
                    interval = []
                    if 'DEL_' in pin[pos_CPX_INTERVALS] or 'DUP_' in pin[pos_CPX_INTERVALS]:
                        interval += [split_cpx_interval(i)+[pin[3]] for i in pin[pos_CPX_INTERVALS].split(",") if "DEL_" in i or "DUP_" in i]
                    if 'DEL_' in pin[pos_SOURCE] or 'DUP_' in pin[pos_SOURCE]:
                        interval += [split_cpx_interval(i)+[pin[3]] for i in pin[pos_SOURCE].split(",") if "DEL_" in i or "DUP_" in i]
                    if len(interval)>0:
                        for i in interval:
                            if i[2]-i[1]>4999:
                                sample_names = pin[5].split(',')
                                if i[3]=="DEL":
                                    for j in sample_names:
                                        CPX_CNV.append(i+[j, sample_depth_hash[j][0]])
                                if i[3]=="DUP":
                                    for j in sample_names:
                                        CPX_CNV.append(i+[j, sample_depth_hash[j][1]])
            f_bed.close()
            return CPX_CNV

        def write_cpx_cnv(CPX_CNV, fileout):
            fo=open(fileout, 'w') 
            for info in CPX_CNV:
                print('\t'.join(info), file=fo)
            fo.close()

        bed_input = "~{bed}"
        fileout = "~{prefix}.lg_CNV.bed"

        sample_depth_hash = sample_depth_calls_readin(sample_depth)
        CPX_CNV = readin_cpx_cnv(bed_input, sample_depth_hash)
        write_cpx_cnv(CPX_CNV, fileout)
        CODE

        bgzip "~{prefix}.lg_CNV.bed"
    >>>

    output{
        File cpx_lg_cnv = "~{prefix}.lg_CNV.bed.gz"
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

        python <<CODE

        import os
        import pysam
        fin=pysam.VariantFile("~{vcf}")
        fo=pysam.VariantFile("~{prefix}.bnd_del.vcf.gz", 'w', header = fin.header)
        for record in fin:
            if record.info['SVTYPE'] in ['BND'] and record.info['STRANDS']=="+-" and record.info['SVLEN']!=-1:
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

task CleanDelBump{
    input{
        File vcf
        File vcf_idx
        File SVID_to_Remove
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix = basename(vcf,"vcf.gz")
    command<<<
        set -euo pipefail
        python <<CODE

        def SVID_to_remove_readin(SVID_to_remove):
            out={}
            fin=open(SVID_to_remove)
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                    out[pin[0]] = pin[1]
            fin.close()
            return out

        def revise_vcf(vcf_input, vcf_output, hash_SVID_to_remove):
            fin=pysam.VariantFile(vcf_input)
            fo=pysam.VariantFile(vcf_output, 'w', header = fin.header)
            for record in fin:
                if record.id in hash_SVID_to_remove.keys(): 
                    if hash_SVID_to_remove[record.id]=='mCNV_under_5Kb':
                        continue
                    elif hash_SVID_to_remove[record.id]=='lg_CNV_failed_manua_review':
                        record.filter.add('FAIL_MANUAL_REVIEW')
                    elif hash_SVID_to_remove[record.id]=="OUTLIER_SAMPLE_ENRICHED":
                        record.filter.add('OUTLIER_SAMPLE_ENRICHED')
                fo.write(record)
            fin.close()
            fo.close()

        import os
        import pysam
        hash_SVID_to_remove = SVID_to_remove_readin("~{SVID_to_Remove}")
        revise_vcf("~{vcf}", "~{prefix}.DEL_bump_cleaned.vcf.gz", hash_SVID_to_remove)

        CODE

        tabix -p vcf "~{prefix}.DEL_bump_cleaned.vcf.gz"

    >>>

    output{
        File cleaned_vcf = "~{prefix}.DEL_bump_cleaned.vcf.gz"
        File cleaned_vcf_idx = "~{prefix}.DEL_bump_cleaned.vcf.gz.tbi"
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


task SplitRawSVsPerChr{
    input{
        File? raw_SVs
        String contig
        String sv_base_mini_docker
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
        touch "raw_SVs_to_add.~{contig}.tsv"
        awk '{if ($1=="~{contig}") print}' ~{raw_SVs} >> "raw_SVs_to_add.~{contig}.tsv"
    >>>

    output{
        File raw_SV_per_chr = "raw_SVs_to_add.~{contig}.tsv"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task CalculateCpxEvidences{
    input{
        File PE_collect_script
        File PE_supp
        File depth_supp
        String prefix
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

        set -eu
        awk '{print $6, $(NF-3)}' ~{PE_collect_script} | grep -v "_CTX_"| uniq > sample_SVID.tsv

        set -euo pipefail


        python <<CODE

        def sample_svid_readin(filein):
            out = {}
            fin=open(filein)
            for line in fin:
                pin=line.strip().split()
                if not pin[1] in out.keys():
                    out[pin[1]] = []
                out[pin[1]].append(pin[0])
            fin.close()
            return out

        def PE_supp_readin(filein):
            fin=os.popen(r'''zcat %s'''%(filein))
            out={}
            for line in fin:
                pin=line.strip().split()
                if not pin[4] in out.keys():
                    out[pin[4]] = {}
                if not pin[3] in out[pin[4]].keys():
                    out[pin[4]][pin[3]] = [[],[]]
                if pin[1]=="+":
                    out[pin[4]][pin[3]][0] = [int(pin[0])]+pin[1:3]
                elif pin[1]=="-":
                    out[pin[4]][pin[3]][1] = [int(pin[0])]+pin[1:3]
            fin.close()
            return out

        def Depth_supp_readin(filein):
            fin=os.popen(r'''zcat %s'''%(filein))
            out={}
            for line in fin:
                pin=line.strip().split()
                if not pin[0][0]=='#':
                    if not pin[4] in out.keys():
                        out[pin[4]] = {}
                    if not pin[5] in out[pin[4]].keys():
                        out[pin[4]][pin[5]] = []
                    out[pin[4]][pin[5]].append(float(pin[7]))
            fin.close()
            return out

        def add_PE_supp(sample_svid, PE_supp):
            out = {}
            for svid in sample_svid.keys():
                out[svid] = {}
                if svid in PE_supp.keys():
                    for sample in sample_svid[svid]:
                        if sample in PE_supp[svid].keys():
                            pe_supp = 'no_PE'
                            if len(PE_supp[svid][sample][0])==0:
                                PE_supp[svid][sample][0].append(0)
                            if len(PE_supp[svid][sample][1])==0:
                                PE_supp[svid][sample][1].append(0)
                            if PE_supp[svid][sample][0][0]>0 or PE_supp[svid][sample][1][0]>0:
                                pe_supp = 'partial_PE'
                            if PE_supp[svid][sample][0][0]>0 and PE_supp[svid][sample][1][0]>0:
                                pe_supp = 'low_PE'
                            if PE_supp[svid][sample][0][0]>2 and PE_supp[svid][sample][1][0]>2:
                                pe_supp = 'high_PE'
                            out[svid][sample] = [pe_supp]
                        else:
                            out[svid][sample] = ['no_PE']
                else:
                    for sample in sample_svid[svid]:
                        out[svid][sample] = ['no_PE']
            return out

        def add_depth_supp(sample_svid_pe, Depth_supp):
            for svid in sample_svid_pe.keys():
                if svid in Depth_supp.keys():
                    for sample in sample_svid_pe[svid].keys():
                        if sample in Depth_supp[svid].keys():
                            depth_supp = 'lack_depth'
                            if Depth_supp[svid][sample][0]>.5:
                                depth_supp = 'depth'
                                if len(Depth_supp[svid][sample])>1:
                                    if not Depth_supp[svid][sample][1]>.5:
                                        depth_supp+=',lack_depth'
                            else:
                                if len(Depth_supp[svid][sample])>1:
                                    if Depth_supp[svid][sample][1]>.5:
                                        depth_supp+=',depth'
                            sample_svid_pe[svid][sample].append(depth_supp)
                        else:
                            sample_svid_pe[svid][sample].append('NA')
                else:
                    for sample in sample_svid_pe[svid].keys():
                        sample_svid_pe[svid][sample].append('NA')
            return sample_svid_pe

        def write_pe_depth_supp(sample_svid_pe_depth, fileout):
            fo=open(fileout,'w')
            print('\t'.join(['VID','sample','supportive_PE_counts','depth_supp']), file=fo)
            for svid in sample_svid_pe_depth.keys():
                for sample in sample_svid_pe_depth[svid].keys():
                    print('\t'.join([svid,sample] + sample_svid_pe_depth[svid][sample]), file=fo)
            fo.close()

        import os
        sample_svid = sample_svid_readin("sample_SVID.tsv")
        PE_supp = PE_supp_readin("~{PE_supp}")
        Depth_supp = Depth_supp_readin("~{depth_supp}")
        sample_svid_pe = add_PE_supp(sample_svid, PE_supp)
        sample_svid_pe_depth = add_depth_supp(sample_svid_pe, Depth_supp)
        write_pe_depth_supp(sample_svid_pe_depth, "~{prefix}.manual_revise.CPX_results")

        CODE
    >>>

    output{
        File manual_revise_CPX_results = "~{prefix}.manual_revise.CPX_results"
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
   
task CalculateCtxEvidences{
    input{
        File PE_collect_script
        File PE_supp
        File depth_supp
        String prefix
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

        set -eu

        awk '{print $6, $(NF-3)}' ~{PE_collect_script} | grep "_CTX_"| uniq > sample_SVID.tsv

        set -eu
        
        zcat ~{PE_supp} | grep "_CTX_"| uniq > PE_supp.tsv

        set -euo pipefail

        python <<CODE

        def sample_svid_readin(filein):
            out = {}
            fin=open(filein)
            for line in fin:
                pin=line.strip().split()
                if not pin[1] in out.keys():
                    out[pin[1]] = []
                out[pin[1]].append(pin[0])
            fin.close()
            return out

        def PE_supp_readin(filein):
            fin=open(filein)
            out={}
            for line in fin:
                pin=line.strip().split()
                if not pin[4] in out.keys():
                    out[pin[4]] = {}
                if not pin[3] in out[pin[4]].keys():
                    out[pin[4]][pin[3]] = [[],[]]
                if pin[1]=="+":
                    out[pin[4]][pin[3]][0] = [int(pin[0])]+pin[1:3]
                elif pin[1]=="-":
                    out[pin[4]][pin[3]][1] = [int(pin[0])]+pin[1:3]
            fin.close()
            return out

        def add_PE_supp(sample_svid, PE_supp):
            out = {}
            for svid in sample_svid.keys():
                out[svid] = {}
                if svid in PE_supp.keys():
                    for sample in sample_svid[svid]:
                        if sample in PE_supp[svid].keys():
                            pe_supp = 'no_PE'
                            if len(PE_supp[svid][sample][0])==0:
                                PE_supp[svid][sample][0].append(0)
                            if len(PE_supp[svid][sample][1])==0:
                                PE_supp[svid][sample][1].append(0)
                            if PE_supp[svid][sample][0][0]>0 or PE_supp[svid][sample][1][0]>0:
                                pe_supp = 'partial_PE'
                            if PE_supp[svid][sample][0][0]>0 and PE_supp[svid][sample][1][0]>0:
                                pe_supp = 'low_PE'
                            if PE_supp[svid][sample][0][0]>2 and PE_supp[svid][sample][1][0]>2:
                                pe_supp = 'high_PE'
                            out[svid][sample] = [pe_supp]
                        else:
                            out[svid][sample] = ['no_PE']
                else:
                    for sample in sample_svid[svid]:
                        out[svid][sample] = ['no_PE']
            return out

        def write_pe_depth_supp(sample_svid_pe_depth, fileout):
            fo=open(fileout,'w')
            print('\t'.join(['VID','sample','supportive_PE_counts']), file=fo)
            for svid in sample_svid_pe_depth.keys():
                for sample in sample_svid_pe_depth[svid].keys():
                    print('\t'.join([svid,sample] + sample_svid_pe_depth[svid][sample]), file=fo)
            fo.close()

        import os
        sample_svid = sample_svid_readin("sample_SVID.tsv")
        PE_supp = PE_supp_readin("PE_supp.tsv")
        sample_svid_pe = add_PE_supp(sample_svid, PE_supp)
        write_pe_depth_supp(sample_svid_pe, "~{prefix}.manual_revise.CTX_results")

        CODE
    >>>

    output{
        File manual_revise_CTX_results = "~{prefix}.manual_revise.CTX_results"
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
   


