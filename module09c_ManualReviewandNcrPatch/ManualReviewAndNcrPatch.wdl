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
import "Utils.wdl" as Utils
import "HailMerge.wdl" as HailMerge
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "RemoveDuplicateEvents.wdl" as remove_duplicate_events
import "EnforceMinNoCallRate.wdl" as enforce_min_no_call_rate
import "ApplyNoCallRateCutoffs.wdl" as apply_no_call_rate_cutoffs

workflow ManualReviewAndNcrPatch{
    input{ 
        File SVID_to_Remove_DEL_bump #SVID list of DELs <5Kb that are enriched for outlier samples; the outlier samples are defined as samples that carry more than average (top 10%) DELs between 400bp and 1Kb, AC>1 and AF<1%
        File SVID_to_Remove_lg_DUP_lack_raw_depth #2 column list of SVID and samples of large DUPs that lack raw depth support

        File vcf
        File vcf_idx
        File? ctx_vcf
        File? ctx_vcf_idx
        File before_clustered_vcf      #the vcf after outlier removal but before SV re-cluster 
        File before_clustered_vcf_idx  #the vcf idx after outlier removal but before SV re-cluster 
        String contig
        String prefix
        Int n_per_split

        Boolean exclude_CTX = true
        Boolean add_CTX = false
        Boolean use_hail = false
        String? gcs_project

        Array[String]? sample_subset_prefixes
        Array[File]? sample_subset_lists

        Array[String]? svtype_list
        Array[String]? ncr_filter_field
        Array[Float]? NCR_cff_list
        String? global_ncr_filter_field
        Float? global_max_ncr
        Float? NCR_clean_del #if provided, use this NCR cutoff for DELs between 400bp and 1Kb
        Float? NCR_clean_dup #if provided, use this NCR cutoff for DUPs between 100bp to 1Kb


        String sv_benchmark_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_additional_filters 
        RuntimeAttr? runtime_attr_additional_filters_Step2
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_preconcat_sharded_cluster
        RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
        RuntimeAttr? runtime_attr_fix_header_sharded_cluster
        RuntimeAttr? runtime_attr_override_split_svid_by_contig
        RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
        RuntimeAttr? runtime_attr_index_vcf
        RuntimeAttr? runtime_attr_split_vcf_by_type
        RuntimeAttr? runtime_attr_apply_ncr_filter
        RuntimeAttr? runtime_attr_concat_filtered_vcfs
        RuntimeAttr? runtime_attr_exclude_type_from_vcf
        RuntimeAttr? runtime_attr_combine_ctx
        }


    #clean up the DEL bump: reapply the SVID_to_Remove_DEL_bump
    call SplitSvidByContig{
        input:
            SVID = SVID_to_Remove_lg_DUP_lack_raw_depth,
            contig = contig,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_override_split_svid_by_contig
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

        call AdditionalFilters{
            input:
                vcf = ScatterVcf.shards[i],
                vcf_idx = ScatterVcf.shards_idx[i],
                SVID_to_Remove_DEL_bump = SVID_to_Remove_DEL_bump,
                SVID_to_Remove_lg_DUP_lack_raw_depth = SplitSvidByContig.SVID_split,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_additional_filters
        }

        call AdditionalFiltersStep2{
            input:
                vcf = AdditionalFilters.revised_vcf,
                vcf_idx = AdditionalFilters.revised_vcf_idx,
                before_clustered_vcf = before_clustered_vcf,
                before_clustered_vcf_idx = before_clustered_vcf_idx,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_additional_filters_Step2
        }

        call enforce_min_no_call_rate.AnnotateNCRs as AnnotateNCRs {
            input:
                vcf=AdditionalFiltersStep2.revised_vcf,
                subset_prefixes=sample_subset_prefixes,
                subset_lists=sample_subset_lists,
                sv_pipeline_base_docker=sv_pipeline_base_docker
        }

        call apply_no_call_rate_cutoffs.ApplyNoCallRateCutoffs{
            input:
                vcf = AnnotateNCRs.annotated_vcf,
                svtype_list  = select_first([svtype_list]),
                ncr_filter_field = select_first([ncr_filter_field]),
                NCR_cff_list = select_first([NCR_cff_list]),
                global_max_ncr = global_max_ncr,
                global_ncr_filter_field = global_ncr_filter_field,
                exclude_CTX = exclude_CTX,
                NCR_clean_del = NCR_clean_del,
                NCR_clean_dup = NCR_clean_dup,

                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker,

                runtime_attr_index_vcf = runtime_attr_index_vcf,
                runtime_attr_split_vcf_by_type = runtime_attr_split_vcf_by_type,
                runtime_attr_apply_ncr_filter = runtime_attr_apply_ncr_filter,
                runtime_attr_concat_filtered_vcfs = runtime_attr_concat_filtered_vcfs, 
                runtime_attr_exclude_type_from_vcf = runtime_attr_exclude_type_from_vcf
            }

        call AdditionalFiltersStep3{
            input:
                vcf = ApplyNoCallRateCutoffs.ncr_filtered_vcf,
                vcf_idx = ApplyNoCallRateCutoffs.ncr_filtered_vcf_idx,
                SVID_to_Remove_DEL_bump = SVID_to_Remove_DEL_bump,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_additional_filters
        }

        File sharded_revised_vcf = AdditionalFiltersStep3.revised_vcf
        File sharded_revised_vcf_idx = AdditionalFiltersStep3.revised_vcf_idx
    }

    if (length(sharded_revised_vcf) == 0) {
        call MiniTasks.GetVcfHeaderWithMembersInfoLine as GetVcfHeader_revised {
            input:
                vcf_gz=vcf,
                prefix="~{prefix}.manual_review_patch",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_vcf_header_with_members_info_line
        }
    }

    if (length(sharded_revised_vcf) > 0) {
        if (use_hail) {
            call HailMerge.HailMerge as ConcatVcfsHail_revised {
                input:
                    vcfs=sharded_revised_vcf,
                    prefix="~{prefix}.manual_review_patch",
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
            call MiniTasks.ConcatVcfs as ConcatVcfs_revised {
                input:
                    vcfs=sharded_revised_vcf,
                    vcfs_idx=sharded_revised_vcf_idx,
                    allow_overlaps=true,
                    outfile_prefix="~{prefix}.manual_review_patch",
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_concat_sharded_cluster
            }
        }
    }

    File reviewed_vcf = select_first([GetVcfHeader_revised.out, ConcatVcfs_revised.concat_vcf, ConcatVcfsHail_revised.merged_vcf])
    File reviewed_vcf_idx = select_first([GetVcfHeader_revised.out_idx, ConcatVcfs_revised.concat_vcf_idx, ConcatVcfsHail_revised.merged_vcf_index])

    if (add_CTX){
        call CombineWithCtx{
            input:
                vcf = reviewed_vcf,
                vcf_idx = reviewed_vcf_idx,
                contig = contig,
                ctx_vcf = ctx_vcf,
                ctx_vcf_idx = ctx_vcf_idx,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_combine_ctx

        }
    }


    output{
        File output_vcf = select_first([reviewed_vcf, CombineWithCtx.combined_vcf])
        File output_vcf_idx = select_first([reviewed_vcf_idx, CombineWithCtx.combined_vcf_idx])
    }
}



# this function is added as a patch to manual review; it'll do: 
#a. tag DELs enriched outliers (top 10%) that carry the most DELs between 400bp and 1Kb, AF<1% and AC>1;  
#b. tag large DUPs (>5Kb) that are covered by SD or SR (>30%); 
#c. MEI DEL: revise end to pos+svlen
#d. remove any SVs that overlaps IGH_MHC
task AdditionalFilters{
    input{
        File vcf
        File vcf_idx
        File SVID_to_Remove_DEL_bump
        File SVID_to_Remove_lg_DUP_lack_raw_depth
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
                    out[pin[0]] = [pin[1]]
                else:
                    out[pin[0]] += [pin[1]]
            fin.close()
            return out

        def SVID_sample_to_remove_readin(SVID_sample_to_remove):
            out={}
            fin=os.popen(r'''zcat %s'''%(SVID_sample_to_remove))
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                    out[pin[0]] = []
                out[pin[0]].append(pin[1])
            fin.close()
            return out

        def count_carrier(record):
            ref = 0
            alt = 0
            for sample in record.samples.keys():
                if record.samples[sample]['GT'] == (0,0):
                    ref +=1
                elif record.samples[sample]['GT'] == (0,1):
                    alt +=1
                elif record.samples[sample]['GT'] == (1,0):
                    alt +=1
                elif record.samples[sample]['GT'] == (1,1):
                    alt +=1
            out = 0
            if alt+ref>0:
                out = alt/(alt+ref)
            return out

        def revise_vcf(vcf_input, vcf_output, hash_SVID_to_remove, hash_SVID_sample_to_remove):

            fin=pysam.VariantFile(vcf_input)
            #### revise header
            header = fin.header            
            header.filters.add('IGH_MHC_OVERLAP', None, None, "SVs that are overlapped by over 50% by IGH or MCH regions, these variants are of low confidence")
            #header.filters.remove_header('OUTLIER_SAMPLE_ENRICHED')
            #header.filters.add('OUTLIER_SAMPLE_ENRICHED', None, None, "Deletion is enriched for non-reference genotypes in outlier samples (5%), likely indicating noisy or unreliable genotypes")
            header.info.add('LIKELY_REFERENCE_ARTIFACT', None, "Flag", "SV sites carried by >99% samples likely indicate reference or technical artifact")
            header.info.add('LOW_CONFIDENCE_REPETITIVE_LARGE_DUP', None, "Flag", "Duplications over 5Kb that are overlapped by segmental duplicates or simple repeats, these variants are not as high confident as others because of the affect by the repetitive sequences")
            header.info.add('OUTLIER_SAMPLE_ENRICHED_LENIENT', None, "Flag", "Deletion is enriched for non-reference genotypes in outlier samples (10%), likely indicating noisy or unreliable genotypes")
            #### readin vcf
            fo=pysam.VariantFile(vcf_output, 'w', header = header)
            for record in fin:
                carrier_prop = count_carrier(record)

                if carrier_prop > .99:                                                  #label sites that are carried by >99% samples 
                    record.info['LIKELY_REFERENCE_ARTIFACT'] = True

                if record.id in hash_SVID_to_remove.keys(): 

                    if 'mCNV_under_5Kb' in hash_SVID_to_remove[record.id]:    #remove mCNVs <5Kb from the callset
                        continue

                    if 'Overlap_SD_SR_over_30perc' in hash_SVID_to_remove[record.id]:   #label DUPs >5Kb if >50% of them are overlapped by either SD or SR sequences 
                        if record.info['SVTYPE']=='DUP' and record.info['SVLEN']>4999:
                            record.info['LOW_CONFIDENCE_REPETITIVE_LARGE_DUP'] = True

                    if 'lg_CNV_failed_manua_review' in hash_SVID_to_remove[record.id]:  #label large CNVs >1Mb that failed manual review
                        record.filter.add('FAIL_MANUAL_REVIEW')

                    if 'OUTLIER_SAMPLE_ENRICHED_5perc' in hash_SVID_to_remove[record.id]:     #label DUPs <5Kb if >50% carriers are outliers defined as the top 5% samples that carry most DUPs between 400bp - 1Kb, AF<1% and AC>1
                        record.filter.add('OUTLIER_SAMPLE_ENRICHED')

                    if 'OUTLIER_SAMPLE_ENRICHED_5to10perc' in hash_SVID_to_remove[record.id]:     #label DUPs <5Kb if >50% carriers are outliers defined as the top 10% samples that carry most DUPs between 400bp - 1Kb, AF<1% and AC>1
                        record.info['OUTLIER_SAMPLE_ENRICHED_LENIENT'] = True

                    if 'IGH_MHC_OVERLAP' in hash_SVID_to_remove[record.id]:       #label SVs if they are overlapped by IGH or MCH sequences
                        record.filter.add('IGH_MHC_OVERLAP')

                if record.id in hash_SVID_sample_to_remove.keys():                      #revise alt GT to ./. for large DUPs >5Kb that are overlapped by SD or SR by <50% if the sample does not have raw depth support
                    for sample in hash_SVID_sample_to_remove[record.id]:
                        if sample in record.samples.keys():
                            if record.samples[sample]['GT']==(0,1) or record.samples[sample]['GT']==(1,1):
                                record.samples[sample]['GT'] = [None,None]
                if record.alts == ('<DEL:ME:LINE1>',) or record.alts == ('<DEL:ME:HERVK>',):    #correct END breakpoint of ME DELs
                    record.stop = record.pos + record.info['SVLEN']
                fo.write(record)
            fin.close()
            fo.close()


        import os
        import pysam
        hash_SVID_to_remove = SVID_to_remove_readin("~{SVID_to_Remove_DEL_bump}")
        hash_SVID_sample_to_remove = SVID_sample_to_remove_readin("~{SVID_to_Remove_lg_DUP_lack_raw_depth}")
        revise_vcf("~{vcf}", "~{prefix}.revised.vcf.gz", hash_SVID_to_remove, hash_SVID_sample_to_remove)

        CODE

        tabix -p vcf "~{prefix}.revised.vcf.gz"

    >>>

    output{
        File revised_vcf = "~{prefix}.revised.vcf.gz"
        File revised_vcf_idx = "~{prefix}.revised.vcf.gz.tbi"
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

# this function is added as a patch to manual review; it'll do:
# remove wham SR only DEL <1Kb, and null out genotypes uniquely contributed by these DELs to the clustered SVs
task AdditionalFiltersStep2{
    input{
        File vcf
        File vcf_idx
        File before_clustered_vcf
        File before_clustered_vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 200,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix = basename(vcf,"vcf.gz")
    command<<<
        set -euo pipefail
        python <<CODE

        def extract_wham_SR_DEL_under1Kb_from_before_clustered_vcf(before_clustered_vcf, SVID_hash):
            SVID_list = []
            for i in SVID_hash.keys():
                SVID_list += SVID_hash[i]
            SVID_sample_hash = {}
            wham_SR_DEL_SVID = []
            fin=pysam.VariantFile(before_clustered_vcf)
            for record in fin:
                if record.id in SVID_list:
                    SVID_sample_hash[record.id] = [sample for sample in record.samples.keys() if record.samples[sample]['GT'] in [(0,1),(1,1)] ]
                    if record.info['SVTYPE']=="DEL" and record.info['ALGORITHMS']==('wham',) and record.info['EVIDENCE']==('SR',) and record.info['SVLEN']<1001:
                        wham_SR_DEL_SVID.append(record.id)
            fin.close()
            out_hash ={}
            for i in SVID_hash.keys():
                wham_sr_del_gt = []
                other_gt = []
                for j in SVID_hash[i]:
                    if j in SVID_sample_hash.keys():
                        if j in wham_SR_DEL_SVID:
                            wham_sr_del_gt += SVID_sample_hash[j]
                        else:
                            other_gt += SVID_sample_hash[j]
                wham_sr_del_unique_gt = [k for k in wham_sr_del_gt if not k in other_gt]
                if len(wham_sr_del_unique_gt)>0:
                    out_hash[i] = wham_sr_del_unique_gt
            return out_hash


        def collect_SVID_from_DEL_clusters(vcf):
            # this function collexcts the members of DELs that were clustered during the SV recluster steps
            out = {}
            fin=pysam.VariantFile(vcf)
            for record in fin:
                if "Redun_Removed." in record.id or "MEMBERS" in record.info.keys():
                    if record.info['SVTYPE']=="DEL":
                        out[record.id] = record.info["MEMBERS"]
            fin.close()
            return out

        def revise_wham_SR_DEL_under1Kb(vcf, vcf_out, hash_SVID_wham_SR_DEL_under1Kb_unique_GT):
            fin=pysam.VariantFile(vcf)
            header = fin.header
            header.filters.add('LOWQUAL_WHAM_SR_DEL', None, None, "deletions under1Kb that are uniquely from wham and have SR-only support")
            fo=pysam.VariantFile(vcf_out, 'w', header = header)
            for record in fin:
                if record.info['SVTYPE']=="DEL" and record.info['ALGORITHMS']==('wham',) and record.info['EVIDENCE']==('SR',) and record.info['SVLEN']<1001:
                    record.filter.add('LOWQUAL_WHAM_SR_DEL')
                if record.id in hash_SVID_wham_SR_DEL_under1Kb_unique_GT.keys():
                    for sample in hash_SVID_wham_SR_DEL_under1Kb_unique_GT[record.id]:
                        record.samples[sample]['GT'] = [None,None]
                fo.write(record)
            fin.close()
            fo.close()

        def write_hash_SVID_wham_SR_DEL_under1Kb_unique_GT(hash_SVID_wham_SR_DEL_under1Kb_unique_GT, fileout):
            fo=open(fileout,'w')
            for i in hash_SVID_wham_SR_DEL_under1Kb_unique_GT.keys():
                for j in hash_SVID_wham_SR_DEL_under1Kb_unique_GT[i]:
                    print('\t'.join([i,j]), file=fo)
            fo.close()

        import os
        import pysam

        hash_SVID_DEL_cluster_MEMBERS = collect_SVID_from_DEL_clusters("~{vcf}")
        hash_SVID_wham_SR_DEL_under1Kb_unique_GT = extract_wham_SR_DEL_under1Kb_from_before_clustered_vcf("~{before_clustered_vcf}", hash_SVID_DEL_cluster_MEMBERS)
        revise_wham_SR_DEL_under1Kb("~{vcf}", "~{prefix}.revised.vcf.gz", hash_SVID_wham_SR_DEL_under1Kb_unique_GT)
        write_hash_SVID_wham_SR_DEL_under1Kb_unique_GT(hash_SVID_wham_SR_DEL_under1Kb_unique_GT,"~{prefix}.SVID_wham_SR_DEL_under1Kb_unique_GT")

        CODE

        tabix -p vcf "~{prefix}.revised.vcf.gz"

    >>>

    output{
        File revised_vcf = "~{prefix}.revised.vcf.gz"
        File revised_vcf_idx = "~{prefix}.revised.vcf.gz.tbi"
        File SVID_wham_SR_DEL_under1Kb_unique_GT = "~{prefix}.SVID_wham_SR_DEL_under1Kb_unique_GT"
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


task SplitSvidByContig{
    input{
        File SVID
        String contig
        String sv_base_mini_docker
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

    String prefix = basename(SVID,".gz")

    command<<<
        set -euo pipefail
        zcat ~{SVID} | grep "_~{contig}_" > ~{prefix}.~{contig}
        bgzip ~{prefix}.~{contig}
    >>>

    output{
        File SVID_split = "~{prefix}.~{contig}.gz"
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

task CombineWithCtx{
    input{
        File vcf
        File vcf_idx
        File? ctx_vcf
        File? ctx_vcf_idx
        String contig
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 200,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    String prefix = basename(vcf,".vcf.gz")

    command<<<
        set -euo pipefail

        bcftools view ~{ctx_vcf} ~{contig} -O z -o ~{contig}.ctx.vcf.gz
        tabix -p vcf ~{contig}.ctx.vcf.gz

        echo "~{vcf}" > vcf_list
        echo "~{contig}.ctx.vcf.gz" >> vcf_list

        bcftools concat -a --allow-overlaps --output-type z --file-list vcf_list --output  "~{prefix}.with_ctx.vcf.gz"
        tabix -p vcf  "~{prefix}.with_ctx.vcf.gz"

    >>>

    output{
        File out_vcf_list = "vcf_list"
        File combined_vcf = "~{prefix}.with_ctx.vcf.gz"
        File combined_vcf_idx = "~{prefix}.with_ctx.vcf.gz.tbi"
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


#this is a patch to fix the SVs that were enriched for 10% DEL outliers but not 5%;
#these SVs will have their filter columns revised to "PASS"
task AdditionalFiltersStep3{
    input{
        File vcf
        File vcf_idx
        File SVID_to_Remove_DEL_bump
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
                    out[pin[0]] = [pin[1]]
                else:
                    out[pin[0]] += [pin[1]]
            fin.close()
            return out

        def SVID_sample_to_remove_readin(SVID_sample_to_remove):
            out={}
            fin=os.popen(r'''zcat %s'''%(SVID_sample_to_remove))
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                    out[pin[0]] = []
                out[pin[0]].append(pin[1])
            fin.close()
            return out

        def revise_vcf(vcf_input, vcf_output, hash_SVID_to_remove):

            fin=pysam.VariantFile(vcf_input)
            header = fin.header
            header.filters.remove_header('MULTIALLELIC')
            if not 'MULTIALLELIC'  in header.info.keys():
                header.info.add('MULTIALLELIC', None, "Flag", "Multiallelic site")

            fo=pysam.VariantFile(vcf_output, 'w', header = header)
            for record in fin:
                if record.id in hash_SVID_to_remove.keys(): 
                    if 'OUTLIER_SAMPLE_ENRICHED_5to10perc' in hash_SVID_to_remove[record.id]:   
                        if "OUTLIER_SAMPLE_ENRICHED" in record.filter.keys():
                            del record.filter['OUTLIER_SAMPLE_ENRICHED']
                        if len(record.filter.keys())==0:
                            record.filter.add("PASS")
                if record.filter.keys()==["MULTIALLELIC"]:
                    record.filter.add("PASS")
                    record.info['MULTIALLELIC'] = True
                fo.write(record)

            fin.close()
            fo.close()


        import os
        import pysam
        hash_SVID_to_remove = SVID_to_remove_readin("~{SVID_to_Remove_DEL_bump}")
        revise_vcf("~{vcf}", "~{prefix}.revised.vcf.gz", hash_SVID_to_remove)

        CODE

        tabix -p vcf "~{prefix}.revised.vcf.gz"

    >>>

    output{
        File revised_vcf = "~{prefix}.revised.vcf.gz"
        File revised_vcf_idx = "~{prefix}.revised.vcf.gz.tbi"
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



