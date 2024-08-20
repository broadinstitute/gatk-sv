version 1.0

import "Structs.wdl"
import "HailMerge.wdl" as HailMerge
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "CollectPEMetricsForCPX.wdl" as collect_pe_metrics_for_cpx
import "CollectLgCnvSupportForCPX.wdl" as collect_lg_cnv_supp_for_cpx

workflow FilterComplex {
    input {
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
        File sample_PE_metrics
        File sample_depth_calls
        File? script_generate_cpx_review_script

        Boolean use_hail = false
        String? gcs_project

        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_collect_pe
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_split_script
        RuntimeAttr? runtime_attr_calcu_pe_stat
        RuntimeAttr? runtime_attr_extract_cpx_ctx
        RuntimeAttr? runtime_attr_concat_evidence
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_preconcat_sharded_cluster
        RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
        RuntimeAttr? runtime_attr_fix_header_sharded_cluster
        RuntimeAttr? runtime_attr_generate_cpx_review_script
        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_extract_cpx_lg_cnv_by_batch
        RuntimeAttr? runtime_attr_seek_depth_supp_for_cpx
        RuntimeAttr? runtime_attr_concat_bed_Step1
        RuntimeAttr? runtime_attr_concat_bed_Step2
        RuntimeAttr? runtime_attr_calcu_cpx_evidences
    }


    call MiniTasks.ScatterVcf {
        input:
          vcf = vcf,
          prefix = prefix,
          records_per_shard = n_per_split,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_scatter_vcf
      }

    scatter (i in range(length(ScatterVcf.shards))) {

        call Vcf2Bed as Vcf2Bed {
            input:
                vcf = ScatterVcf.shards[i],
                vcf_idx = ScatterVcf.shards_idx[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_vcf2bed
        }

        call SplitCpxCtx {
            input:
                bed = Vcf2Bed.bed,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_extract_cpx_ctx
        }

        call collect_lg_cnv_supp_for_cpx.CollectLgCnvSupportForCPX {
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

        call GenerateCpxReviewScript {
            input:
                bed = SplitCpxCtx.cpx_ctx_bed,
                sample_PE_metrics = sample_PE_metrics,
                script = script_generate_cpx_review_script,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_generate_cpx_review_script
        }

        call collect_pe_metrics_for_cpx.CollectPEMetricsForCPX {
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

        call CalculateCpxEvidences {
            input:
                PE_collect_script = GenerateCpxReviewScript.pe_evidence_collection_script,
                PE_supp = CollectPEMetricsForCPX.evi_stat,
                depth_supp = CollectLgCnvSupportForCPX.lg_cnv_depth_supp,
                prefix = "~{prefix}.~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_calcu_cpx_evidences
        }

        call ReviseVcf {
            input:
                vcf_file = ScatterVcf.shards[i],
                vcf_index = ScatterVcf.shards_idx[i],
                CPX_manual = CalculateCpxEvidences.manual_revise_CPX_results,
                prefix = "~{prefix}.~{i}",
                sv_pipeline_docker = sv_pipeline_docker
        }
    }

    if (length(ReviseVcf.revised_vcf) == 0) {
        call MiniTasks.GetVcfHeaderWithMembersInfoLine as GetVcfHeader {
            input:
                vcf_gz=vcf,
                prefix="~{prefix}.manual_review",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_vcf_header_with_members_info_line
        }
    }

    if (length(ReviseVcf.revised_vcf) > 0) {
        if (use_hail) {
            call HailMerge.HailMerge as ConcatVcfsHail_annotated {
                input:
                    vcfs=ReviseVcf.revised_vcf,
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
            call MiniTasks.ConcatVcfs {
                input:
                    vcfs=ReviseVcf.revised_vcf,
                    vcfs_idx=ReviseVcf.revised_vcf_idx,
                    allow_overlaps=true,
                    outfile_prefix="~{prefix}.manual_review",
                    sv_base_mini_docker=sv_base_mini_docker,
                    runtime_attr_override=runtime_attr_concat_sharded_cluster
            }
        }
    }

    call MiniTasks.ConcatFiles as concat_cpx_evidences {
        input:
            shard_bed_files = CalculateCpxEvidences.manual_revise_CPX_results,
            prefix = "~{prefix}.CPX_evidence" ,
            sv_base_mini_docker = sv_base_mini_docker

    }

    File reviewed_vcf = select_first([GetVcfHeader.out, ConcatVcfs.concat_vcf, ConcatVcfsHail_annotated.merged_vcf])
    File reviewed_vcf_idx = select_first([GetVcfHeader.out_idx, ConcatVcfs.concat_vcf_idx, ConcatVcfsHail_annotated.merged_vcf_index])

    output {
        File revised_output_vcf = reviewed_vcf
        File revised_output_vcf_idx = reviewed_vcf_idx
        File cpx_evidences = concat_cpx_evidences.merged_file
    }
}


task ReviseVcf {
    input {
        File vcf_file
        File vcf_index
        File CPX_manual
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 7.5,
        disk_gb: ceil(5.0 +  size(vcf_file, "GB")*3),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail

        python <<CODE

        def recal_qual_score(record):
            """
            Recalibrate quality score for a single variant
            """
            quals = []
            for s in [s for s in record.samples]:
                GT = record.samples[s]['GT']
                if GT in NULL_and_REF_GTs:
                    continue
                elif GT in HET_GTs:
                    quals.append(record.samples[s]['GQ'])
                else:
                    quals.append(99)

            if len(quals) > 0:
                return int(median(quals))

        def CPX_manual_readin(CPX_manual):
            out={}
            fin=open(CPX_manual)
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                    out[pin[0]] = {}
                out[pin[0]][pin[1]] = pin[2:]
            fin.close()
            return out

        def revise_vcf(vcf_input, vcf_output, hash_CPX_manual):
            fin=pysam.VariantFile(vcf_input)
            #revise vcf header
            header = fin.header
            fo=pysam.VariantFile(vcf_output, 'w', header = header)
            for record in fin:
                #label CPX with manual review results:
                if record.id in hash_CPX_manual.keys():
                    unresolve_rec = 0
                    for sample in hash_CPX_manual[record.id].keys():
                        if sample in record.samples.keys():
                            if hash_CPX_manual[record.id][sample][0] == 'no_PE':
                                if hash_CPX_manual[record.id][sample][1] in ["NA","lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                            elif hash_CPX_manual[record.id][sample][0] == 'low_PE':
                               if hash_CPX_manual[record.id][sample][1] in ["NA","lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                            elif hash_CPX_manual[record.id][sample][0] == 'partial_PE':
                               if hash_CPX_manual[record.id][sample][1] in ["NA","lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                                    unresolve_rec+=1
                            elif hash_CPX_manual[record.id][sample][0] == 'high_PE':
                               if hash_CPX_manual[record.id][sample][1] in ["lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                    if not unresolve_rec/len(hash_CPX_manual[record.id].keys())<.5:
                        if 'PASS' in record.filter:
                            record.filter.clear()
                        record.filter.add('UNRESOLVED')
                    #Recalibrate QUAL score after CPX filtering
                    newQUAL = recal_qual_score(record)
                    if newQUAL is not None:
                        record.qual = newQUAL
                fo.write(record) # write out every record that was in the input - NCR will remove ones with no carriers left
            fin.close()
            fo.close()

        import os
        import sys
        from numpy import median
        import pysam
        import argparse

        #Define global variables
        NULL_GTs = [(None, None), (None, )]
        REF_GTs = [(0, 0), (0, ), (None, 2)]
        NULL_and_REF_GTs = NULL_GTs + REF_GTs
        HET_GTs = [(0, 1), (None, 1), (None, 3)]

        hash_CPX_manual =  CPX_manual_readin("~{CPX_manual}")
        print(len(hash_CPX_manual.keys()))
        revise_vcf("~{vcf_file}", "~{prefix}.Manual_Revised.vcf.gz", hash_CPX_manual)

CODE

        tabix ~{prefix}.Manual_Revised.vcf.gz
    >>>

    output {
        File revised_vcf = "~{prefix}.Manual_Revised.vcf.gz"
        File revised_vcf_idx = "~{prefix}.Manual_Revised.vcf.gz.tbi"
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

task Vcf2Bed {
    input {
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
    command <<<
        set -euo pipefail
        svtk vcf2bed -i ALL --include-filters ~{vcf} ~{prefix}.bed
        bgzip ~{prefix}.bed
    >>>

    output {
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

task SplitCpxCtx {
    input {
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

    command <<<
        set -eu

        zcat ~{bed} | head -1 > ~{prefix}.cpx_ctx.bed

        set -o pipefail

        filterColumn=$(zcat ~{bed} | head -1 | tr "\t" "\n" | awk '$1=="FILTER" {print NR}')

        zcat ~{bed} | awk 'NR > 1' | { grep CPX || true; } | awk -v filter_column=${filterColumn} '$filter_column !~ /UNRESOLVED/' >> ~{prefix}.cpx_ctx.bed

        bgzip ~{prefix}.cpx_ctx.bed
    >>>

    output {
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
    input {
        File bed
        File sample_PE_metrics
        File? script
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

    command <<<
        set -euo pipefail

        python ~{default="/opt/sv-pipeline/scripts/manual_review/reformat_CPX_bed_and_generate_script.py" script} \
        -i ~{bed} \
        -s ~{sample_PE_metrics} \
        -p CPX_CTX_disINS.PASS.PE_evidences \
        -c collect_PE_evidences.CPX_CTX_disINS.PASS.sh \
        -r ~{prefix}.svelter

    >>>

    output {
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
    input {
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

    output {
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

task CalculateCpxEvidences{
    input {
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

    command <<<

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
                            if PE_supp[svid][sample][0][0]>1 and PE_supp[svid][sample][1][0]>1:
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

    output {
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
