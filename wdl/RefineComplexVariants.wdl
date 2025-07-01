version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "CollectPEMetricsForCPX.wdl" as collect_pe_metrics_for_cpx
import "CollectLargeCNVSupportForCPX.wdl" as collect_lg_cnv_supp_for_cpx
import "Utils.wdl" as util

workflow RefineComplexVariants {
    input {
        File vcf
        String prefix

        Array[String] batch_name_list
        Array[File] batch_sample_lists
        Array[File] PE_metrics
        Array[File] PE_metrics_indexes
        Array[File] Depth_DEL_beds
        Array[File] Depth_DUP_beds

        Int n_per_split
        Int min_pe_cpx = 3
        Int min_pe_ctx = 3

        String sv_base_mini_docker
        String sv_pipeline_docker
        String linux_docker

        RuntimeAttr? runtime_attr_sample_batch
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_collect_pe
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_split_script
        RuntimeAttr? runtime_attr_calcu_pe_stat
        RuntimeAttr? runtime_attr_split_cpx_ctx
        RuntimeAttr? runtime_attr_concat_evidence
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_preconcat
        RuntimeAttr? runtime_attr_fix_header
        RuntimeAttr? runtime_attr_generate_cpx_review_script
        RuntimeAttr? runtime_attr_generate_cnv_segments_from_cpx
        RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
        RuntimeAttr? runtime_attr_extract_cpx_lg_cnv_by_batch
        RuntimeAttr? runtime_attr_seek_depth_supp_for_cpx
        RuntimeAttr? runtime_attr_concat_bed_Step1
        RuntimeAttr? runtime_attr_concat_bed_Step2
        RuntimeAttr? runtime_attr_calcu_cpx_evidences
    }

    call GetSampleBatchPEMap {
        input:
            batch_name_list = batch_name_list,
            batch_sample_lists = batch_sample_lists,
            batch_pe_files = write_lines(PE_metrics),
            prefix = prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_sample_batch
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

        call util.VcfToBed {
            input:
                vcf_file = ScatterVcf.shards[i],
                args = "-i ALL --include-filters",
                variant_interpretation_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_vcf2bed
        }

        call SplitCpxCtx {
            input:
                bed = VcfToBed.bed_output,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_split_cpx_ctx
        }

        call collect_lg_cnv_supp_for_cpx.CollectLargeCNVSupportForCPX {
            input:
                cpx_ctx_bed = SplitCpxCtx.cpx_ctx_bed,
                prefix = "~{prefix}.~{i}",

                sample_batch_pe_map = GetSampleBatchPEMap.sample_batch_pe_map,
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
                sample_batch_pe_map = GetSampleBatchPEMap.sample_batch_pe_map,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_generate_cpx_review_script
        }

        call collect_pe_metrics_for_cpx.CollectPEMetricsForCPX {
            input:
                batch_name_list = batch_name_list,
                PE_metrics = PE_metrics,
                PE_metrics_idxes = PE_metrics_indexes,
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
                depth_supp = CollectLargeCNVSupportForCPX.lg_cnv_depth_supp,
                min_pe_cpx = min_pe_cpx,
                prefix = "~{prefix}.~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_calcu_cpx_evidences
        }

        call CalculateCtxEvidences{
            input:
                PE_collect_script = GenerateCpxReviewScript.pe_evidence_collection_script,
                PE_supp = CollectPEMetricsForCPX.evi_stat,
                min_pe_ctx = min_pe_ctx,
                prefix = "~{prefix}.~{i}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_calcu_cpx_evidences
        }

        call ReviseVcf {
            input:
                vcf_file = ScatterVcf.shards[i],
                CPX_manual = CalculateCpxEvidences.manual_revise_CPX_results,
                CTX_manual = CalculateCtxEvidences.manual_revise_CTX_results,
                unresolved_svids = GenerateCpxReviewScript.unresolved_svids,
                prefix = "~{prefix}.~{i}",
                sv_pipeline_docker = sv_pipeline_docker
        }
    }

    call MiniTasks.ConcatVcfs {
        input:
            vcfs=ReviseVcf.revised_vcf,
            vcfs_idx=ReviseVcf.revised_vcf_idx,
            allow_overlaps=true,
            outfile_prefix="~{prefix}.cpx_refined",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_concat
    }

    call MiniTasks.ConcatHeaderedTextFiles {
        input:
            text_files = CalculateCpxEvidences.manual_revise_CPX_results,
            output_filename = "~{prefix}.CPX_evidence.txt" ,
            linux_docker = linux_docker

    }

    output {
        File cpx_refined_vcf = ConcatVcfs.concat_vcf
        File cpx_refined_vcf_index = ConcatVcfs.concat_vcf_idx
        File cpx_evidences = ConcatHeaderedTextFiles.out
    }
}


task GetSampleBatchPEMap {
    input {
        Array[File] batch_sample_lists
        Array[String] batch_name_list
        File batch_pe_files
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10 + 2*ceil(size(flatten([batch_sample_lists]), "GB")),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
        set -euo pipefail

        python <<CODE
import os
batch_sample_lists = ["~{sep='", "' batch_sample_lists}"]
batch_name_list = ["~{sep='", "' batch_name_list}"]
batch_pe_files = []
with open("~{batch_pe_files}", 'r') as pe:
    for line in pe:
        local_file = os.path.basename(line.strip())
        batch_pe_files.append(local_file)
with open("~{prefix}.sample_batch_pe_map.tsv", 'w') as out:
    for i in range(len(batch_name_list)):
        with open(batch_sample_lists[i], 'r') as inp:
            for line in inp:
                sample = line.strip()
                batch = batch_name_list[i]
                pe_file = batch_pe_files[i]
                out.write(f"{sample}\t{batch}\t{pe_file}\n")
CODE

    >>>

    output {
        File sample_batch_pe_map = "~{prefix}.sample_batch_pe_map.tsv"
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


task ReviseVcf {
    input {
        File vcf_file
        File CPX_manual
        File CTX_manual
        File unresolved_svids
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: ceil(5.0 +  size(vcf_file, "GB")*3),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail

        python <<CODE

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


def CTX_manual_readin(CTX_manual):
    fin=open(CTX_manual)
    out={}
    for line in fin:
        pin=line.strip().split('\t')
        if not pin[0] in out.keys():
            out[pin[0]] = {}
        if not pin[1] in out[pin[0]].keys():
            out[pin[0]][pin[1]] = pin[2:]
    fin.close()
    return out


def unresolved_readin(unresolved_svids):
    svids = set()
    with open(unresolved_svids, 'r') as inp:
        for line in inp:
            svids.add(line.strip())
    return svids


def revise_vcf(vcf_input, vcf_output, hash_CPX_manual, unresolved_svids, hash_CTX_manual):
    fin=pysam.VariantFile(vcf_input)
    #revise vcf header
    header = fin.header
    fo=pysam.VariantFile(vcf_output, 'w', header = header)
    for record in fin:
        if record.id in unresolved_svids:
            if 'PASS' in record.filter:
                record.filter.clear()
            record.filter.add('UNRESOLVED')
        #label CPX with manual review results:
        elif record.id in hash_CPX_manual.keys():
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
        #label CTX with manual review results:
        elif record.id in hash_CTX_manual.keys():
            if 'NA' in hash_CTX_manual[record.id].keys():
                record.filter.add('UNRESOLVED')
            else:
                for sample in hash_CTX_manual[record.id].keys():
                    if sample in record.samples.keys():
                        if hash_CTX_manual[record.id][sample][0] in ['no_PE', 'low_PE', 'partial_PE']:
                            record.samples[sample]['GT'] = [None,None]
        # revisions to insertion-type complex events
        if record.info['SVTYPE']=="CPX" and record.info['CPX_TYPE'] in ['dDUP','dDUP_iDEL','INS_iDEL']:
            if record.info['CPX_TYPE']=='INS_iDEL':
                record.info['CPX_INTERVALS'] = ','.join([x for x in record.info['CPX_INTERVALS']] + [record.info['SOURCE']])
            del record.info['SOURCE']
        elif record.info['SVTYPE']=="INS" and 'SOURCE' in record.info.keys() and "INV" in record.info['SOURCE']: 
            record.info['SVTYPE']="CPX"
            record.alts = ('<CPX>',)
            del_section = record.stop - record.pos
            if del_section<50:
                record.info['CPX_TYPE']="dDUP"
                record.info['CPX_INTERVALS'] = record.info['SOURCE'].replace('INV','DUP')+','+record.info['SOURCE']
            else:
                record.info['CPX_TYPE']="dDUP_iDEL"
                record.info['CPX_INTERVALS'] = record.info['SOURCE'].replace('INV','DUP')+','+record.info['SOURCE']+','+"DEL_"+record.chrom+":"+str(record.pos)+'-'+str(record.stop)
            del record.info['SOURCE']
        fo.write(record) # write out every record that was in the input - NCR will remove ones with no carriers left
    fin.close()
    fo.close()

import pysam

unresolved_svids = unresolved_readin("~{unresolved_svids}")
hash_CPX_manual =  CPX_manual_readin("~{CPX_manual}")
hash_CTX_manual = CTX_manual_readin("~{CTX_manual}")
print(len(hash_CPX_manual.keys()))
revise_vcf("~{vcf_file}", "~{prefix}.Manual_Revised.vcf.gz", hash_CPX_manual, unresolved_svids, hash_CTX_manual)

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

task SplitCpxCtx {
    input {
        File bed
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(bed, ".bed.gz")

    command <<<
        set -eu

        zcat ~{bed} | head -1 > ~{prefix}.cpx_ctx.bed

        filterColumn=$(zcat ~{bed} | head -1 | tr "\t" "\n" | awk '$1=="FILTER" {print NR}')

        set -o pipefail

        zcat ~{bed} | awk 'NR > 1' | { grep CPX || true; } | awk -v filter_column=${filterColumn} '$filter_column !~ /UNRESOLVED/' >> ~{prefix}.cpx_ctx.bed

        zcat ~{bed} | awk 'NR > 1' | { grep CTX || true; } >> ~{prefix}.cpx_ctx.bed

        # INS with INV in SOURCE - will be converted to CPX later so need to evaluate evidence
        zcat ~{bed} | awk 'NR > 1' | { grep INS || true; } | { grep INV || true; } >> ~{prefix}.cpx_ctx.bed

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

task GenerateCpxReviewScript {
    input {
        File bed
        File sample_batch_pe_map
        File? script
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(bed, ".bed.gz")

    command <<<
        set -euo pipefail

        cut -f1,3 ~{sample_batch_pe_map} > sample_PE_metrics.tsv

        bgzip -cd '~{bed}' \
          | awk -F'\t' '/^#/{print; next} {$2++; print}' OFS='\t' - \
          | bgzip -c > shifted.bed.gz

        python ~{default="/opt/sv-pipeline/scripts/reformat_CPX_bed_and_generate_script.py" script} \
        -i shifted.bed.gz \
        -s sample_PE_metrics.tsv \
        -p CPX_CTX_disINS.PASS.PE_evidences \
        -c collect_PE_evidences.CPX_CTX_disINS.PASS.sh \
        -r ~{prefix}.svelter \
        -u ~{prefix}.unresolved_svids.txt

    >>>

    output {
        File pe_evidence_collection_script = "collect_PE_evidences.CPX_CTX_disINS.PASS.sh"
        File svelter = "~{prefix}.svelter"
        File unresolved_svids = "~{prefix}.unresolved_svids.txt"
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

task CalculateCpxEvidences {
    input {
        File PE_collect_script
        File PE_supp
        File depth_supp
        Int min_pe_cpx
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
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
                            if PE_supp[svid][sample][0][0]>=~{min_pe_cpx} and PE_supp[svid][sample][1][0]>=~{min_pe_cpx}:
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

task CalculateCtxEvidences {
    input {
        File PE_collect_script
        File PE_supp
        Int min_pe_ctx
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 3.75,
        disk_gb: 10,
        boot_disk_gb: 10,
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
                            if PE_supp[svid][sample][0][0]>=~{min_pe_ctx} and PE_supp[svid][sample][1][0]>=~{min_pe_ctx}:
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

    output {
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

