version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge

# Workflow to shard a filtered vcf & run vcfcluster (sub-sub-sub workflow)
workflow ManualRevise {
  input {
    File vcf
    String prefix
    File? exclude_list
    Float merging_shard_scale_factor = 30000000

    File SVID_to_Remove
    File MEI_DEL_Rescue
    File CPX_manual
    File CTX_manual

    Boolean use_hail
    String? gcs_project

    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_base_mini_docker
    String sv_benchmark_docker

    # overrides for local tasks
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_revise_vcf 
    RuntimeAttr? runtime_attr_shard_clusters
    RuntimeAttr? runtime_attr_shard_vids
    RuntimeAttr? runtime_attr_pull_vcf_shard
    RuntimeAttr? runtime_attr_svtk_vcf_cluster
    RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line

    RuntimeAttr? runtime_attr_preconcat_sharded_cluster
    RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
    RuntimeAttr? runtime_attr_fix_header_sharded_cluster

    # overrides for merge subworkflow
    RuntimeAttr? runtime_attr_merge_clusters
    RuntimeAttr? runtime_attr_concat_inner_shards

    # overrides for MiniTasks
    RuntimeAttr? runtime_attr_concat_sharded_cluster
    RuntimeAttr? runtime_attr_sort_merged_vcf
    RuntimeAttr? runtime_attr_count_samples
    RuntimeAttr? runtime_attr_get_vids
    RuntimeAttr? runtime_attr_cat_vid_lists_sharded
    RuntimeAttr? runtime_attr_make_sites_only
    RuntimeAttr? runtime_attr_clean_up_formats_revised_vcf
  }


  File vcf_idx = vcf + ".tbi"
  if (defined(exclude_list)) {
    File exclude_list_idx = exclude_list + ".tbi"
  }

  #Run vcfcluster per shard
  call ReviseVcf{
      input:
          vcf_file = vcf,
          vcf_index = vcf_idx,
          SVID_to_Remove = SVID_to_Remove,
          MEI_DEL_Rescue = MEI_DEL_Rescue,
          CPX_manual = CPX_manual,
          CTX_manual = CTX_manual,
          sv_benchmark_docker = sv_benchmark_docker,
          runtime_attr_override = runtime_attr_revise_vcf
  }

  call CleanUpFormats as clean_up_formats_revised_vcf{
    input:
      vcf_file = ReviseVcf.out_manual_revised_vcf,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_clean_up_formats_revised_vcf
  }

  call CleanUpFormats as clean_up_formats_cpx_ctx_vcf{
    input:
      vcf_file = ReviseVcf.out_cpx_ctx_vcf,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_clean_up_formats_revised_vcf
  }

  call MiniTasks.SortVcf as sort_manual_revised_vcf{
    input:
      vcf = clean_up_formats_revised_vcf.reformatted_vcf,
      outfile_prefix = "~{prefix}.manual_revised.sorted.shard",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_sort_merged_vcf
  }

  call MiniTasks.SortVcf as sort_cpx_ctx_vcf{
    input:
      vcf = clean_up_formats_cpx_ctx_vcf.reformatted_vcf,
      outfile_prefix = "~{prefix}.cpx_ctx.sorted.shard",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_sort_merged_vcf
  }


  #Output
  output {
    File manual_revised_vcf = sort_manual_revised_vcf.out
    File manual_revised_vcf_idx = sort_manual_revised_vcf.out_index
    File cpx_ctx_vcf = sort_cpx_ctx_vcf.out
    File cpx_ctx_vcf_idx = sort_cpx_ctx_vcf.out_index
    }
}

# Adds MEMBERS definition to header (workaround for when VIDs_list is empty)
task GetVcfHeaderWithMembersInfoLine {
  input {
    File vcf_gz
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr runtime_default = object {
    mem_gb: 1,
    disk_gb: 10,
    cpu_cores: 1,
    preemptible_tries: 3,
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

  command <<<
    set -euxo pipefail
    zgrep "^##" ~{vcf_gz} > ~{prefix}.vcf
    echo "##INFO=<ID=MEMBERS,Number=.,Type=String,Description=\"IDs of cluster's constituent records.\">" >> ~{prefix}.vcf
    zgrep "^#CHROM" ~{vcf_gz} >> ~{prefix}.vcf
    bgzip ~{prefix}.vcf
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_idx = "~{prefix}.vcf.gz.tbi"
  }
}

#Do fast cluster on vcf (sample_overlap = 0) to generate shards
task ShardClusters {
  input {
    File vcf
    String prefix
    Int dist
    Float frac
    File? exclude_list
    File? exclude_list_idx
    Int svsize
    Array[String] sv_types
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 1.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
    cpu_cores: 1,
    preemptible_tries: 3,
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
    svtk vcfcluster <(echo "~{vcf}") ~{prefix}.vcf.gz \
      -d ~{dist} \
      -f ~{frac} \
      ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
      -z ~{svsize} \
      -p ~{prefix} \
      -t ~{sep=',' sv_types} \
      -o 0 \
      --preserve-header \
      --preserve-ids \
      --skip-merge
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}

task CleanUpFormats{
    input{
      File vcf_file
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

    String prefix  = basename(vcf_file, ".vcf.gz")

    command<<<
        set -euo pipefail

        tabix -p vcf ~{vcf_file}
        python <<CODE

        import os
        import pysam
        
        formats_to_keep = ["GT","GQ","EV","CN","CNQ","RD_CN","RD_GQ","PE_GT","PE_GQ","SR_GT","SR_GQ"]

        fin=pysam.VariantFile("~{vcf_file}")
        header = fin.header
        for item in header.formats.keys():
          if not item in formats_to_keep:
            header.formats.remove_header(item)

        fo=pysam.VariantFile("test.vcf.gz", 'w', header = header)
        for record in fin:
          if record.info['SVTYPE']=="CPX" and record.info['CPX_TYPE'] in ['dDUP','dDUP_iDEL','INS_iDEL']:
            if record.info['CPX_TYPE']=='INS_iDEL':
              record.info['CPX_INTERVALS'] = ','.join([x for x in record.info['CPX_INTERVALS']] + [record.info['SOURCE']])
              del record.info['SOURCE']
            else:
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
          for i in record.format.keys():
            if not i in formats_to_keep:
              del(record.format[i])
          fo.write(record)
        fo.close()
        CODE

        zcat test.vcf.gz | grep -v "##bcftools_" | grep -v "##hailversion=" | grep -v "##source=cleanvcf" | bgzip > "~{prefix}.reformatted.vcf.gz"
        tabix -p vcf "~{prefix}.reformatted.vcf.gz"
    >>>

    output{
        File reformatted_vcf = "~{prefix}.reformatted.vcf.gz"
        File reformatted_vcf_idx = "~{prefix}.reformatted.vcf.gz.tbi"
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

task ReviseVcf{
    input{
        File vcf_file
        File vcf_index
        File SVID_to_Remove
        File MEI_DEL_Rescue
        File CPX_manual
        File CTX_manual
        String sv_benchmark_docker
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


    String prefix  = basename(vcf_file, ".vcf.gz")
    command<<<
        set -euo pipefail

        #python /src/revise_vcf_with_manual_review_results.py ~{vcf_file} ~{prefix}.Manual_Revised.vcf.gz --cpx_vcf ~{prefix}.CPX_CTX.vcf.gz --SVID_to_remove  ~{SVID_to_Remove} --MEI_DEL_rescue  ~{MEI_DEL_Rescue} --CPX_manual  ~{CPX_manual} --CTX_manual  ~{CTX_manual}

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
                    quals.append(999)

            if len(quals) > 0:
                return int(median(quals))

        def SVID_to_remove_readin(SVID_to_remove):
            out={}
            fin=open(SVID_to_remove)
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                    out[pin[0]] = pin[1]
            fin.close()
            return out

        def SVID_MEI_DEL_readin(MEI_DEL_reset):
            out={}
            fin=open(MEI_DEL_reset)
            for line in fin:
                pin=line.strip().split()
                if not pin[0] in out.keys():
                    out[pin[0]] = pin[3]
            fin.close()
            return out

        def SVID_duplicated_readin(duplicated_SVID_manual):
            duplicated_SVID_list = []
            duplicated_member_merge={}
            duplicated_info = {}
            fin=open(duplicated_SVID_manual)
            for line in fin:
                pin=line.strip().split()
                duplicated_SVID_list+=pin[1:3]
                duplicated_member_merge[pin[1]] = pin[0]
                duplicated_member_merge[pin[2]] = pin[0]
                duplicated_info[pin[0]] = pin[3:]
            fin.close()
            return([duplicated_member_merge, duplicated_info])

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

        def select_filter(filters):
            filter_uni = []
            for i in filters.split(','):
                if not i in filter_uni:
                    filter_uni.append(i)
            if len(filter_uni) == 1:
                out = filter_uni[0]
            elif 'PASS' in filter_uni:
                out = 'PASS'
            elif 'UNRESOLVED' in filter_uni:
                out = 'UNRESOLVED'
            elif 'HIGH_PCRMINUS_NOCALL_RATE' in filter_uni:
                out = 'HIGH_PCRMINUS_NOCALL_RATE'
            elif 'HIGH_PCRPLUS_NOCALL_RATE' in filter_uni:
                out = 'HIGH_PCRPLUS_NOCALL_RATE'
            elif 'PESR_GT_OVERDISPERSION' in filter_uni:
                out = 'PESR_GT_OVERDISPERSION'
            elif 'BOTHSIDES_SUPPORT' in filter_uni:
                out = 'BOTHSIDES_SUPPORT'
            else:
                out = 'LowQual'
            return out

        def sample_merge(record1, record2):
            filter = select_filter(','.join(record1.filter.keys() + record2.filter.keys()))
            record1.filter.add(filter)
            for i in record1.samples.keys():
                if record1.samples[i]['GT'] != record2.samples[i]['GT']:
                    if record1.samples[i]['GT']==(0, 0): 
                        for j in record1.samples[i].keys():
                            record1.samples[i][j] = record2.samples[i][j]
            return record1

        def revise_vcf(vcf_input, vcf_output, ctx_output, hash_SVID_to_remove, hash_MEI_DEL_reset, hash_CPX_manual, hash_CTX_manual):
            fin=pysam.VariantFile(vcf_input)
            #revise vcf header
            header = fin.header
            for filt in filts_to_remove:
                if filt in header.filters.keys():
                    header.filters.remove_header(filt)
            header.filters.add('REDUNDANT_LG_CNV', None, None, "Multiple large CNVs called at the same locus likely indicates unreliable clustering and/or low-quality multiallelic locus")
            header.filters.add('FAIL_MANUAL_REVIEW', None, None, "Low-quality variant that did not pass manual review of supporting evidence")
            header.filters.add('OUTLIER_SAMPLE_ENRICHED',None, None, "Deletion is enriched for non-reference genotypes in outlier samples, likely indicating noisy or unreliable genotypes")
            header.filters.remove_header('BOTHSIDES_SUPPORT')
            header.filters.remove_header('HIGH_SR_BACKGROUND')
            header.filters.remove_header('PESR_GT_OVERDISPERSION')
            header.formats.add('MANUAL', '1', 'String', 'Reason for a failure from manual review')
            NEW_INFOS = ['##INFO=<ID=PESR_GT_OVERDISPERSION,Number=0,Type=Flag,Description="PESR genotyping data is overdispersed. Flags sites where genotypes are likely noisier.">',
                         '##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description="Suspicious accumulation of split reads in predicted non-carrier samples. Flags sites more prone to false discoveries and where breakpoint precision is reduced.">',
                         '##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description="Variant has read-level support for both sides of breakpoint.Indicates higher-confidence variants.">']
            for info in NEW_INFOS:
                header.add_line(info)
            fo=pysam.VariantFile(vcf_output, 'w', header = header)
            fo2 = pysam.VariantFile(ctx_output, 'w', header = header)
            for record in fin:
                #revise filter columns:
                for filt in filts_for_info:
                    if filt in record.filter:
                        record.info[filt] = True
                newfilts = [filt for filt in record.filter if filt not in filts_to_remove]
                record.filter.clear()
                for filt in newfilts:
                    record.filter.add(filt)
                if len(record.filter) == 0:
                    record.filter.add('PASS')
                #Recalibrate QUAL score
                newQUAL = recal_qual_score(record)
                if newQUAL is not None:
                    record.qual = newQUAL
                ###remove mCNVs <5Kb, large bi-allelic CNVS that failed manual review, and DEL bump (400bp-1Kb, AF<1, AC>1) enriched for outlier samples
                if record.id in hash_SVID_to_remove.keys(): 
                    if hash_SVID_to_remove[record.id]=='mCNV_under_5Kb':
                        continue
                    elif hash_SVID_to_remove[record.id]=='lg_CNV_failed_manua_review':
                        record.filter.add('FAIL_MANUAL_REVIEW')
                    elif hash_SVID_to_remove[record.id]=="OUTLIER_SAMPLE_ENRICHED":
                        record.filter.add('OUTLIER_SAMPLE_ENRICHED')
                #rescue MEI DELs
                elif record.id in hash_MEI_DEL_reset.keys(): 
                    record.filter.add('PASS')
                    record.info['SVTYPE'] = 'DEL'
                    if hash_MEI_DEL_reset[record.id] == 'overlap_LINE1':
                        record.alts = ('<DEL:ME:LINE1>',)
                    if hash_MEI_DEL_reset[record.id] == 'overlap_HERVK':
                        record.alts = ('<DEL:ME:HERVK>',)
                #label CPX with manual review results:
                elif record.id in hash_CPX_manual.keys():
                    unresolve_rec = 0
                    for sample in hash_CPX_manual[record.id].keys():
                        if sample in record.samples.keys():
                            if hash_CPX_manual[record.id][sample][0] == 'no_PE':
                                if hash_CPX_manual[record.id][sample][1] in ["NA","lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'NO_PE'
                            elif hash_CPX_manual[record.id][sample][0] == 'low_PE':
                               if hash_CPX_manual[record.id][sample][1] in ["NA","lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'LOW_PE'
                            elif hash_CPX_manual[record.id][sample][0] == 'partial_PE':
                               if hash_CPX_manual[record.id][sample][1] in ["NA","lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'PARTIAL_PE'
                                    unresolve_rec+=1
                            elif hash_CPX_manual[record.id][sample][0] == 'high_PE':
                               if hash_CPX_manual[record.id][sample][1] in ["lack_depth","lack_depth,lack_depth", "lack_depth,depth","depth,lack_depth"]:
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'LACK_DEPTH'
                    if not unresolve_rec/len(hash_CPX_manual[record.id].keys())<.5:
                        record.filter.add('UNRESOLVED')
                    else:
                        record.filter.add('PASS')
                    fo2.write(record)
                    continue
                #label CTX with manual review results:
                elif record.id in hash_CTX_manual.keys():
                    if 'NA' in hash_CTX_manual[record.id].keys():
                        record.filter.add('UNRESOLVED')
                    else:
                        for sample in hash_CTX_manual[record.id].keys():
                            if sample in record.samples.keys():
                                if hash_CTX_manual[record.id][sample][1] == 'Not_Enough_PE_Pairs':
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'NO_PE'
                                elif hash_CTX_manual[record.id][sample][1] == 'PARTIAL_PE':
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'PARTIAL_PE'
                                elif hash_CTX_manual[record.id][sample][1] == 'unbalanced_paternal_transmission_of_tloc,2nd_junction_is_missing':
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'PARTIAL_PE'
                                elif hash_CTX_manual[record.id][sample][1] == 'NON_SPECIFIC':
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'NON_SPECIFIC'
                                    record.filter.add('FAIL_MANUAL_REVIEW')
                                elif hash_CTX_manual[record.id][sample][1] == 'Not_Enough_Coordinate_Span':
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'UNDER_DISPERSED_PE'
                                elif hash_CTX_manual[record.id][sample][1] == 'Interrupted_PE_Patterns_when_sorted_on_2nd_breakpoint':
                                    record.samples[sample]['GT'] = [None,None]
                                    record.samples[sample]['MANUAL'] = 'INTERRUPTED_PE'
                                elif hash_CTX_manual[record.id][sample][1] == 'CTX_with_DEL':
                                    record.stop = int(hash_CTX_manual[record.id][sample][0].split(':')[1].split('-')[1])
                    fo2.write(record)
                    continue
                ref_count = len([s for s in record.samples if record.samples[s]['GT'] in NULL_and_REF_GTs])
                alt_count = len(record.samples) - ref_count
                if alt_count>0 or record.info['SVTYPE']=="CNV":
                    fo.write(record)
            fin.close()
            fo.close()
            fo2.close()

        import os
        import sys
        from numpy import median
        import pysam
        import argparse

        #Define global variables
        filts_for_info = 'PESR_GT_OVERDISPERSION HIGH_SR_BACKGROUND BOTHSIDES_SUPPORT VARIABLE_ACROSS_BATCHES'.split(' ')
        filts_to_remove = 'HIGH_PCRPLUS_NOCALL_RATE HIGH_PCRMINUS_NOCALL_RATE'.split(' ')
        filts_to_remove = filts_to_remove + filts_for_info
        NULL_GTs = [(None, None), (None, )]
        REF_GTs = [(0, 0), (0, ), (None, 2)]
        NULL_and_REF_GTs = NULL_GTs + REF_GTs
        HET_GTs = [(0, 1), (None, 1), (None, 3)]

        hash_SVID_to_remove = SVID_to_remove_readin("~{SVID_to_Remove}")
        print(len(hash_SVID_to_remove.keys()))
        hash_MEI_DEL_reset = SVID_MEI_DEL_readin("~{MEI_DEL_Rescue}")
        print(len(hash_MEI_DEL_reset.keys()))
        hash_CPX_manual =  CPX_manual_readin("~{CPX_manual}")
        print(len(hash_CPX_manual.keys()))
        hash_CTX_manual =  CTX_manual_readin("~{CTX_manual}")
        print(len(hash_CTX_manual.keys()))
        revise_vcf("~{vcf_file}", "~{prefix}.Manual_Revised.vcf.gz", "~{prefix}.CPX_CTX.vcf.gz", hash_SVID_to_remove, hash_MEI_DEL_reset, hash_CPX_manual, hash_CTX_manual)

        CODE

    >>>

    output{
        File out_manual_revised_vcf = "~{prefix}.Manual_Revised.vcf.gz"
        File out_cpx_ctx_vcf = "~{prefix}.CPX_CTX.vcf.gz"
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

task FixEndsRescaleGQ {
  input {
    File vcf
    String prefix

    Boolean? fix_ends
    Boolean? rescale_gq

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String outfile = "~{prefix}.vcf.gz"
  Boolean fix_ends_ = select_first([fix_ends, true])
  Boolean rescale_gq_ = select_first([rescale_gq, true])

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<

    set -euo pipefail

    python <<CODE
    import pysam
    import argparse
    from math import floor


    GQ_FIELDS = ["GQ", "PE_GQ", "SR_GQ", "RD_GQ"]


    def fix_bad_end(record):
      # pysam converts to 0-based half-open intervals by subtracting 1 from start, but END is unaltered
      if record.info["SVTYPE"] == "BND"     
        if record.stop < record.start + 2:
          record.info["END2"] = record.stop  # just in case it is not already set. not needed for INS or CPX
        record.stop = record.start + 1
      if record.info["SVTYPE"] == "CTX":  
        record.stop = record.start + 1
      if record.info["SVTYPE"] in ['INS','CPX']:
        if record.stop < record.start + 2:
          record.stop = record.start + 1
        if "END2" in record.info.keys():
          del record.info["END2"]


    def rescale_gq(record):
      for sample in record.samples:
        for gq_field in GQ_FIELDS:
          if gq_field in record.samples[sample] and record.samples[sample][gq_field] is not None:
            record.samples[sample][gq_field] = floor(record.samples[sample][gq_field] / 10)


    with pysam.VariantFile("~{vcf}", 'r') as f_in, pysam.VariantFile("~{outfile}", 'w', header=f_in.header) as f_out:
      for record in f_in:
        if "~{fix_ends_}" == "true":
          fix_bad_end(record)
        if "~{rescale_gq_}" == "true":
          rescale_gq(record)
        f_out.write(record)

    CODE
    tabix ~{outfile}

  >>>
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

