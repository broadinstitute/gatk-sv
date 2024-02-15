#########################################################################################

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
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "SplitCtxForReformat.wdl" as split_ctx_for_reformat
import "RemoveDuplicateEvents.wdl" as remove_duplicate_events
import "ReviseVcfWithManualResults.wdl" as revise_vcf_with_manual_results

### this script is designed to split CTX event from vcf, so they can be manually reformatted and inserted back
workflow SplitCtxForReformatVcfList{
    input{ 
        Array[File] vcf_list
        Array[File] vcf_idx_list
        Array[File] contig_list
        String prefix
        File? CTX_manual

        File? raw_SVs

        Boolean reformat_ctx = true
        Boolean exclude_ctx = false

        String sv_pipeline_docker
        String sv_benchmark_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_extract_ctx_calls 
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_add_raw_SVs
        RuntimeAttr? runtime_attr_exclude_ctx_calls
        RuntimeAttr? runtime_attr_override_reformat_ctx
        RuntimeAttr? runtime_attr_split_raw_SVs_per_chr
        RuntimeAttr? runtime_attr_remove_duplicate_events_task
    }

    scatter(i in range(length(vcf_list))){
        call ExtractCtxCallsV2{
            input:  
                vcf = vcf_list[i],
                vcf_idx = vcf_idx_list[i],
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_extract_ctx_calls
        }

        if (defined(CTX_manual)){
            call ReviseCtxVcf{
                input:
                    vcf_file = ExtractCtxCallsV2.ctx_vcf,
                    vcf_index = ExtractCtxCallsV2.ctx_vcf_idx,
                    CTX_manual = CTX_manual,
                    sv_benchmark_docker = sv_benchmark_docker
            }
        }

        File ctx_vcf_Step1 = select_first([ReviseCtxVcf.out_manual_revised_vcf, ExtractCtxCallsV2.ctx_vcf])
        File ctx_vcf_idx_Step1 = select_first([ReviseCtxVcf.out_manual_revised_vcf_idx, ExtractCtxCallsV2.ctx_vcf_idx])

        if (reformat_ctx){
            call ReformatCtx{
                input:
                    vcf = ctx_vcf_Step1,
                    vcf_idx = ctx_vcf_idx_Step1,
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_override_reformat_ctx
            }
        }


        File ctx_vcf_Step2 = select_first([ReformatCtx.reformatted_vcf, ctx_vcf_Step1])
        File ctx_vcf_idx_Step2 = select_first([ReformatCtx.reformatted_vcf_idx, ctx_vcf_idx_Step1])

        if (defined(raw_SVs)){
            call SplitRawSVsPerChr{
                input:
                    raw_SVs = raw_SVs,
                    contig = contig_list[i],
                    sv_base_mini_docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_attr_split_raw_SVs_per_chr
            }

            call revise_vcf_with_manual_results.AddRawSVs{
                input:
                    prefix = "~{prefix}.~{contig_list[i]}",
                    batch_name = contig_list[i],
                    vcf_file = ctx_vcf_Step2,
                    raw_SVs = SplitRawSVsPerChr.raw_SV_per_chr,
                    sv_benchmark_docker = sv_benchmark_docker,
                    runtime_attr_override = runtime_attr_add_raw_SVs
            }
        }

        File ctx_vcf_Step3 = select_first([AddRawSVs.vcf_with_raw_SVs, ctx_vcf_Step2])
        File ctx_vcf_idx_Step3 = select_first([AddRawSVs.vcf_idx_with_raw_SVs, ctx_vcf_idx_Step2])
     }

    call MiniTasks.ConcatVcfs as ConcatVcfs_CTX {
        input:
            vcfs = ctx_vcf_Step3,
            vcfs_idx = ctx_vcf_idx_Step3,
            allow_overlaps = true,
            outfile_prefix = "~{prefix}.ctx",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_sharded_cluster
    }

    call remove_duplicate_events.RemoveDuplicateEventsTask as RemoveDuplicateEvents{
        input:
            vcf = ConcatVcfs_CTX.concat_vcf,
            vcf_index = ConcatVcfs_CTX.concat_vcf_idx,
            prefix = "~{prefix}.ctx.manual_review",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_remove_duplicate_events_task
    }


    if(exclude_ctx){
        scatter(i in range(length(vcf_list))){
            call ExcludeCtxCallsV2{
                input:  
                    vcf = vcf_list[i],
                    vcf_idx = vcf_idx_list[i],
                    sv_pipeline_docker = sv_pipeline_docker,
                    runtime_attr_override = runtime_attr_exclude_ctx_calls
            }
        }
     }

   output{
        File concat_cpx_vcf = RemoveDuplicateEvents.deduplicated_vcf
        File concat_cpx_vcf_idx = RemoveDuplicateEvents.deduplicated_vcf_index
        Array[File] cpx_vcf_list = ctx_vcf_Step2
        Array[File] ctx_vcf_idx_list = ctx_vcf_idx_Step2
        Array[File]? no_ctx_vcf_list = ExcludeCtxCallsV2.no_ctx_vcf
        Array[File]? no_ctx_vcf_idx_list = ExcludeCtxCallsV2.no_ctx_vcf_idx
    }
}

task ExtractCtxCallsV2{

    input{
        File vcf
        File vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 200,
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix  = basename(vcf, ".vcf.gz")

    command<<<
        set -euo pipefail

        bcftools view -Oz --include 'INFO/SVTYPE="CTX"' -o ~{prefix}.ctx.vcf.gz ~{vcf}

        tabix -p vcf "~{prefix}.ctx.vcf.gz"

    >>>

    output{
        File ctx_vcf = "~{prefix}.ctx.vcf.gz"
        File ctx_vcf_idx ="~{prefix}.ctx.vcf.gz.tbi"
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


task ExcludeCtxCallsV2{
    input{
        File vcf
        File vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
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
        fin2=pysam.VariantFile("~{vcf}")
        fo2=pysam.VariantFile("~{prefix}.no_ctx.vcf.gz", 'w', header = fin2.header)
        for record in fin2:
            if record.info['SVTYPE']!="CTX":
                fo2.write(record)
        fin2.close()
        fo2.close()
        CODE

        tabix -p vcf "~{prefix}.no_ctx.vcf.gz"

    >>>

    output{
        File no_ctx_vcf = "~{prefix}.no_ctx.vcf.gz"
        File no_ctx_vcf_idx ="~{prefix}.no_ctx.vcf.gz.tbi"
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

task ReformatCtx{
    
    input{
        File vcf
        File vcf_idx
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 20,
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

        def extract_info_char(info,char):
            out = ''
            for i in info.split(';'):
                if i.split('=')[0] == char:
                    out = i.split('=')[1]
                    break
            return out

        def update_info_char(info, char, value):
            tmp = []
            for i in info.split(';'):
                if i.split('=')[0] == char:
                    i_new =  i.split('=')[0]+'='+str(value)
                else:
                    i_new = i
                tmp.append(i_new)
            return ';'.join(tmp)

        fin = os.popen(r'''zcat %s'''%("~{vcf}"))
        fo = open("~{prefix}.reformatted.vcf", 'w')
        for line in fin:
            pin=line.strip().split('\t')
            if pin[0][0]=='#':
                print('\t'.join(pin), file=fo)
            else:
                svtype = extract_info_char(pin[7],'SVTYPE')
                if svtype=="CTX" and 'END2' in pin[7]:
                    end = int(extract_info_char(pin[7],'END'))
                    end2 = int(extract_info_char(pin[7],'END2'))
                    if end==end2:
                        end_new = int(pin[1])+1
                    elif end==int(pin[1]):
                        end_new = int(pin[1])+1
                    info_new = update_info_char(pin[7], 'END', end_new)
                    pin[7] = info_new
                print('\t'.join(pin), file=fo)
        fin.close()
        fo.close()
        CODE

        bgzip "~{prefix}.reformatted.vcf"
        tabix -p vcf "~{prefix}.reformatted.vcf.gz"

    >>>

    output{
        File reformatted_vcf = "~{prefix}.reformatted.vcf.gz"
        File reformatted_vcf_idx ="~{prefix}.reformatted.vcf.gz.tbi"
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

task ReviseCtxVcf{
    input{
        File vcf_file
        File vcf_index
        File? CTX_manual
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


        python <<CODE

        def CTX_manual_readin(CTX_manual):
            fin=open(CTX_manual)
            out={}
            for line in fin:
                pin=line.strip().split('\t')
                if not pin[0] in out.keys():
                    out[pin[0]] = {}
                if not pin[1] in out[pin[0]].keys():
                    out[pin[0]][pin[1]] = pin[3:]
            fin.close()
            return out

        def revise_vcf(vcf_input, vcf_output, hash_CTX_manual):
            fin=pysam.VariantFile(vcf_input)
            header = fin.header
            header.formats.add('MANUAL', '1', 'String', 'Reason for a failure from manual review')
            fo=pysam.VariantFile(vcf_output, 'w', header = header)
            for record in fin:
                if record.id in hash_CTX_manual.keys():
                    if 'NA' in hash_CTX_manual[record.id].keys():
                        record.filter.add('UNRESOLVED')
                    else:
                        for sample in hash_CTX_manual[record.id].keys():
                            if sample in record.samples.keys():
                                if hash_CTX_manual[record.id][sample][1] in ['Not_Enough_PE_Pairs', 'No_PE_Pairs']:
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
                ref_count = len([s for s in record.samples if record.samples[s]['GT'] in NULL_and_REF_GTs])
                alt_count = len(record.samples) - ref_count
                if alt_count>0 or record.info['SVTYPE']=="CNV":
                    fo.write(record)
            fin.close()
            fo.close()

        import os
        import sys
        from numpy import median
        import pysam
        import argparse

        NULL_GTs = [(None, None), (None, )]
        REF_GTs = [(0, 0), (0, ), (None, 2)]
        NULL_and_REF_GTs = NULL_GTs + REF_GTs
        HET_GTs = [(0, 1), (None, 1), (None, 3)]

        #Define global variables
        hash_CTX_manual =  CTX_manual_readin("~{CTX_manual}")
        revise_vcf("~{vcf_file}", "~{prefix}.Manual_Revised.vcf.gz", hash_CTX_manual)

        CODE

        tabix -p vcf "~{prefix}.Manual_Revised.vcf.gz"
    >>>

    output{
        File out_manual_revised_vcf = "~{prefix}.Manual_Revised.vcf.gz"
        File out_manual_revised_vcf_idx = "~{prefix}.Manual_Revised.vcf.gz.tbi"
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


