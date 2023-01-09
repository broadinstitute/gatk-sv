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
import "ManualRevise.wdl" as manual_revise

workflow ManualReviseCTX{
    input{ 
        File vcf
        File vcf_idx
        File? CTX_manual

        String sv_base_mini_docker
     }

     call ReviseCtxVcf{
        input:
            vcf_file = vcf,
            vcf_index = vcf_idx,
            CTX_manual = CTX_manual,
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_revise_ctx_vcf
     }


    output{
        File revised_vcf = ReviseCtxVcf.out_manual_revised_vcf
        File revised_vcf_idx = ReviseCtxVcf.out_manual_revised_vcf_idx
    }
}


task ReviseCtxVcf{
    input{
        File vcf_file
        File vcf_index
        File? CTX_manual
        String sv_base_mini_docker
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
            fo=pysam.VariantFile(vcf_output, 'w', header = fin.header)
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
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


