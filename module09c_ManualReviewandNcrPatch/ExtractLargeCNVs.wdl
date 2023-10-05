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
workflow ExtractLargeCNVs{
    input{ 
        Array[File] vcf_files
        Array[File] vcf_indexes
        String mid_fix
        Int min_size
        String sv_pipeline_base_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_extract_lg_cnv
    }

    scatter(i in range(length(vcf_files))){
        call ExtractLargeCNVsFromVcf{
            input:
                vcf = vcf_files[i],
                vcf_idx = vcf_indexes[i],
                mid_fix = mid_fix,
                min_size = min_size,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_extract_lg_cnv
        }

        call Vcf2Bed{
            input:
                vcf = ExtractLargeCNVsFromVcf.output_vcf,
                sv_pipeline_docker = sv_pipeline_docker
        }
    }

    output{
        Array[File] output_vcf_list = ExtractLargeCNVsFromVcf.output_vcf
        Array[File] output_idx_list = ExtractLargeCNVsFromVcf.output_vcf_idx
        Array[File] output_bed_list = Vcf2Bed.output_bed
    }
 }
   

task ExtractLargeCNVsFromVcf{
    input{
        File vcf
        File vcf_idx
        String mid_fix
        Int min_size
        String sv_pipeline_base_docker
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

    String prefix  = basename(vcf, ".vcf.gz")
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail

        python <<CODE

        import os
        import pysam
        
        fin=pysam.VariantFile("~{vcf}")
        fo=pysam.VariantFile("~{prefix}.~{mid_fix}.vcf.gz",'w', header = fin.header)
        for record in fin:
            if record.info['SVLEN']>~{min_size}:
                if not record.info['SVTYPE'] in ['CPX','INV']:
                    fo.write(record)
        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.~{mid_fix}.vcf.gz"

    >>>

    output{
        File output_vcf = "~{prefix}.~{mid_fix}.vcf.gz"
        File output_vcf_idx = "~{prefix}.~{mid_fix}.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_base_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


task Vcf2Bed{
    input{
        File vcf
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

    String prefix  = basename(vcf, ".vcf.gz")
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail

        svtk vcf2bed -i SVTYPE -i ALL ~{vcf} ~{prefix}.bed
        bgzip ~{prefix}.bed

    >>>

    output{
        File output_bed = "~{prefix}.bed.gz"
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




