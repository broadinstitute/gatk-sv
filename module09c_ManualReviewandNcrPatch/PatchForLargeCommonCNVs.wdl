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
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow PatchForLargeCommonCNVs{
    #this workflow is developed to fix the bug where false positive large common CNVs reside in the vcf and overlapping smaller CNVs got down-genotyped
    input{ 
        File before_cleanvcf_vcf   #the vcf file after complex regenotype and before cleanvc
        File before_cleanvcf_vcf_idx
        File after_filter_vcf      #the vcf after cleanvcf, machine learning filterings, outlier removal, and recluster of small SVs
        File after_filter_vcf_idx
        String svtype
        String mid_fix
        Int min_size
        String contig   #chromosome of the SV
        Int min_pos     #start coordinate of the SV
        Int max_pos     #end coordinate of the SV
        File ped_file
        String prefix
        String vcf_prefix

        String sv_pipeline_docker
        String sv_pipeline_base_docker
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_extract_lg_cnv
        RuntimeAttr? runtime_attr_extract_overlap_SVs
        RuntimeAttr? runtime_attr_svtk_cluster_lg_cnvs
        RuntimeAttr? runtime_attr_strip_out_overlap_SVs
        RuntimeAttr? runtime_attr_concat_vcf
        RuntimeAttr? runtime_attr_remove_outlier_samples
        RuntimeAttr? runtime_attr_revise_evidence_tag
    }

    
    call ExtractLargeCNVsFromVcf{
        input:
            vcf = before_cleanvcf_vcf, 
            vcf_idx = before_cleanvcf_vcf_idx,
            mid_fix = mid_fix,
            min_size = min_size,
            sv_pipeline_base_docker = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_extract_lg_cnv
        }

    call ExtractOverlappingSVsChrXDel{
        input:
            vcf = ExtractLargeCNVsFromVcf.output_vcf,
            vcf_idx = ExtractLargeCNVsFromVcf.output_vcf_idx,
            ped_file = ped_file,
            contig = contig, 
            svtype = svtype,
            min_size = min_size,
            min_pos = min_pos,
            max_pos = max_pos,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_extract_overlap_SVs
       }

    call SvtkClusterLgCNVs{
        input:
            vcf = ExtractOverlappingSVsChrXDel.output_vcf,
            vcf_idx = ExtractOverlappingSVsChrXDel.output_idx,
            dist = 1000000000,
            frac = 0.5,
            svsize = 0,
            sv_types = [svtype],
            sample_overlap = 0.5,
            vid_prefix = "~{prefix}.~{contig}.final_cleanup_~{svtype}_patched_reclustered_",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_svtk_cluster_lg_cnvs
        }

    call StripOutOverlapingSVs{
        input:
            vcf = after_filter_vcf,
            vcf_idx = after_filter_vcf_idx,
            contig = contig,
            svtype = svtype,
            min_size = min_size,
            min_pos = min_pos,
            max_pos = max_pos,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_strip_out_overlap_SVs
        }

    call RemoveOutlierSamples{
        input:
            vcf = SvtkClusterLgCNVs.vcf_out,
            vcf_idx = SvtkClusterLgCNVs.vcf_out_idx,
            vcf_outlier_removed = after_filter_vcf,
            vcf_idx_outlier_removed = after_filter_vcf_idx,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_remove_outlier_samples
    }

    call ReviseEvidenceTags{
        input:
            vcf = RemoveOutlierSamples.output_vcf,
            vcf_idx = RemoveOutlierSamples.output_vcf_idx,
            sv_pipeline_base_docker = sv_pipeline_base_docker,
            runtime_attr_override = runtime_attr_revise_evidence_tag
    }
    
    call MiniTasks.ConcatVcfs{
        input:
            vcfs = [ReviseEvidenceTags.output_vcf, StripOutOverlapingSVs.output_vcf],
            vcfs_idx = [ReviseEvidenceTags.output_vcf_idx , StripOutOverlapingSVs.output_vcf_idx],
            allow_overlaps = true,
            outfile_prefix = "~{vcf_prefix}.patch",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcf
        }

    output{
        File patched_vcf = ConcatVcfs.concat_vcf
        File patched_vcf_idx = ConcatVcfs.concat_vcf_idx

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

task ExtractOverlappingSVsChrXDel{
    input{
        File vcf
        File vcf_idx
        File ped_file
        String contig
        String svtype
        Int min_pos
        Int max_pos
        Int min_size
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    String prefix  = basename(vcf, ".vcf.gz")
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail

        bcftools view ~{vcf} ~{contig}:~{min_pos}-~{max_pos} | bgzip > "~{prefix}.target.vcf.gz"
        tabix -p vcf "~{prefix}.target.vcf.gz"
        bcftools view -h ~{prefix}.target.vcf.gz | awk '{if ($1=="#CHROM") print}' | cut -f10- | sed -e 's/\t/\n/g' > sample_list.tsv

        python <<CODE

        def readin_male_samples(ped, sample_list):
            out = []
            fin=open(ped)
            for line in fin:
                pin=line.strip().split()
                if pin[4]=='1':
                    if pin[1] in sample_list:
                        out.append(pin[1])
            fin.close()
            return  out

        def readin_female_samples(ped,sample_list):
            out = []
            fin=open(ped)
            for line in fin:
                pin=line.strip().split()
                if pin[4]=='2':
                    if pin[1] in sample_list:
                        out.append(pin[1])
            fin.close()
            return  out

        def readin_sample_names(samp_all):
            fin=open(samp_all)
            out=[]
            for line in fin:
                pin=line.strip().split()
                out+=pin
            fin.close()
            return   out

        import os
        import pysam

        sample_list = readin_sample_names("sample_list.tsv")
        male_list = readin_male_samples("~{ped_file}", sample_list)
        female_list = readin_female_samples("~{ped_file}", sample_list)

        fin=pysam.VariantFile("~{prefix}.target.vcf.gz")
        fo=pysam.VariantFile("~{prefix}.regenotyped.vcf.gz", 'w', header = fin.header)
        for record in fin:
            if record.info["SVTYPE"] == "~{svtype}" and not record.info["SVLEN"]<~{min_size}: 
                for sample in male_list:
                    if record.samples[sample]['RD_CN']==0:
                        record.samples[sample]['GT']=(0, 1)
                    else:
                        record.samples[sample]['GT']=(0, 0)
                for sample in female_list:
                    if record.samples[sample]['RD_CN']==0:
                        record.samples[sample]['GT']=(1, 1)
                    elif record.samples[sample]['RD_CN']==1:
                        record.samples[sample]['GT']=(0, 1)
                    else:
                        record.samples[sample]['GT']=(0, 0)
                fo.write(record)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.regenotyped.vcf.gz"

    >>>

    output{
        File target_vcf = "~{prefix}.target.vcf.gz"
        File target_idx = "~{prefix}.target.vcf.gz.tbi"
        File output_vcf = "~{prefix}.regenotyped.vcf.gz"
        File output_idx = "~{prefix}.regenotyped.vcf.gz.tbi"
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

task SvtkClusterLgCNVs {
    input {
        File vcf
        File vcf_idx
        Int dist
        Float frac
        Float sample_overlap
        File? exclude_list
        File? exclude_list_idx
        Int svsize
        Array[String] sv_types
        String vid_prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float default_mem_gb = 100
    RuntimeAttr runtime_default = object {
        mem_gb: default_mem_gb,
        disk_gb: ceil(50.0 + size(vcf, "GiB") * 50.0),
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
        }

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    String prefix = basename(vcf,".vcf.gz")

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
        svtk vcfcluster <(echo "~{vcf}") - \
                -d ~{dist} \
                -f ~{frac} \
                ~{if defined(exclude_list) then "-x ~{exclude_list}" else ""} \
                -z ~{svsize} \
                -p ~{vid_prefix} \
                -t ~{sep=',' sv_types} \
                -o ~{sample_overlap} \
                --preserve-ids \
                --preserve-genotypes \
                --preserve-header \
                | bcftools sort -O z -o ~{prefix}.vcf.gz - 

        tabix -p vcf ~{prefix}.vcf.gz
    >>>

    output {
        File vcf_out = "~{prefix}.vcf.gz"
        File vcf_out_idx = "~{prefix}.vcf.gz.tbi"
    }
}

task StripOutOverlapingSVs{
    input{
        File vcf
        File vcf_idx
        String svtype
        String contig
        Int min_size
        Int min_pos
        Int max_pos
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

        bcftools view ~{vcf} ~{contig}:~{min_pos}-~{max_pos} | bgzip > "~{prefix}.target.vcf.gz"
        tabix -p vcf "~{prefix}.target.vcf.gz"

        python <<CODE

        import os
        import pysam

        target_SVID = []
        fin=pysam.VariantFile("~{prefix}.target.vcf.gz")
        for record in fin:
            if record.info["SVTYPE"] == "~{svtype}" and not record.info["SVLEN"]<~{min_size}: 
                target_SVID.append(record.id)
        fin.close()

        fin=pysam.VariantFile("~{vcf}")
        fo=pysam.VariantFile("~{prefix}.target_removed.vcf.gz", 'w', header = fin.header)
        for record in fin:
            if not record.id in target_SVID:
                fo.write(record)
        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.target_removed.vcf.gz"

    >>>

    output{
        File output_vcf = "~{prefix}.target_removed.vcf.gz"
        File output_vcf_idx = "~{prefix}.target_removed.vcf.gz.tbi"
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

task RemoveOutlierSamples{
    input{
        File vcf
        File vcf_idx
        File vcf_outlier_removed
        File vcf_idx_outlier_removed
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

        bcftools view -h ~{vcf_outlier_removed} | awk '{if ($1=="#CHROM") print}' | cut -f10- | sed -e 's/\t/\n/g' > sample_list.tsv
        bcftools view -S sample_list.tsv ~{vcf} -O z -o ~{prefix}.outlier_removed.vcf.gz
        tabix -p vcf "~{prefix}.outlier_removed.vcf.gz"
    >>>

    output{
        File output_vcf = "~{prefix}.outlier_removed.vcf.gz"
        File output_vcf_idx = "~{prefix}.outlier_removed.vcf.gz.tbi"
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

task ReviseEvidenceTags{
    input{
        File vcf
        File vcf_idx
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

        evidence_tag_hash = {}
        evidence_tag_hash[('1',)] = ("RD",)
        evidence_tag_hash[('2',)] = ("PE",)
        evidence_tag_hash[('3',)] = ("RD","PE")
        evidence_tag_hash[('4',)] = ("SR",)
        evidence_tag_hash[('5',)] = ("RD","SR")
        evidence_tag_hash[('6',)] = ("PE","SR")
        evidence_tag_hash[('7',)] = ("RD","PE","SR")
        evidence_tag_hash[(None,)] = (".",)

        fin=pysam.VariantFile("~{vcf}")
        header = fin.header
        header.formats.remove_header('EV')
        header.add_line('##FORMAT=<ID=EV,Number=.,Type=String,Description="Classes of evidence supporting final genotype">')

        fo=pysam.VariantFile("~{prefix}.evidence_tag_revised.vcf.gz", 'w', header = header)
        for record in fin:
            for sample in record.samples:
                new_EV = evidence_tag_hash[record.samples[sample]['EV']]
                record.samples[sample]['EV'] = new_EV
            fo.write(record)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.evidence_tag_revised.vcf.gz"
    >>>

    output{
        File output_vcf = "~{prefix}.evidence_tag_revised.vcf.gz"
        File output_vcf_idx = "~{prefix}.evidence_tag_revised.vcf.gz.tbi"
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




