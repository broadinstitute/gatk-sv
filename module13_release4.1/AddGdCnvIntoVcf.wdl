#this script is developed to fix false postiive inversions in gnomad v3 vcf that are >1Mb in size but have PE_GT=0
#INVs with no passing samples have filter column revised from "PASS" to "UNRESOLVED"

#developing workdir on erisone: /data/talkowski/xuefang/data/gnomad_V3/module08/step9_sm_depyh_only_dup_fix

version 1.0

import "Structs.wdl"
import "FinalTouchAnnotationFixes.wdl" as final_touch

workflow AddGdCnvIntoVcf {

    input {
        Array[String] contigs
        File SVID_to_remove

        Array[File] vcfs
        Array[File] vcf_idxes
        #Array[File] vcfs_with_5_subsets_af
        #Array[File] vcf_idxes_with_5_subsets_af
        Array[File] vcfs_with_2_subsets_af
        Array[File] vcf_idxes_with_2_subsets_af


        File GD_CNV_vcf
        File GD_CNV_vcf_idx

        File GD_5_subsets_af_header_lines
        File GD_2_subsets_af_header_lines

        File GD_5_subset_af
        File GD_2_subset_af

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_revise_SV_info
    }

    call final_touch.AddAfOfSubsets as Add_af_of_2_subsets_GD{
            input:
                vcf = GD_CNV_vcf,
                vcf_idx = GD_CNV_vcf_idx,
                subsets_af_header_lines = GD_2_subsets_af_header_lines,
                subsets_af_info = GD_2_subset_af,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

    call final_touch.AddAfOfSubsets as Add_af_of_5_subsets_GD{
            input:
                vcf = GD_CNV_vcf,
                vcf_idx = GD_CNV_vcf_idx,
                subsets_af_header_lines = GD_5_subsets_af_header_lines,
                subsets_af_info = GD_5_subset_af,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }


    scatter (i in range(length(vcfs))){

        call AddGdCnv as add_gd_cnv{
            input:
                vcf = vcfs[i],
                vcf_idx = vcf_idxes[i],
                contig = contigs[i],
                GD_CNV_vcf = GD_CNV_vcf,
                GD_CNV_vcf_idx = GD_CNV_vcf_idx,
                SVID_to_remove = SVID_to_remove,

                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_revise_SV_info
        }

        call AddGdCnv as add_gd_cnv_with_2_subset_af {
            input:
                vcf = vcfs_with_2_subsets_af[i],
                vcf_idx = vcf_idxes_with_2_subsets_af[i],
                contig = contigs[i],
                GD_CNV_vcf = Add_af_of_2_subsets_GD.subset_af_annotated_vcf,
                GD_CNV_vcf_idx = Add_af_of_2_subsets_GD.subset_af_annotated_vcf_idx,
                SVID_to_remove = SVID_to_remove,

                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_revise_SV_info
        }

        call RemoveAllRefSVs as remove_ac_equals_0{
            input:
                vcf = add_gd_cnv.with_gd_cnv_vcf,
                vcf_idx = add_gd_cnv.with_gd_cnv_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call RemoveAllRefSVs as remove_ac_equals_0_with_2_subset_af{
            input:
                vcf = add_gd_cnv_with_2_subset_af.with_gd_cnv_vcf,
                vcf_idx =  add_gd_cnv_with_2_subset_af.with_gd_cnv_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call final_touch.ExtractVcfSites as extract_vcf_sites{
            input:
                vcf = remove_ac_equals_0.ac_examined_vcf,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call final_touch.ExtractVcfSites as extract_vcf_sites_with_2_subsets_af{
            input:
                vcf = remove_ac_equals_0_with_2_subset_af.ac_examined_vcf,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call final_touch.Vcf2Bed as vcf_to_bed{
            input:
                vcf = remove_ac_equals_0.ac_examined_vcf,
                vcf_idx = remove_ac_equals_0.ac_examined_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

    }


    output{
        Array[File] vcf_list          = add_gd_cnv.with_gd_cnv_vcf
        Array[File] vcf_idx_list     = add_gd_cnv.with_gd_cnv_vcf_idx
        Array[File] vcf_with_2_subsets_af_list     = add_gd_cnv_with_2_subset_af.with_gd_cnv_vcf
        Array[File] vcf_with_2_subsets_af_idx_list = add_gd_cnv_with_2_subset_af.with_gd_cnv_vcf_idx

        Array[File] sites_list         = extract_vcf_sites.vcf_sites
        Array[File] sites_idx_list     = extract_vcf_sites.vcf_sites_idx
        Array[File] sites_with_2_subsets_af_list     = extract_vcf_sites_with_2_subsets_af.vcf_sites
        Array[File] sites_with_2_subsets_af_idx_list = extract_vcf_sites_with_2_subsets_af.vcf_sites_idx

        Array[File] bed_list         = vcf_to_bed.zipped_bed

    }
}

task AddGdCnv{
    input{
        File vcf
        File vcf_idx
        String contig
        File SVID_to_remove
        File GD_CNV_vcf
        File GD_CNV_vcf_idx

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(100.0 + size(vcf, "GiB") * 6),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail

        python <<CODE

        import os
        import pysam

        def svid_to_remove_readin(svid_to_remove):
            out = []
            fin = open(svid_to_remove)
            for line in fin:
                pin=line.strip().split()
                out += pin
            fin.close()
            return out


        svid_to_remove=svid_to_remove_readin("~{SVID_to_remove}")

        fin = pysam.VariantFile("~{vcf}")
        fo = pysam.VariantFile("~{contig}.temp.vcf.gz", 'w', header = fin.header)
        for record in fin:
            if record.id in svid_to_remove: 
                continue
            else:
                fo.write(record)
        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{contig}.temp.vcf.gz"
        
        bcftools view ~{GD_CNV_vcf} ~{contig} -O z -o ~{contig}.GD_CNV.vcf.gz
        tabix -p vcf "~{contig}.GD_CNV.vcf.gz"

        echo "~{contig}.temp.vcf.gz" > vcf_list.tsv
        echo "~{contig}.GD_CNV.vcf.gz" >> vcf_list.tsv
        bcftools concat --allow-overlaps --file-list vcf_list.tsv -O z -o ~{contig}.with_GD_cnv.vcf.gz
        tabix -p vcf "~{contig}.with_GD_cnv.vcf.gz"

    >>>

    output{
        File with_gd_cnv_vcf = "~{contig}.with_GD_cnv.vcf.gz"
        File with_gd_cnv_vcf_idx = "~{contig}.with_GD_cnv.vcf.gz.tbi"
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


task RemoveAllRefSVs{
    input{
        File vcf
        File vcf_idx

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(100.0 + size(vcf, "GiB") * 6),
        boot_disk_gb: 30,
        preemptible_tries: 1,
        max_retries: 1
    }

    String prefix = basename(vcf,".vcf.gz")
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command<<<
        set -euo pipefail
        bcftools view -i 'SVTYPE=="CNV" || AC>0' ~{vcf} -O z -o "~{prefix}.AC_examined.vcf.gz"
        tabix -p vcf "~{prefix}.AC_examined.vcf.gz"

    >>>

    output{
        File ac_examined_vcf = "~{prefix}.AC_examined.vcf.gz"
        File ac_examined_vcf_idx = "~{prefix}.AC_examined.vcf.gz.tbi"
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


