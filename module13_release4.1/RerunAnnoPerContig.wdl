#this script is developed to fix false postiive inversions in gnomad v3 vcf that are >1Mb in size but have PE_GT=0
#INVs with no passing samples have filter column revised from "PASS" to "UNRESOLVED"

#developing workdir on erisone: /data/talkowski/xuefang/data/gnomad_V3/module08/step9_sm_depyh_only_dup_fix

version 1.0

import "Structs.wdl"
import "FinalTouchAnnotationFixes.wdl" as final_touch
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge
import "AnnotateFunctionalConsequences.wdl" as func
import "PruneAndAddVafs.wdl" as pav
import "AddGdCnvIntoVcf.wdl" as add_gd_cnv_into_vcf

workflow RerunAnnoPerContig {

    input {
        File CNV_to_Fix #1 column file including the two CNVs that were revised from DEL / DUP and have GT format other than ./.
        File anno_id
        File ped_file
        File par_bed
        File? allosomes_list
        File sample_pop_assignments
        File two_subsets_af_header

        File vcf
        File vcf_idx
        String contig
        String prefix
        File two_subsets_af_info

        Int  sv_per_shard
        Boolean fix_mCNV_formats = false

        File protein_coding_gtf
        File? noncoding_bed
        Int? promoter_window
        Int? max_breakend_as_cnv_length
        String? svannotate_additional_args

        Boolean use_hail
        String? gcs_project

        String gatk_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_base_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_shard_vcf
        RuntimeAttr? runtime_attr_compute_AFs
        RuntimeAttr? runtime_attr_combine_vcfs
        RuntimeAttr? runtime_attr_concat_vcfs
        RuntimeAttr? runtime_attr_svannotate
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_revise_cpx_func_anno
        RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
        RuntimeAttr? runtime_attr_fix_header_sharded_cluster
        RuntimeAttr? runtime_attr_concat_sharded_cluster
        RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
        RuntimeAttr? runtime_attr_preconcat_sharded_cluster
    }

    call MiniTasks.ScatterVcf{
        input:
            vcf = vcf,
            prefix = prefix,
            records_per_shard = sv_per_shard,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_scatter_vcf
    }

    scatter (i in range(length(ScatterVcf.shards))) {

        call RemoveAnno{
            input:
                vcf = ScatterVcf.shards[i],
                vcf_idx = ScatterVcf.shards_idx[i],
                anno_id = anno_id,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        if(fix_mCNV_formats){
            call FixMcnvGT{
                input:
                    vcf = RemoveAnno.anno_removed_vcf,
                    vcf_idx = RemoveAnno.anno_removed_vcf_idx,
                    mCNV_SVID = CNV_to_Fix,
                    sv_pipeline_base_docker = sv_pipeline_base_docker
            }
        }

        File annotation_removed_vcf = select_first([FixMcnvGT.mcnv_gt_fixed_vcf, RemoveAnno.anno_removed_vcf])
        File annotation_removed_vcf_idx = select_first([FixMcnvGT.mcnv_gt_fixed_vcf_idx, RemoveAnno.anno_removed_vcf_idx])

        call func.AnnotateFunctionalConsequences {
            input:
                vcf = annotation_removed_vcf,
                vcf_index = annotation_removed_vcf_idx,
                prefix = "~{prefix}.~{i}",

                protein_coding_gtf = protein_coding_gtf,
                noncoding_bed = noncoding_bed,
                promoter_window = promoter_window,
                max_breakend_as_cnv_length = max_breakend_as_cnv_length,
                additional_args = svannotate_additional_args,
                gatk_docker = gatk_docker,
                runtime_attr_svannotate = runtime_attr_svannotate
        }

        call pav.PruneAndAddVafs as PruneAndAddVafs {
            input:
                vcf                    = AnnotateFunctionalConsequences.annotated_vcf,
                vcf_idx                = AnnotateFunctionalConsequences.annotated_vcf_index,
                prefix                 = prefix,
                contig                 = contig,
                ped_file               = ped_file,
                par_bed                = par_bed,
                allosomes_list         = allosomes_list,
                sample_pop_assignments = sample_pop_assignments,

                sv_base_mini_docker     = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_shard_vcf    = runtime_attr_shard_vcf,
                runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
                runtime_attr_combine_vcfs = runtime_attr_combine_vcfs,
                runtime_attr_concat_vcfs  = runtime_attr_concat_vcfs
        }

        call final_touch.ReviseCpxFunctionalAnnotations{
            input:
                vcf     = PruneAndAddVafs.output_vcf,
                vcf_idx = PruneAndAddVafs.output_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_revise_cpx_func_anno
        }
    }


    Array[File] sharded_annotated_vcf = select_first([ReviseCpxFunctionalAnnotations.cpx_anno_revised_vcf]) 
    Array[File] sharded_annotated_vcf_idx = select_first([ReviseCpxFunctionalAnnotations.cpx_anno_revised_vcf_idx]) 


    if (length(sharded_annotated_vcf) == 0) {
        call MiniTasks.GetVcfHeaderWithMembersInfoLine as GetVcfHeader_annotated {
            input:
                vcf_gz=vcf,
                prefix="~{prefix}.annotated",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_vcf_header_with_members_info_line
       }
    }

    if (length(sharded_annotated_vcf) > 0) {
        if (use_hail) {
            call HailMerge.HailMerge as ConcatVcfsHail_annotated {
                input:
                    vcfs=sharded_annotated_vcf,
                    prefix="~{prefix}.annotated",
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
                outfile_prefix="~{prefix}.annotated",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_concat_sharded_cluster
            }
        }
    }

    File concat_vcf = select_first([GetVcfHeader_annotated.out, ConcatVcfs_annotated.concat_vcf, ConcatVcfsHail_annotated.merged_vcf])
    File concat_vcf_idx = select_first([GetVcfHeader_annotated.out_idx, ConcatVcfs_annotated.concat_vcf_idx, ConcatVcfsHail_annotated.merged_vcf_index])

    call add_gd_cnv_into_vcf.RemoveAllRefSVs as remove_ac_equals_0{
        input:
            vcf = concat_vcf,
            vcf_idx = concat_vcf_idx,
            sv_pipeline_base_docker = sv_pipeline_base_docker

    }

    call RemoveNoncodingAnno{
        input:
            vcf = remove_ac_equals_0.ac_examined_vcf,
            vcf_idx = remove_ac_equals_0.ac_examined_vcf_idx,
            sv_pipeline_base_docker = sv_pipeline_base_docker
    }

    call final_touch.ExtractVcfSites as extract_vcf_sites{
        input:
            vcf = RemoveNoncodingAnno.wo_non_coding_anno_vcf,
            sv_pipeline_base_docker = sv_pipeline_base_docker
    }

    call final_touch.Vcf2Bed as vcf_to_bed{
        input:
            vcf = RemoveNoncodingAnno.wo_non_coding_anno_vcf,
            vcf_idx = RemoveNoncodingAnno.wo_non_coding_anno_vcf_idx,
            sv_pipeline_base_docker = sv_pipeline_base_docker
    }

    call final_touch.AddAfOfSubsets as Add_af_of_2_subsets{
        input:
            vcf = RemoveNoncodingAnno.wo_non_coding_anno_vcf,
            vcf_idx = RemoveNoncodingAnno.wo_non_coding_anno_vcf_idx,
            subsets_af_header_lines = two_subsets_af_header,
            subsets_af_info = two_subsets_af_info,
            sv_pipeline_base_docker = sv_pipeline_base_docker
    }

    call final_touch.ExtractVcfSites as extract_vcf_sites_with_2_subset{
        input:
            vcf = Add_af_of_2_subsets.subset_af_annotated_vcf,
            sv_pipeline_base_docker = sv_pipeline_base_docker
    }


    output {
        File output_vcf = RemoveNoncodingAnno.wo_non_coding_anno_vcf
        File output_vcf_idx = RemoveNoncodingAnno.wo_non_coding_anno_vcf_idx

        File output_site = extract_vcf_sites.vcf_sites
        File output_site_idx = extract_vcf_sites.vcf_sites_idx

        File output_bed = vcf_to_bed.zipped_bed

        File with_two_subset_af_vcf = Add_af_of_2_subsets.subset_af_annotated_vcf
        File with_two_subset_af_vcf_idx = Add_af_of_2_subsets.subset_af_annotated_vcf_idx

        File with_two_subset_af_site = extract_vcf_sites_with_2_subset.vcf_sites
        File with_two_subset_af_site_idx = extract_vcf_sites_with_2_subset.vcf_sites_idx

    }
}



#remove existing annotations
task RemoveAnno{
    input{
        File vcf
        File vcf_idx
        File anno_id

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

        python <<CODE

        import os
        import pysam

        anno = []
        fin=open("~{anno_id}")
        for line in fin:
          pin=line.strip().split()
          anno+=pin
        fin.close()

        fin = pysam.VariantFile('~{vcf}')
        header=fin.header
        for i in anno:
           header.info.remove_header(i)

        fo = pysam.VariantFile("~{prefix}.anno_removed.vcf.gz", 'w', header = header)
        for record in fin:
          for i in anno:
            if i in record.info.keys():
              del record.info[i]
          fo.write(record)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.anno_removed.vcf.gz"

    >>>

    output{
        File anno_removed_vcf = "~{prefix}.anno_removed.vcf.gz"
        File anno_removed_vcf_idx = "~{prefix}.anno_removed.vcf.gz.tbi"
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

task FixMcnvGT{
    input{
        File vcf
        File vcf_idx
        File mCNV_SVID

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

        python <<CODE

        import os
        import pysam

        def svid_readin(svid_to_remove):
            out = []
            fin = open(svid_to_remove)
            for line in fin:
                pin=line.strip().split()
                out += pin
            fin.close()
            return out


        mCNV_svid=svid_readin("~{mCNV_SVID}")

        fin = pysam.VariantFile("~{vcf}")
        fo = pysam.VariantFile("~{prefix}.mCNV_gt_fixed.vcf.gz", 'w', header = fin.header)
        for record in fin:
            if record.id in mCNV_svid and record.info['SVTYPE'] == 'CNV': 
                for sample in record.samples.keys():
                    record.samples[sample]['GT']=(None,None)
                record.info['AC'] = (0,)
                record.info['AN'] = 0
            fo.write(record)
        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.mCNV_gt_fixed.vcf.gz"

    >>>

    output{
        File mcnv_gt_fixed_vcf = "~{prefix}.mCNV_gt_fixed.vcf.gz"
        File mcnv_gt_fixed_vcf_idx = "~{prefix}.mCNV_gt_fixed.vcf.gz.tbi"
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

task RemoveNoncodingAnno{

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

        python <<CODE

        import os
        import pysam

        fin = pysam.VariantFile('~{vcf}')
        header=fin.header
        for i in ['PREDICTED_NONCODING_BREAKPOINT', 'PREDICTED_NONCODING_SPAN']:
           header.info.remove_header(i)

        fo = pysam.VariantFile("~{prefix}.wo_non_coding_anno.vcf.gz", 'w', header = header)
        for record in fin:
          for i in ['PREDICTED_NONCODING_BREAKPOINT', 'PREDICTED_NONCODING_SPAN']:
            if i in record.info.keys():
              del record.info[i]
          fo.write(record)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.wo_non_coding_anno.vcf.gz"

    >>>

    output{
        File wo_non_coding_anno_vcf = "~{prefix}.wo_non_coding_anno.vcf.gz"
        File wo_non_coding_anno_vcf_idx = "~{prefix}.wo_non_coding_anno.vcf.gz.tbi"
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















