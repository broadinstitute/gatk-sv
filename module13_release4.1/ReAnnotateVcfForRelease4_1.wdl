version 1.0

import "Structs.wdl"
import "ShardedAnnotateVcf.wdl" as sharded_annotate_vcf
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge
import "FinalTouchAnnotationFixes.wdl" as final_touch


workflow ReAnnotateVcfForRelease4_1{
    input{
        Array[File] vcf_list
        Array[File] vcf_idx_list
        Array[String] contig_list
        Array[String] prefix_list
        File? protein_coding_gtf
        File? noncoding_bed
        Int? promoter_window
        Int? max_breakend_as_cnv_length
        String? svannotate_additional_args

        File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
        File? sample_keep_list              # List of samples to be retained from the output vcf
        File? ped_file                # Used for M/F AF calculations
        File? par_bed
        File? allosomes_list
        Int   sv_per_shard

        File? ref_bed              # File with external allele frequencies
        String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
        Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

        Boolean use_hail
        String? gcs_project

        String sv_pipeline_docker
        String sv_pipeline_hail_docker
        String sv_base_mini_docker
        String gatk_docker
        String sv_pipeline_base_docker

        RuntimeAttr? runtime_attr_clean_annotations
        RuntimeAttr? runtime_attr_svannotate
        RuntimeAttr? runtime_attr_scatter_vcf
        RuntimeAttr? runtime_attr_subset_vcf_by_samples_list
        RuntimeAttr? runtime_attr_compute_AFs
        RuntimeAttr? runtime_attr_modify_vcf
        RuntimeAttr? runtime_attr_split_ref_bed
        RuntimeAttr? runtime_attr_split_query_vcf
        RuntimeAttr? runtime_attr_bedtools_closest
        RuntimeAttr? runtime_attr_select_matched_svs
        RuntimeAttr? runtime_attr_concat
        RuntimeAttr? runtime_attr_preconcat
        RuntimeAttr? runtime_attr_hail_merge
        RuntimeAttr? runtime_attr_fix_header
    }

    scatter (i in range(length(contig_list))){
        call CleanAnnotation{
            input:
                vcf = vcf_list[i],
                vcf_idx = vcf_idx_list[i],
                sv_pipeline_base_docker = sv_pipeline_base_docker,
                runtime_attr_override = runtime_attr_clean_annotations
        }

        call sharded_annotate_vcf.ShardedAnnotateVcf {
            input:
                vcf = CleanAnnotation.anno_removed_vcf,
                vcf_idx = CleanAnnotation.anno_removed_vcf_idx,
                contig = contig_list[i],
                prefix = prefix_list[i],
                protein_coding_gtf = protein_coding_gtf,
                noncoding_bed = noncoding_bed,
                promoter_window = promoter_window,
                svannotate_additional_args = svannotate_additional_args,
                max_breakend_as_cnv_length = max_breakend_as_cnv_length,

                sample_pop_assignments = sample_pop_assignments,
                ped_file = ped_file,
                par_bed = par_bed,
                sv_per_shard = sv_per_shard,
                allosomes_list = allosomes_list,

                min_records_per_shard_step1 = 200,
                max_shards_per_chrom_step1 = 5000,


                ref_bed = ref_bed,
                ref_prefix = ref_prefix,
                population = population,

                use_hail = use_hail,
                gcs_project = gcs_project,

                gatk_docker = gatk_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_hail_docker = sv_pipeline_hail_docker,
                sv_pipeline_base_docker = sv_pipeline_base_docker,

                runtime_attr_svannotate = runtime_attr_svannotate,
                runtime_attr_scatter_vcf = runtime_attr_scatter_vcf,
                runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
                runtime_attr_modify_vcf = runtime_attr_modify_vcf,
                runtime_attr_split_ref_bed  = runtime_attr_split_ref_bed,
                runtime_attr_split_query_vcf  = runtime_attr_split_query_vcf,
                runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
                runtime_attr_select_matched_svs = runtime_attr_select_matched_svs
        }

        call AddEvidenceForGD{
            input:
                vcf = ShardedAnnotateVcf.output_vcf,
                vcf_idx = ShardedAnnotateVcf.output_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call final_touch.ExtractVcfSites as extract_vcf_sites{
            input:
                vcf = AddEvidenceForGD.GD_Evi_added_vcf,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }

        call final_touch.Vcf2Bed as vcf_to_bed{
            input:
                vcf = AddEvidenceForGD.GD_Evi_added_vcf,
                vcf_idx = AddEvidenceForGD.GD_Evi_added_vcf_idx,
                sv_pipeline_base_docker = sv_pipeline_base_docker
        }
   }

    output {
        Array[File] annotated_vcfs = ShardedAnnotateVcf.output_vcf
        Array[File] annotated_vcf_idxes = ShardedAnnotateVcf.output_vcf_idx
        Array[File] annotated_sites = extract_vcf_sites.vcf_sites
        Array[File] annotated_site_idxes = extract_vcf_sites.vcf_sites_idx
        Array[File] annotated_beds = vcf_to_bed.zipped_bed
    }

}


task CleanAnnotation{
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

        fin=pysam.VariantFile("~{vcf}")
        header = fin.header
        key_to_remove = []
        for key in header.info.keys():
          if 'PREDICTED_' in key:
            header.info.remove_header(key)
            key_to_remove.append(key)
          if 'AF' in key or 'AC' in key or 'AN' in key:
            header.info.remove_header(key)
            key_to_remove.append(key)
          if 'N_BI_GENOS' in key:
            header.info.remove_header(key)
            key_to_remove.append(key)
          if 'N_HOMREF' in key or 'N_HET' in key or 'N_HOMALT' in key:
            header.info.remove_header(key)
            key_to_remove.append(key)
          if 'FREQ_HOMREF' in key or 'FREQ_HET' in key or 'FREQ_HOMALT' in key:
            header.info.remove_header(key)
            key_to_remove.append(key)
          if 'CN_NUMBER' in key or 'CN_COUNT' in key or 'CN_FREQ' in key or 'CN_NONREF_COUNT' in key or 'CN_NONREF_FREQ' in key:
            header.info.remove_header(key)
            key_to_remove.append(key)
          if 'HEMIREF' in key or 'HEMIALT' in key:
            header.info.remove_header(key)
            key_to_remove.append(key)

        fo = pysam.VariantFile("~{prefix}.AF_Func_annotation_cleaned.vcf.gz", 'w', header =  header)
        for rec in fin:
          for key in rec.info.keys():
            if key in key_to_remove:
              del rec.info[key]
          fo.write(rec)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.AF_Func_annotation_cleaned.vcf.gz"

    >>>

    output{
        File anno_removed_vcf = "~{prefix}.AF_Func_annotation_cleaned.vcf.gz"
        File anno_removed_vcf_idx = "~{prefix}.AF_Func_annotation_cleaned.vcf.gz.tbi"
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

task AddEvidenceForGD{
    input{
        File vcf
        File vcf_idx

        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(100.0 + size(vcf, "GiB") * 3),
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

        fin=pysam.VariantFile("~{vcf}")
        header = fin.header
        fo = pysam.VariantFile("~{prefix}.GD_Evi_fixed.vcf.gz", 'w', header =  header)
        for rec in fin:
            if not "EVIDENCE" in rec.info.keys():
                print(rec.id)
                rec.info['EVIDENCE'] = ('RD')
            fo.write(rec)

        fin.close()
        fo.close()

        CODE

        tabix -p vcf "~{prefix}.GD_Evi_fixed.vcf.gz"

    >>>

    output{
        File GD_Evi_added_vcf = "~{prefix}.GD_Evi_fixed.vcf.gz"
        File GD_Evi_added_vcf_idx = "~{prefix}.GD_Evi_fixed.vcf.gz.tbi"
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


