version 1.0

import "Utils.wdl" as utils
import "PatchSRBothsidePass.wdl" as patch
import "Structs.wdl"

workflow PatchSRBothsidePassScatter {
    input {
        Array[File] batch_vcfs
        Array[File] cohort_contig_vcfs
        Array[File] updated_bothside_pass_lists
        String cohort_name
        File contig_list

        File patch_script

        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_get_sample_ids
        RuntimeAttr? runtime_attr_get_non_ref_vids
        RuntimeAttr? runtime_attr_calculate_support_frac
    }

    scatter (i in range(length(batch_vcfs))) {
        call utils.GetSampleIdsFromVcf {
            input:
                vcf=batch_vcfs[i],
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_get_sample_ids
        }
    }

    Array[String] contigs = transpose(read_tsv(contig_list))[0]
    scatter ( i in range(length(contigs)) ) {
        call patch.PatchSRBothsidePass {
            input:
                batch_sample_lists=GetSampleIdsFromVcf.out_file,
                cohort_vcf=cohort_contig_vcfs[i],
                updated_bothside_pass_list=updated_bothside_pass_lists[i],
                cohort_name=cohort_name,
                contig=contigs[i],
                patch_script=patch_script,
                sv_base_mini_docker=sv_base_mini_docker,
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_get_non_ref_vids=runtime_attr_get_non_ref_vids,
                runtime_attr_calculate_support_frac=runtime_attr_calculate_support_frac
        }
    }

    output {
        Array[File] out = PatchSRBothsidePass.out
    }
}
