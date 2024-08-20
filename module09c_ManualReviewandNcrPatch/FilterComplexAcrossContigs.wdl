version 1.0

import "FilterComplex.wdl" as fc

workflow FilterComplexAcrossContigs {
	input {
		Array[File] vcfs
		Array[File] vcf_idxs
		File primary_contigs_list

		Array[String] batch_name_list
        Array[File] PE_metrics
        Array[File] PE_metrics_idxes
        Array[File] Depth_DEL_beds
        Array[File] Depth_DUP_beds

        String prefix
        Int n_per_split
        File sample_PE_metrics
        File sample_depth_calls
        File? script_generate_cpx_review_script

        Boolean use_hail = false
        String? gcs_project

        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker
	}

	Array[String] contigs = read_lines(primary_contigs_list)

	scatter (i in range(length(vcfs))) {
		call fc.FilterComplex {
			input:
				vcf = vcfs[i],
				vcf_idx = vcf_idxs[i],
				contig = contigs[i],
				batch_name_list = batch_name_list,
				PE_metrics = PE_metrics,
				PE_metrics_idxes = PE_metrics_idxes,
				Depth_DUP_beds = Depth_DUP_beds,
				Depth_DEL_beds = Depth_DEL_beds,
				prefix = "~{prefix}.~{contigs[i]}",
				n_per_split = n_per_split,
				sample_PE_metrics = sample_PE_metrics,
				sample_depth_calls = sample_depth_calls,
				script_generate_cpx_review_script = script_generate_cpx_review_script,
				use_hail = use_hail,
				gcs_project = gcs_project,
				sv_base_mini_docker = sv_base_mini_docker,
				sv_pipeline_docker = sv_pipeline_docker,
				sv_pipeline_hail_docker = sv_pipeline_hail_docker
		}
	}

	output {
		Array[File] cpx_filtered_vcfs = FilterComplex.revised_output_vcf
		Array[File] cpx_filtered_vcf_idxs = FilterComplex.revised_output_vcf_idx
		Array[File] cpx_evidences = FilterComplex.cpx_evidences
	}
}
