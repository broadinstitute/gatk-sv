version 1.0

import "Utils.wdl" as Utils
import "ExtractVcfGqFilterProperties.wdl" as ExtractVcfGqFilterProperties
import "RecalibrateGq.wdl" as RecalibrateGq
import "BenchmarkGqFilter.wdl" as BenchmarkGqFilter

workflow PreprocessVcfForGqRecalibrator {
    input {
        File input_vcf
        File input_vcf_index
        File? annotations_vcf
        File? annotations_vcf_index
        Array[String] annotations_to_transfer = [
            "INFO/STATUS", "INFO/NON_REF_GENOTYPE_CONCORDANCE", "INFO/VAR_PPV",
            "INFO/VAR_SENSITIVITY", "INFO/TRUTH_AF", "FORMAT/CONC_ST"
        ]
        Boolean standardize_vcf = true
        String sv_utils_docker
        String samtools_cloud_docker
        # optional arguments for overriding default tool behaviors
        Array[String] standardize_vcf_args = []
        Array[String] annotate_gq_args = []
    }

    if(standardize_vcf) {
        call Utils.StandardizeVcfForGatk {
            input:
                vcf=input_vcf,
                standardize_vcf_args=standardize_vcf_args,
                sv_utils_docker=sv_utils_docker
        }
    }
    if(defined(annotations_vcf)) {
        call Utils.TransferVcfAnnotations {
            input:
                vcf_to_annotate=select_first([StandardizeVcfForGatk.fixed_vcf, input_vcf]),
                vcf_to_annotate_index=select_first(
                    [StandardizeVcfForGatk.fixed_vcf_index, input_vcf_index]
                ),
                vcf_with_annotations=select_first([annotations_vcf]),
                vcf_with_annotations_index=select_first([annotations_vcf_index]),
                annotations_to_transfer=annotations_to_transfer,
                samtools_cloud_docker=samtools_cloud_docker
        }
    }

    output {
        File preprocessed_vcf = select_first(
            [TransferVcfAnnotations.annotated_vcf, StandardizeVcfForGatk.fixed_vcf, input_vcf]
        )
        File preprocessed_vcf_index = select_first(
            [
                TransferVcfAnnotations.annotated_vcf_index,
                StandardizeVcfForGatk.fixed_vcf_index,
                input_vcf_index
            ]
        )
    }
}