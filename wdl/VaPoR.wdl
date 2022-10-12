version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as mini_tasks
import "VaPoRBedPerContig.wdl" as vapor_bed

workflow VaPoR {
    input {
        String prefix
        String bam_or_cram_file
        String bam_or_cram_index
        File old_vapor_output
        File? bed_file
        String? sample_to_extract_from_bed  # provide sample ID to extact if bed_file contains multiple samples
        File ref_fasta
        File ref_fai
        File ref_dict
        Array[String] contigs
        String vapor_docker
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_Vapor 
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
    }

    if (defined(bed_file)) {
        call vapor_bed.VaPoRBed as VaPoR_bed{
            input:
                prefix = prefix,
                old_vapor_output = old_vapor_output,
                bam_or_cram_file = bam_or_cram_file,
                bam_or_cram_index = bam_or_cram_index,
                bed_file = select_first([bed_file]),
                sample_to_extract = sample_to_extract_from_bed,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                contigs = contigs,
                vapor_docker = vapor_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_Vapor    = runtime_attr_Vapor,
                runtime_attr_bcf2vcf = runtime_attr_bcf2vcf,
                runtime_attr_vcf2bed = runtime_attr_vcf2bed,
                runtime_attr_SplitVcf = runtime_attr_SplitVcf,
                runtime_attr_ConcatBeds = runtime_attr_ConcatBeds
        }
    }
    output {
            File? bed_out = VaPoR_bed.bed
            File? bed_plots = VaPoR_bed.plots
    }
}


