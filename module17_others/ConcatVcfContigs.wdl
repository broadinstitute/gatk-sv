version 1.0
import "Structs.wdl"
import "Utils.wdl" as util

workflow ConcatVcfContigs {
    input {
        Array[File] vcf_list
        Array[File] vcf_list_idx
        outfile_prefix = outfile_prefix,
        sv_base_mini_docker = sv_base_mini_docker,
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    call util.ConcatVcfs as concat_vcfs{
        input:
            vcfs = vcf_list,
            vcfs_idx = vcf_list_idx,
            outfile_prefix = outfile_prefix, 
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    output {
        File vcf_subset = concat_vcfs.concat_vcf
        File vcf_subset_index = concat_vcfs.concat_vcf_idx
    }
 
 }


