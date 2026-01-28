version 1.0

import "Structs.wdl"
import "Utils.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow SubsetVcfToSamples {
    input {
        File vcf
        File vcf_idx
        File samples_list
        Array[String] contigs
        String prefix
        
        String sv_base_mini_docker

        RuntimeAttr? runtime_attr_subset_vcf
        RuntimeAttr? runtime_attr_subset_samples
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (contig in contigs) {
        call Utils.SubsetVcfToContig {
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                contig = contig,
                prefix = prefix,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_subset_vcf
        }

        call Utils.SubsetVcfBySamplesList {
            input:
                vcf = SubsetVcfToContig.subset_vcf,
                vcf_idx = SubsetVcfToContig.subset_vcf_idx,
                list_of_samples = samples_list,
                outfile_name = "~{prefix}.contig.subset_samples",
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_subset_samples
        }
    }

    call MiniTasks.ConcatVcfs {
        input:
            vcfs = SubsetVcfBySamplesList.vcf_subset,
            vcfs_idx = SubsetVcfBySamplesList.vcf_subset_index,
            allow_overlaps = true,
            outfile_prefix = "~{prefix}.subset_samples",
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    output {
        File subset_vcf = ConcatVcfs.concat_vcf
        File subset_vcf_idx = ConcatVcfs.concat_vcf_idx
    }
}
