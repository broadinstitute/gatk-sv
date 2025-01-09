version 1.0
import "SubsetVcfBySamples.wdl" as subset_vcf_by_samples
import "Structs.wdl"
import "Utils.wdl" as util

workflow SubsetVcfBySamplesWholeGenome {
    input {
        File vcf
        File vcf_index
        File list_of_samples    # List of samples to keep (default, remove_samples = false) or remove (remove_samples = true)
        Array[String] contigs
        String? outfile_prefix
        Boolean? remove_samples    # If false (default), keep samples in provided list. If true, remove them.
        Boolean? remove_private_sites    # If true (default), remove sites that are private to excluded samples. If false, keep sites even if no remaining samples are non-ref.

        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_subset_vcf_by_contig
        RuntimeAttr? runtime_attr_subset_by_samples
        RuntimeAttr? runtime_attr_concat_vcfs
    }

    scatter (contig in contigs){
        call util.SubsetVcfByContig as subset_vcf_by_contig{
            input:
                vcf = vcf,
                vcf_idx = vcf_index,
                contig = contig,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_subset_vcf_by_contig
        }

        call util.SubsetVcfBySamplesList as subset_vcf_by_sample_list {
            input:
                vcf = subset_vcf_by_contig.vcf_subset_by_contig,
                vcf_idx = subset_vcf_by_contig.vcf_subset_by_contig_index,
                list_of_samples = list_of_samples,
                remove_samples = remove_samples,
                remove_private_sites = remove_private_sites,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_subset_by_samples
        }
    }

    call util.ConcatVcfs as concat_vcfs{
        input:
            vcfs = subset_vcf_by_sample_list.vcf_subset_by_sample,
            vcfs_idx = subset_vcf_by_sample_list.vcf_subset_by_sample_index,
            outfile_prefix = outfile_prefix, 
            sv_base_mini_docker = sv_base_mini_docker,
            runtime_attr_override = runtime_attr_concat_vcfs
    }

    output {
        File vcf_subset = concat_vcfs.concat_vcf
        File vcf_subset_index = concat_vcfs.concat_vcf_idx
    }
 
 }


