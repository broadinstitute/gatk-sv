version 1.0

import "Structs.wdl"
import "CalculateFstGnomad.wdl" as calculate_fst
import "CalculateFstLargeVcf.wdl" as calculate_fst_large_vcf

workflow CalculateFstWholeGenomeGnomad {
    input{
        Array[File] vcf_list
        Array[File] vcf_idx_list
        File samp_pop
        String prefix
        String variant_type 
        String sv_fst_docker
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_fst
    }

    scatter(i in range(length(vcf_list))){
        call calculate_fst.CalculateFst{
            input:
                vcf = vcf_list[i],
                vcf_idx = vcf_idx_list[i],
                samp_pop = samp_pop,
                variant_type = variant_type,
                sv_fst_docker = sv_fst_docker
        }
    }

    call calculate_fst_large_vcf.ConcatBeds{
            input:
                shard_bed_files = CalculateFst.output_fst_sites,
                prefix = prefix,
                sv_base_mini_docker = sv_base_mini_docker
        }

    output {
        File output_fst_sites = ConcatBeds.merged_file
    }
}

