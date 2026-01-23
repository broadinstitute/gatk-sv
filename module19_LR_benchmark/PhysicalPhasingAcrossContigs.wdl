version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import "PhysicalPhasingPerContig.wdl" as PhysicalPhasingPerContig

workflow PhysicalPhasingAcrossContigs {

    input  {
        File all_chr_bam
        File all_chr_bai
        File reference_fasta
        File reference_fasta_fai
        File small_vcf
        File small_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        File? trgt_vcf
        File? trgt_vcf_idx
        Int hiphase_memory
        String hiphase_extra_args
        String sample_id
        Array[String] region_list

        String sv_base_mini_docker
        String sv_pipeline_base_docker
    }

    scatter (region in region_list){
        call PhysicalPhasingPerContig.PhysicalPhasingPerContig{
            input:
                all_chr_bam = all_chr_bam,
                all_chr_bai = all_chr_bai,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                small_vcf = small_vcf,
                small_vcf_idx = small_vcf_idx,
                sv_vcf = sv_vcf,
                sv_vcf_idx = sv_vcf_idx,
                trgt_vcf = trgt_vcf,
                trgt_vcf_idx = trgt_vcf_idx,
                hiphase_memory = hiphase_memory,
                hiphase_extra_args = hiphase_extra_args,
                sample_id = sample_id,
                region = region,
                sv_pipeline_base_docker = sv_pipeline_base_docker      
            }
        }

   call LongReadGenotypeTasks.ConcatVcfs as concat_hiphase_vcf{
      input:
        vcfs     = PhysicalPhasingPerContig.hiphase_vcf,
        vcfs_idx = PhysicalPhasingPerContig.hiphase_idx,
        outfile_prefix = "~{sample_id}.phased",
        sv_base_mini_docker = sv_base_mini_docker
    }


    output {
        File hiphase_vcf = concat_hiphase_vcf.concat_vcf
        File hiphase_idx = concat_hiphase_vcf.concat_vcf_idx
    }
}







