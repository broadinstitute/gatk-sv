version 1.0

import "Structs.wdl"
import "LongReadGenotypeTasks.wdl" as LongReadGenotypeTasks
import PhysicalPhasingPerContig.wdl as PhysicalPhasingPerContig

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
    }

    scatter (region in region_list){
        call PhysicalPhasingPerContig.PhysicalPhasingPerContig{
            input:
                all_chr_bam = all_chr_bam,
                all_chr_bai = all_chr_bai,
                reference_fasta = reference_fasta,
                reference_fasta_fai = reference_fasta_fai,
                small_vcf   = small_vcf,
                small_vcf_idx   = small_vcf_idx,
                sv_vcf  = sv_vcf,
                sv_vcf_idx  = sv_vcf_idx,
                trgt_vcf    = trgt_vcf,
                trgt_vcf_idx    = trgt_vcf_idx,
                hiphase_memory  = hiphase_memory,
                hiphase_extra_args  = hiphase_extra_args,
                sample_id   = sample_id,
                region  = region        
            }
        }

   call LongReadGenotypeTasks.ConcatVcfs as concat_small_vcf{
      input:
        vcfs     = PhysicalPhasingPerContig.hiphase_short_vcf,
        vcfs_idx = PhysicalPhasingPerContig.hiphase_short_idx,
        outfile_prefix = "~{sample_id}.small_variants",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LongReadGenotypeTasks.ConcatVcfs as concat_sv_vcf{
      input:
        vcfs     = PhysicalPhasingPerContig.hiphase_sv_vcf,
        vcfs_idx = PhysicalPhasingPerContig.hiphase_sv_idx,
        outfile_prefix = "~{sample_id}.SVs",
        sv_base_mini_docker = sv_base_mini_docker
    }   

    if (defined(trgt_vcf)){
        call LongReadGenotypeTasks.ConcatVcfs as concat_trgt_vcf{
            input:
                vcfs     = select_first([PhysicalPhasingPerContig.hiphase_trgt_vcf]),
                vcfs_idx = select_first([PhysicalPhasingPerContig.hiphase_trgt_idx]),
                outfile_prefix = "~{sample_id}.TRGT",
                sv_base_mini_docker = sv_base_mini_docker
        }
    }

    call LongReadGenotypeTasks.ConcatVcfs as concat_all_vcf{
        input:
            vcfs = [concat_small_vcf.concat_vcf, concat_sv_vcf.concat_vcf, select_first([concat_trgt_vcf.concat_vcf])],
            vcfs_idx = [concat_small_vcf.concat_vcf_idx, concat_sv_vcf.concat_vcf_idx, select_first([concat_trgt_vcf.concat_vcf_idx])],
            outfile_prefix = "~{sample_id}.phased",
            sv_base_mini_docker = sv_base_mini_docker
    }

    output {
        File hiphase_vcf = concat_all_vcf.concat_vcf
        File hiphase_idx = concat_all_vcf.concat_vcf_idx
    }
}







