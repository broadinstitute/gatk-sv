version 1.0

import "Structs.wdl"
import "ManualRevise.wdl" as manual_revise
import "TasksBenchmark.wdl" as Tasks


workflow ReviseVcfWithCpxManualResultsOnly {
    input{
        Array[File] vcf_file_list
        Array[File] vcf_index_list
        File SVID_to_Remove
        File MEI_DEL_Rescue
        Map[String, Array[File]] CPX_manual
        File CTX_manual

        String prefix
        Array[String] contig_list

        String sv_base_mini_docker
        String sv_pipeline_docker
        String sv_pipeline_hail_docker

        RuntimeAttr? runtime_attr_override_concat_revise_files
        RuntimeAttr? runtime_attr_override_revise_vcf
        RuntimeAttr? runtime_attr_override_concat_vcfs
    }

    scatter (i in range(length(vcf_file_list))) {

        String contig = contig_list[i]

        call Tasks.CatUncompressedFiles as CatReviseCpxFiles {
            input:
                shards = CPX_manual[contig],
                outfile_name="~{prefix}.manual_revise_files.~{contig}.txt",
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_override_concat_revise_files
        }

        call manual_revise.ReviseVcf as ReviseVcf {
            input:
                vcf_file = vcf_file_list[i],
                vcf_index = vcf_index_list[i],
                SVID_to_Remove = SVID_to_Remove,
                MEI_DEL_Rescue = MEI_DEL_Rescue,
                CPX_manual = CatReviseCpxFiles.outfile,
                CTX_manual = CTX_manual,
                sv_pipeline_base_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_override_revise_vcf
        }

        call Tasks.ConcatVcfs as ConcatVcfs{
            input:
                vcfs = [ReviseVcf.out_manual_revised_vcf, ReviseVcf.out_cpx_ctx_vcf],
                outfile_prefix = "~{prefix}.~{contig}.manually_revised",
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_override_concat_vcfs
        }
    }

    output{
        Array[File] revised_vcf = ConcatVcfs.concat_vcf
        Array[File] revised_vcf_idx = ConcatVcfs.concat_vcf_idx

    }
}


