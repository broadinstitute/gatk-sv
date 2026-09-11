version 1.0

import "VcfToBed.wdl" as bed
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow VcfToBedAcrossContigs {
  input {
    Array[File] vcfs
    Boolean concat_beds = false
    String prefix
    File primary_contigs_list
    String sv_pipeline_docker
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_svtk_vcf2bed
    RuntimeAttr? runtime_attr_concat_bed
  }

  Array[String] contigs = read_lines(primary_contigs_list)

  scatter (i in range(length(vcfs))) {
    call bed.VcfToBed {
      input:
        vcf = vcfs[i],
        prefix = "~{prefix}.~{contigs[i]}",
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_svtk_vcf2bed=runtime_attr_svtk_vcf2bed
    }
  }

  if (concat_beds) {
    call MiniTasks.ConcatBeds {
      input:
        shard_bed_files = VcfToBed.bed,
        prefix=prefix,
        index_output=true,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_concat_bed
    }
  }

  output {
    Array[File] beds = VcfToBed.bed
    Array[File] bed_idxs = VcfToBed.bed_index
    File? concat_bed = ConcatBeds.merged_bed_file
    File? concat_bed_idx = ConcatBeds.merged_bed_idx
  }
}
