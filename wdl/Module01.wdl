version 1.0

import "GATKSVTools.wdl" as gatk
import "Tasks0506.wdl" as tasks0506

workflow Module01 {
  input {
    File vcf
    String batch
    File contig_list
    File ref_dict

    File pesr_include_intervals
    Float pesr_min_inclusion_interval_overlap
    Int pesr_min_size

    File depth_include_intervals
    Float depth_min_inclusion_interval_overlap
    Int depth_min_size

    String cluster1_algorithm = "SINGLE_LINKAGE"
    String cluster1_breakpoint_summary_strategy = "MEDIAN_START_MEDIAN_END"
    Float cluster1_depth_overlap_fraction = 0.9
    Float cluster1_mixed_overlap_fraction = 0.9
    Float cluster1_pesr_overlap_fraction = 0.9
    Int cluster1_depth_breakend_window = 50
    Int cluster1_mixed_breakend_window = 50
    Int cluster1_pesr_breakend_window = 50
    Float cluster1_depth_sample_overlap = 0
    Float cluster1_mixed_sample_overlap = 0
    Float cluster1_pesr_sample_overlap = 0

    String cluster2_algorithm = "MAX_CLIQUE"
    String cluster2_breakpoint_summary_strategy = "MEDIAN_START_MEDIAN_END"
    Float cluster2_depth_overlap_fraction = 0.5
    Float cluster2_mixed_overlap_fraction = 0.5
    Float cluster2_pesr_overlap_fraction = 0.5
    Int cluster2_depth_breakend_window = 1000
    Int cluster2_mixed_breakend_window = 500
    Int cluster2_pesr_breakend_window = 500
    Float cluster2_depth_sample_overlap = 0
    Float cluster2_mixed_sample_overlap = 0
    Float cluster2_pesr_sample_overlap = 0

    String gatk_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_select_pesr
    RuntimeAttr? runtime_attr_select_depth
    RuntimeAttr? runtime_attr_annotate_regions_pesr
    RuntimeAttr? runtime_attr_annotate_regions_depth
    RuntimeAttr? runtime_attr_select_regions_pesr
    RuntimeAttr? runtime_attr_select_regions_depth
    RuntimeAttr? runtime_attr_cluster_1
    RuntimeAttr? runtime_attr_cluster_2
    RuntimeAttr? runtime_attr_select_size
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  call gatk.SelectVariants as SelectPESR {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      output_name = "~{batch}.pesr",
      select_expression="ALGORITHMS!=\"depth\"",
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_select_pesr
  }

  call gatk.SelectVariants as SelectDepth {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      output_name = "~{batch}.depth",
      select_expression = "ALGORITHMS==\"depth\"",
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_select_depth
  }

  String pesr_region_name = "INCL_PESR"
  call gatk.SVAnnotateOverlappingRegions as AnnotateIncludedRegionsPESR {
    input:
      vcf = SelectPESR.out,
      vcf_index = SelectPESR.out_index,
      output_name = "~{batch}.pesr.ann_regions",
      region_names = [pesr_region_name],
      region_files = [pesr_include_intervals],
      require_breakend_overlap = true,
      min_overlap_fraction = pesr_min_inclusion_interval_overlap,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_annotate_regions_pesr
  }

  String depth_region_name = "INCL_DEPTH"
  call gatk.SVAnnotateOverlappingRegions as AnnotateIncludedRegionsDepth {
    input:
      vcf = SelectDepth.out,
      vcf_index = SelectDepth.out_index,
      output_name = "~{batch}.depth.ann_regions",
      region_names = [depth_region_name],
      region_files = [depth_include_intervals],
      require_breakend_overlap = false,
      min_overlap_fraction = depth_min_inclusion_interval_overlap,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_annotate_regions_depth
  }

  call gatk.SelectVariants as SelectRegionsPESR {
    input:
      vcf = AnnotateIncludedRegionsPESR.out,
      vcf_index = AnnotateIncludedRegionsPESR.out_index,
      output_name = "~{batch}.pesr.select_regions",
      select_expression = pesr_region_name + ">=" + pesr_min_inclusion_interval_overlap,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_select_regions_pesr
  }

  call gatk.SelectVariants as SelectRegionsDepth {
    input:
      vcf = AnnotateIncludedRegionsDepth.out,
      vcf_index = AnnotateIncludedRegionsDepth.out_index,
      output_name = "~{batch}.depth.select_regions",
      select_expression = depth_region_name + ">=" + depth_min_inclusion_interval_overlap,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_select_regions_depth
  }

  Array[String] contigs = read_lines(contig_list)
  scatter (contig in contigs) {

    call gatk.SVCluster as Cluster1 {
      input:
        vcfs = [SelectRegionsPESR.out, SelectRegionsDepth.out],
        vcf_indexes = [SelectRegionsPESR.out_index, SelectRegionsDepth.out_index],
        output_name="~{batch}.~{contig}.cluster_1",
        ref_dict=ref_dict,
        contig=contig,
        vid_prefix="~{batch}_~{contig}_c1_",
        omit_members=true,
        algorithm=cluster1_algorithm,
        breakpoint_summary_strategy=cluster1_breakpoint_summary_strategy,
        depth_overlap_fraction=cluster1_depth_overlap_fraction,
        mixed_overlap_fraction=cluster1_mixed_overlap_fraction,
        pesr_overlap_fraction=cluster1_pesr_overlap_fraction,
        depth_breakend_window=cluster1_depth_breakend_window,
        mixed_breakend_window=cluster1_mixed_breakend_window,
        pesr_breakend_window=cluster1_pesr_breakend_window,
        depth_sample_overlap=cluster1_depth_sample_overlap,
        mixed_sample_overlap=cluster1_mixed_sample_overlap,
        pesr_sample_overlap=cluster1_pesr_sample_overlap,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_cluster_1
    }

    call gatk.SVCluster as Cluster2 {
      input:
        vcfs = [Cluster1.out],
        vcf_indexes = [Cluster1.out_index],
        output_name="~{batch}.~{contig}.cluster_2",
        ref_dict=ref_dict,
        vid_prefix="~{batch}_~{contig}_",
        omit_members=true,
        algorithm=cluster2_algorithm,
        breakpoint_summary_strategy=cluster2_breakpoint_summary_strategy,
        depth_overlap_fraction=cluster2_depth_overlap_fraction,
        mixed_overlap_fraction=cluster2_mixed_overlap_fraction,
        pesr_overlap_fraction=cluster2_pesr_overlap_fraction,
        depth_breakend_window=cluster2_depth_breakend_window,
        mixed_breakend_window=cluster2_mixed_breakend_window,
        pesr_breakend_window=cluster2_pesr_breakend_window,
        depth_sample_overlap=cluster2_depth_sample_overlap,
        mixed_sample_overlap=cluster2_mixed_sample_overlap,
        pesr_sample_overlap=cluster2_pesr_sample_overlap,
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_cluster_2
    }

    call gatk.SelectVariants as SelectSize {
      input:
        vcf = Cluster2.out,
        vcf_index = Cluster2.out_index,
        output_name = "~{batch}.select_size",
        select_expression = "(ALGORITHMS==\"depth\" && SVLEN>=" + depth_min_size + ")||(ALGORITHMS!=\"depth\" && SVLEN>=" + pesr_min_size + ")",
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_select_size
    }
  }

  call tasks0506.ConcatVcfs {
    input:
      vcfs = SelectSize.out,
      vcfs_idx = SelectSize.out_index,
      outfile_prefix = "~{batch}.clustered",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_vcfs
  }

  output {
    File out = ConcatVcfs.concat_vcf
    File out_index = ConcatVcfs.concat_vcf_idx
  }
}
