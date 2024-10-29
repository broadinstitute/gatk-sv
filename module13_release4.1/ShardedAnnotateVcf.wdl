version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "HailMerge.wdl" as HailMerge
import "AnnotateFunctionalConsequences.wdl" as func
import "PruneAndAddVafs.wdl" as pav
import "AnnotateExternalAF.wdl" as eaf

workflow ShardedAnnotateVcf {

  input {
    File vcf
    File vcf_idx
    String prefix
    String contig

    File? protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? sample_list
    File? ped_file                # Used for M/F AF calculations
    File? par_bed
    File? allosomes_list
    Int   sv_per_shard

    File? ref_bed              # File with external allele frequencies
    String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
    Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

    Boolean use_hail
    String? gcs_project
    Boolean run_fix_ends = false
    Boolean fin_add_pos2 = true

    String sv_pipeline_docker
    String sv_pipeline_hail_docker
    String sv_pipeline_base_docker
    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_fix_ends_rescale_GQ
    RuntimeAttr? runtime_attr_concat_sharded_cluster
    RuntimeAttr? runtime_attr_preconcat_sharded_cluster
    RuntimeAttr? runtime_attr_hail_merge_sharded_cluster
    RuntimeAttr? runtime_attr_fix_header_sharded_cluster
    RuntimeAttr? runtime_attr_get_vcf_header_with_members_info_line
  }

  call MiniTasks.ScatterVcf{
    input:
      vcf = vcf,
      prefix = prefix,
      records_per_shard = sv_per_shard,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_scatter_vcf
  }

  scatter (i in range(length(ScatterVcf.shards))) {

    if (run_fix_ends){
      call FixEndsRescaleGQ {
        input:
          vcf = ScatterVcf.shards[i],
          prefix = "~{prefix}.~{i}",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_fix_ends_rescale_GQ
        }
      }

    File scattered_vcf = select_first([FixEndsRescaleGQ.out, ScatterVcf.shards[i]])
    File scattered_vcf_idx = select_first([FixEndsRescaleGQ.out_idx, ScatterVcf.shards_idx[i]])

    if (fin_add_pos2){
      call FixEndsAndAddPos2{
        input:
          vcf = scattered_vcf, 
          prefix = "~{prefix}.~{i}",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_fix_ends_rescale_GQ
      }
    }

    File scattered_Step2_vcf = select_first([FixEndsAndAddPos2.out,  scattered_vcf])
    File scattered_Step2_vcf_idx = select_first([FixEndsAndAddPos2.out_idx, scattered_vcf_idx])

    call SyncVcfAndIdx{
      input:
        vcf = scattered_Step2_vcf,
        vcf_idx = scattered_Step2_vcf_idx,
        sv_base_mini_docker = sv_base_mini_docker
    }

    if (defined(protein_coding_gtf)){
        call func.AnnotateFunctionalConsequences {
          input:
            vcf = SyncVcfAndIdx.output_vcf,
            vcf_index = SyncVcfAndIdx.output_vcf_idx,
            prefix = "~{prefix}.~{i}",
            protein_coding_gtf = protein_coding_gtf,
            noncoding_bed = noncoding_bed,
            promoter_window = promoter_window,
            max_breakend_as_cnv_length = max_breakend_as_cnv_length,
            additional_args = svannotate_additional_args,
            gatk_docker = gatk_docker,
            runtime_attr_svannotate = runtime_attr_svannotate
        }
    }

    File anno_vcf = select_first([AnnotateFunctionalConsequences.annotated_vcf, SyncVcfAndIdx.output_vcf])
    File anno_vcf_idx = select_first([AnnotateFunctionalConsequences.annotated_vcf_index, SyncVcfAndIdx.output_vcf_idx])

    call pav.PruneAndAddVafs as PruneAndAddVafs {
      input:
        vcf                    = anno_vcf,
        vcf_idx                = anno_vcf_idx,
        prefix                 = prefix,
        contig                 = contig,
        ped_file               = ped_file,
        par_bed                = par_bed,
        allosomes_list         = allosomes_list,
        sample_pop_assignments = sample_pop_assignments,

        sv_base_mini_docker     = sv_base_mini_docker,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_shard_vcf    = runtime_attr_shard_vcf,
        runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
        runtime_attr_combine_vcfs = runtime_attr_combine_vcfs,
        runtime_attr_concat_vcfs  = runtime_attr_concat_vcfs
    }

    if (defined(ref_bed)) {
      call eaf.AnnotateExternalAF as AnnotateExternalAF {
        input:
          vcf     = PruneAndAddVafs.output_vcf,
          vcf_idx = PruneAndAddVafs.output_vcf_idx,
          ref_bed = select_first([ref_bed]),
          population = select_first([population]),
          ref_prefix = select_first([ref_prefix]),
          prefix = prefix,
          contigs = [contig],
          max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
          min_records_per_shard_step1 = min_records_per_shard_step1,
          sv_base_mini_docker = sv_base_mini_docker,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_modify_vcf = runtime_attr_modify_vcf,
          runtime_attr_split_vcf = runtime_attr_split_vcf,
          runtime_attr_combine_vcfs = runtime_attr_combine_vcfs,
          runtime_attr_split_ref_bed = runtime_attr_split_ref_bed,
          runtime_attr_split_query_vcf = runtime_attr_split_query_vcf,
          runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
          runtime_attr_select_matched_svs = runtime_attr_select_matched_svs
      }
    }

  }

  #Array[File?] sharded_annotated_vcf = select_first([AnnotateExternalAF.annotated_vcf, PruneAndAddVafs.output_vcf])
  #Array[File?] sharded_annotated_vcf_idx = select_first([AnnotateExternalAF.annotated_vcf_tbi, PruneAndAddVafs.output_vcf_idx])
  Array[File] sharded_annotated_vcf = PruneAndAddVafs.output_vcf
  Array[File] sharded_annotated_vcf_idx = PruneAndAddVafs.output_vcf_idx


  if (length(sharded_annotated_vcf) == 0) {
    call MiniTasks.GetVcfHeaderWithMembersInfoLine as GetVcfHeader_annotated {
      input:
        vcf_gz=vcf,
        prefix="~{prefix}.annotated",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_get_vcf_header_with_members_info_line
    }
  }

  if (length(sharded_annotated_vcf) > 0) {
    if (use_hail) {
      call HailMerge.HailMerge as ConcatVcfsHail_annotated {
        input:
          vcfs=sharded_annotated_vcf,
          prefix="~{prefix}.annotated",
          gcs_project=gcs_project,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker,
          sv_pipeline_hail_docker=sv_pipeline_hail_docker,
          runtime_attr_preconcat=runtime_attr_preconcat_sharded_cluster,
          runtime_attr_hail_merge=runtime_attr_hail_merge_sharded_cluster,
          runtime_attr_fix_header=runtime_attr_fix_header_sharded_cluster
      }
    }

    if (!use_hail) {
      call MiniTasks.ConcatVcfs as ConcatVcfs_annotated {
        input:
          vcfs=sharded_annotated_vcf,
          vcfs_idx=sharded_annotated_vcf_idx,
          allow_overlaps=true,
          outfile_prefix="~{prefix}.annotated",
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_attr_concat_sharded_cluster
      }
    }

  }


  output {
    File output_vcf = select_first([GetVcfHeader_annotated.out, ConcatVcfs_annotated.concat_vcf, ConcatVcfsHail_annotated.merged_vcf])
    File output_vcf_idx = select_first([GetVcfHeader_annotated.out_idx, ConcatVcfs_annotated.concat_vcf_idx, ConcatVcfsHail_annotated.merged_vcf_index])
  }
}


#function to fix BND, CTX, CPX, INS that have END and END2 represent the breakpoint on the 2nd chromosome
#Note: this is a temp function for the first beta version of gnomad SV callset. It'll be revised and added as part of the manunal revise / clean up script
task FixEndsRescaleGQ {
  input {
    File vcf
    String prefix

    Boolean? fix_ends
    Boolean? rescale_gq

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String outfile = "~{prefix}.vcf.gz"
  Boolean fix_ends_ = select_first([fix_ends, true])
  Boolean rescale_gq_ = select_first([rescale_gq, true])

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<

    set -euo pipefail

    python <<CODE
    import pysam
    import argparse
    from math import floor


    GQ_FIELDS = ["GQ", "PE_GQ", "SR_GQ", "RD_GQ"]

    filts_for_info = 'PESR_GT_OVERDISPERSION HIGH_SR_BACKGROUND BOTHSIDES_SUPPORT VARIABLE_ACROSS_BATCHES'.split(' ')
    filts_to_remove = 'HIGH_PCRPLUS_NOCALL_RATE HIGH_PCRMINUS_NOCALL_RATE'.split(' ')
    filts_to_remove = filts_to_remove + filts_for_info

    def fix_bad_end(record):
      # pysam converts to 0-based half-open intervals by subtracting 1 from start, but END is unaltered
      if record.stop < record.start + 2:
        if record.info["SVTYPE"] == "BND" or record.info["SVTYPE"] == "CTX":
          record.info["END2"] = record.stop  # just in case it is not already set. not needed for INS or CPX
        record.stop = record.start + 1

    def rescale_gq(record):
      for sample in record.samples:
        for gq_field in GQ_FIELDS:
          if gq_field in record.samples[sample] and record.samples[sample][gq_field] is not None:
            record.samples[sample][gq_field] = floor(record.samples[sample][gq_field] / 10)


    with pysam.VariantFile("~{vcf}", 'r') as f_in, pysam.VariantFile("~{outfile}", 'w', header=f_in.header) as f_out:
      for record in f_in:
        newfilts = [filt for filt in record.filter if filt not in filts_to_remove]
        record.filter.clear()
        for filt in newfilts:
            record.filter.add(filt)
        if len(record.filter) == 0:
            record.filter.add('PASS')
        if "~{fix_ends_}" == "true":
          fix_bad_end(record)
        if "~{rescale_gq_}" == "true":
          rescale_gq(record)
        f_out.write(record)

    CODE
    tabix ~{outfile}

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task FixEndsAndAddPos2 {
  input {
    File vcf
    String prefix

    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String outfile = "~{prefix}.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<

    set -euo pipefail

    python <<CODE

    import os
    import pysam
    import argparse

    def fix_bad_end(record):
      if record.info["SVTYPE"] == "BND":
        record.info['POS2'] = record.info['END2'] -1 
      elif record.info["SVTYPE"] == "CTX":
        record.info["POS2"] = record.info["END2"]-1
      elif record.info["SVTYPE"] == "INS":
        ins_section = [i for i in record.info['SOURCE'].split(',') if "INS_" in i or "DUP_" in i][0]
        [pos2, end2] = [int(i) for i in ins_section.split(':')[1].split('-')]
        record.info["POS2"] = pos2
        record.info["END2"] = end2
      elif record.info["SVTYPE"] == "CPX":
        dup_section = record.info['CPX_INTERVALS'][0]
        [pos2, end2] = [int(i) for i in dup_section.split(':')[1].split('-')]
        record.info["POS2"] = pos2
        record.info["END2"] = end2
      else:
        print(record.info["SVTYPE"]) 
    
    f_in = pysam.VariantFile("~{vcf}", 'r')
    header = f_in.header
    if not "POS2" in header.info.keys():
        info = ('##INFO=<ID=POS2,Number=1,Type=Integer,'
                'Description="Start position of the structural variant on CHR2">')
        header.add_line(info)

    f_out = pysam.VariantFile("~{outfile}", 'w', header = header)
    
    for record in f_in:
      if "END2" in record.info.keys():
        fix_bad_end(record)
      f_out.write(record)

    f_in.close()
    f_out.close()

    CODE
    tabix ~{outfile}

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task SyncVcfAndIdx {
  input {
    File vcf
    File vcf_idx
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String prefix = basename(vcf, ".vcf.gz")

  output {
    File output_vcf = "~{prefix}.vcf.gz"
    File output_vcf_idx = "~{prefix}.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail

    cp ~{vcf} "~{prefix}.vcf.gz"
    cp ~{vcf_idx} "~{prefix}.vcf.gz.tbi"
    echo "nothing else to do here"

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

