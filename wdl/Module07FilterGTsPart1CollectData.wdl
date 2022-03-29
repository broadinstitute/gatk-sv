version 1.0

import "Structs.wdl"
import "CalcAF.wdl" as calcAF
import "MinGQRocOpt.wdl" as roc_opt_sub
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Module07MinGQTasks.wdl" as minGQTasks
import "ReviseSVtypeINStoMEI.wdl" as ReviseSVtype

workflow Module07FilterGTsPart1 {
  input {
    File vcf
    File vcf_idx
    String prefix
    File contiglist
    File trios_famfile
    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1
    Int records_per_shard_AF_annotation = 1000
    String filter_metric = "GQ"
    Array[String] gather_trio_geno_options = []
    File? pcrplus_samples_list

    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_base_docker
    String sv_pipeline_updates_docker

    # overrides for local tasks
    RuntimeAttr? runtime_attr_SplitVcfPerContig
    RuntimeAttr? runtime_attr_GatherTrioData
    RuntimeAttr? runtime_attr_ReviseSVtypeMEI
    RuntimeAttr? runtime_override_combine_step_1_vcfs
    RuntimeAttr? runtime_override_split_vcf_to_clean
    RuntimeAttr? runtime_attr_collect_trio_svdat_pcrminus
    RuntimeAttr? runtime_attr_collect_trio_svdat_pcrplus
  }

  Array[Array[String]] contigs = read_tsv(contiglist)


  # Differenciate MEI from INS
  call ReviseSVtype.ReviseSVtypeINStoMEI as ReviseSVtypeMEI {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      sv_pipeline_updates_docker = sv_pipeline_updates_docker,
      prefix = prefix,
      contiglist = contiglist,
      max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
      min_records_per_shard_step1 = min_records_per_shard_step1,
      runtime_attr_ReviseSVtypeMEI = runtime_attr_ReviseSVtypeMEI,
      runtime_override_split_vcf_to_clean=runtime_override_split_vcf_to_clean,
      runtime_override_combine_step_1_vcfs = runtime_override_combine_step_1_vcfs
  }

  # extract sample names from the vcf, and output a 2-column file with SVID in 1st column and PCR status in 2nd column
  ## we should consider revising this function to use PCR+ and PCR- sample list as input to save mem and disk usage, as vcf can be large
  call minGQTasks.GetSampleLists {
    input:
      vcf = ReviseSVtypeMEI.updated_vcf,
      vcf_idx = ReviseSVtypeMEI.updated_vcf_idx,
      pcrplus_samples_list = pcrplus_samples_list,
      prefix = prefix,
      sv_base_mini_docker = sv_base_mini_docker
  }

  # Shard VCF
  call MiniTasks.ScatterVcf as SplitVcfStep1{
    input:
      vcf = ReviseSVtypeMEI.updated_vcf,
      prefix = prefix,
      records_per_shard = min_records_per_shard_step1,
      sv_pipeline_docker = sv_pipeline_updates_docker,
      runtime_attr_override = runtime_attr_SplitVcfPerContig
  }

  # Split VCF by PCR status and collect AFs separately for each PCR status
  scatter (i in range(length(SplitVcfStep1.shards))) {

    # Shard VCF per-chromosome and add AF annotation
    call calcAF.CalcAF as getAFs {
      input:
        vcf = SplitVcfStep1.shards[i],
        vcf_idx = SplitVcfStep1.shards_idx[i],
        sv_per_shard = records_per_shard_AF_annotation,
        prefix = prefix,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_updates_docker = sv_pipeline_updates_docker
    }

    #Split VCF into PCR+ and PCR-
    call minGQTasks.SplitPcrVcf {
      input:
        vcf = getAFs.vcf_wAFs,
        prefix = prefix,
        pcrplus_samples_list = pcrplus_samples_list,
        sv_base_mini_docker = sv_base_mini_docker
    }

    # Annotate PCR- specific AFs
    call calcAF.CalcAF as getAFs_byPCR {
      input:
        vcf = SplitPcrVcf.PCRMINUS_vcf,
        vcf_idx = SplitPcrVcf.PCRMINUS_vcf_idx,
        sv_per_shard = records_per_shard_AF_annotation,
        prefix = prefix,
        sample_pop_assignments = GetSampleLists.sample_PCR_labels,
        sv_pipeline_docker = sv_pipeline_docker,
        sv_pipeline_updates_docker = sv_pipeline_updates_docker
    }

    # Gather table of AC/AN/AF for PCRPLUS and PCRMINUS samples
    call minGQTasks.GetAfTables {
      input:
        vcf = getAFs_byPCR.vcf_wAFs,
        pcrplus_samples_list = pcrplus_samples_list,
        vcf_idx = getAFs_byPCR.vcf_wAFs_idx,
        prefix = prefix,
        sv_pipeline_docker = sv_pipeline_docker
    }

  }

  # Build tables of AF per variant by PCR status prior to GT filtering
  # This is useful for downstream batch effect detection & other filters
  call minGQTasks.CombineRocOptResults as cat_AF_table_PCRMINUS {
    input:
      shards = GetAfTables.PCRMINUS_AF_table,
      outfile = "~{prefix}.PCRMINUS.AF_preGTFiltering.txt",
      sv_base_mini_docker = sv_base_mini_docker
  }
  if (defined(pcrplus_samples_list)){
    call minGQTasks.CombineRocOptResults as cat_AF_table_PCRPLUS {
      input:
        shards = GetAfTables.PCRPLUS_AF_table,
        outfile = "~{prefix}.PCRPLUS.AF_preGTFiltering.txt",
        sv_base_mini_docker = sv_base_mini_docker,
    }
  }

  # Collect training data for each PCR- trio
  call minGQTasks.SplitFamfile as SplitFamfile_PCRMINUS {
    input:
      vcf = SplitPcrVcf.PCRMINUS_vcf[0],
      vcf_idx = SplitPcrVcf.PCRMINUS_vcf_idx[0],
      famfile = trios_famfile,
      max_count_famfile_shards = 1000,
      fams_per_shard = 1,
      prefix = "~{prefix}.PCRMINUS",
      sv_pipeline_docker = sv_pipeline_docker
  }
  scatter ( fam in SplitFamfile_PCRMINUS.famfile_shards ) {
    call minGQTasks.CollectTrioSVdat as CollectTrioSVdat_PCRMINUS {
      input:
        vcf_shards = SplitPcrVcf.PCRMINUS_vcf,
        famfile = fam,
        filter_metric = filter_metric,
        gather_trio_geno_options=gather_trio_geno_options,
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_collect_trio_svdat_pcrminus
    }
  }
  call minGQTasks.GatherTrioData as GatherTrioData_PCRMINUS {
    input:
      files = CollectTrioSVdat_PCRMINUS.trio_SVdata,
      prefix = "~{prefix}.PCRMINUS",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_GatherTrioData
  }

  # Collect training data for each PCR- trio
  if (defined(pcrplus_samples_list)){
    call minGQTasks.SplitFamfile as SplitFamfile_PCRPLUS {
      input:
        vcf = SplitPcrVcf.PCRPLUS_vcf[0],
        vcf_idx = SplitPcrVcf.PCRPLUS_vcf_idx[0],
        famfile = trios_famfile,
        fams_per_shard = 1,
        max_count_famfile_shards = 1000,
        prefix = "~{prefix}.PCRPLUS",
        sv_pipeline_docker = sv_pipeline_docker
    }
    scatter ( fam in SplitFamfile_PCRPLUS.famfile_shards ) {
      call minGQTasks.CollectTrioSVdat as CollectTrioSVdat_PCRPLUS {
        input:
          vcf_shards = SplitPcrVcf.PCRPLUS_vcf,
          famfile = fam,
          filter_metric = filter_metric,
          gather_trio_geno_options=gather_trio_geno_options,
          sv_pipeline_base_docker = sv_pipeline_base_docker,
          runtime_attr_override = runtime_attr_collect_trio_svdat_pcrplus
      }
    }
    call minGQTasks.GatherTrioData as GatherTrioData_PCRPLUS {
      input:
        files = CollectTrioSVdat_PCRPLUS.trio_SVdata,
        prefix = "~{prefix}.PCRPLUS",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_GatherTrioData
    }
  }


  output {

    Array[File] PCRMINUS_vcfs = SplitPcrVcf.PCRMINUS_vcf
    Array[File] PCRMINUS_vcf_idxs = SplitPcrVcf.PCRMINUS_vcf_idx
    File PCRMINUS_trio_tarball = GatherTrioData_PCRMINUS.tarball
    File PCRMINUS_cleaned_trios_famfile = SplitFamfile_PCRMINUS.cleaned_trios_famfile
    File AF_table_preGTFiltering_PCRMINUS = cat_AF_table_PCRMINUS.merged_file
    
    Array[File]? PCRPLUS_vcfs = SplitPcrVcf.PCRPLUS_vcf
    Array[File]? PCRPLUS_vcf_idxs = SplitPcrVcf.PCRPLUS_vcf_idx
    File? PCRPLUS_trio_tarball = GatherTrioData_PCRPLUS.tarball
    File? PCRPLUS_cleaned_trios_famfile = SplitFamfile_PCRPLUS.cleaned_trios_famfile
    File? AF_table_preGTFiltering_PCRPLUS = cat_AF_table_PCRPLUS.merged_file

  }
}
