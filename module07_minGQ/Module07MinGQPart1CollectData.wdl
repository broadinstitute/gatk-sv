version 1.0

import "Structs.wdl"
import "CalcAF.wdl" as calcAF
import "MinGQRocOpt.wdl" as roc_opt_sub
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Module07MinGQTasks.wdl" as minGQTasks
import "ReviseSVtypeINStoMEI.wdl" as ReviseSVtype

workflow Module07MinGQPart1 {
  input {
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_updates_docker
    File vcf
    File vcf_idx
    String prefix
    File contiglist
    File trios_famfile
    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1
    File? pcrminus_samples_list
    File? pcrplus_samples_list
    File? pcrminus_filter_lookup_table
    File? pcrplus_filter_lookup_table
    Boolean MingqTraining=! defined(pcrminus_filter_lookup_table)

    # overrides for local tasks
    RuntimeAttr? runtime_attr_SplitVcfPerContig
    RuntimeAttr? runtime_attr_CombineVcfs
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


  # collect variants on a specific chromosome, and scatter them into smaller shards
  call MiniTasks.ScatterVcf as SplitVcfStep1{
    input:
      vcf = ReviseSVtypeMEI.updated_vcf,
      prefix = prefix,
      records_per_shard = min_records_per_shard_step1,
      sv_pipeline_docker = sv_pipeline_updates_docker,
      runtime_attr_override = runtime_attr_SplitVcfPerContig
  }

  scatter (i in range(length(SplitVcfStep1.shards))) {
    # Shard VCF per-chromosome and add AF annotation
    call calcAF.CalcAF as getAFs {
      input:
        vcf     = SplitVcfStep1.shards[i],
        vcf_idx = SplitVcfStep1.shards_idx[i],
        sv_per_shard=1000,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_updates_docker=sv_pipeline_updates_docker
    }

    #Split VCF into PCR+ and PCR-
    call minGQTasks.SplitPcrVcf {
      input:
        vcf=getAFs.vcf_wAFs,
        prefix=prefix,
        pcrplus_samples_list=pcrplus_samples_list,
        sv_base_mini_docker=sv_base_mini_docker
    }

    # Dev note Feb 18 2021: the output from cat_AF_table_PCRMINUS is a required
    # input to Module07XfBatchEffect.wdl, so the subsequent three tasks always 
    # need to be generated (even if passing a precomputed minGQ cutoff table)

    # Annotate PCR- specific AFs
    call calcAF.CalcAF as getAFs_byPCR {
      input:
        vcf     = SplitPcrVcf.PCRMINUS_vcf,
        vcf_idx = SplitPcrVcf.PCRMINUS_vcf_idx,
        sv_per_shard=1000,
        prefix=prefix,
        sample_pop_assignments=GetSampleLists.sample_PCR_labels,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_updates_docker=sv_pipeline_updates_docker
    }
    # Gather table of AC/AN/AF for PCRPLUS and PCRMINUS samples
    call minGQTasks.GetAfTables {
      input:
        vcf=getAFs_byPCR.vcf_wAFs,
        pcrplus_samples_list=pcrplus_samples_list,
        vcf_idx=getAFs_byPCR.vcf_wAFs_idx,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  call minGQTasks.CombineRocOptResults as cat_AF_table_PCRMINUS {
    input:
      shards=GetAfTables.PCRMINUS_AF_table,
      outfile="~{prefix}.PCRMINUS.AF_preMinGQ.txt",
      sv_base_mini_docker=sv_base_mini_docker,
  }

  call MiniTasks.ConcatVcfs as CombineVcfs_PCRMINUS {
    input:
      vcfs=SplitPcrVcf.PCRMINUS_vcf,
      naive=true,
      outfile_prefix="~{prefix}.PCRMINUS",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_CombineVcfs
  }

  if (defined(pcrplus_samples_list)){
    call minGQTasks.CombineRocOptResults as cat_AF_table_PCRPLUS {
      input:
        shards=GetAfTables.PCRPLUS_AF_table,
        outfile="~{prefix}.PCRPLUS.AF_preMinGQ.txt",
        sv_base_mini_docker=sv_base_mini_docker,
    }

    call MiniTasks.ConcatVcfs as CombineVcfs_PCRPLUS {
      input:
        vcfs=SplitPcrVcf.PCRPLUS_vcf,
        naive=true,
        outfile_prefix="~{prefix}.PCRPLUS",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_CombineVcfs
    }
  }

  if (MingqTraining) {
    ### PCRMINUS
    call minGQTasks.SplitFamfile as SplitFamfile_PCRMINUS {
      input:
        vcf=SplitPcrVcf.PCRMINUS_vcf[0],
        vcf_idx=SplitPcrVcf.PCRMINUS_vcf_idx[0],
        famfile=trios_famfile,
        max_count_famfile_shards=1000,
        fams_per_shard=1,
        prefix="~{prefix}.PCRMINUS",
        sv_pipeline_docker=sv_pipeline_docker
    }
    scatter ( fam in SplitFamfile_PCRMINUS.famfile_shards ) {
      call minGQTasks.CollectTrioSVdat as CollectTrioSVdat_PCRMINUS {
        input:
          vcf_shards=SplitPcrVcf.PCRMINUS_vcf,
          famfile=fam,
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override = runtime_attr_collect_trio_svdat_pcrminus
      }
    }
    call minGQTasks.GatherTrioData as GatherTrioData_PCRMINUS {
      input:
        files=CollectTrioSVdat_PCRMINUS.trio_SVdata,
        prefix="~{prefix}.PCRMINUS",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_GatherTrioData
    }

    ### PCRPLUS
    if (defined(pcrplus_samples_list)){
      call minGQTasks.SplitFamfile as SplitFamfile_PCRPLUS {
        input:
          vcf=SplitPcrVcf.PCRPLUS_vcf[0],
          vcf_idx=SplitPcrVcf.PCRPLUS_vcf_idx[0],
          famfile=trios_famfile,
          fams_per_shard=1,
          max_count_famfile_shards = 1000,
          prefix="~{prefix}.PCRPLUS",
          sv_pipeline_docker=sv_pipeline_docker
      }
      scatter ( fam in SplitFamfile_PCRPLUS.famfile_shards ) {
        call minGQTasks.CollectTrioSVdat as CollectTrioSVdat_PCRPLUS {
          input:
            vcf_shards=SplitPcrVcf.PCRPLUS_vcf,
            famfile=fam,
            sv_pipeline_docker=sv_pipeline_docker,
            runtime_attr_override = runtime_attr_collect_trio_svdat_pcrplus
        }
      }
      call minGQTasks.GatherTrioData as GatherTrioData_PCRPLUS {
        input:
          files=CollectTrioSVdat_PCRPLUS.trio_SVdata,
          prefix="~{prefix}.PCRPLUS",
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_attr_GatherTrioData
      }
    }
  }


  # Final output
  output {

    Array[File] PCRMINUS_vcf_lists = SplitPcrVcf.PCRMINUS_vcf
    Array[File] PCRMINUS_vcf_idx_lists = SplitPcrVcf.PCRMINUS_vcf_idx
    Array[File] PCRPLUS_vcf_lists = SplitPcrVcf.PCRPLUS_vcf
    Array[File] PCRPLUS_vcf_idx_lists = SplitPcrVcf.PCRPLUS_vcf_idx

    File? AF_table_preMinGQ_PCRMINUS = cat_AF_table_PCRMINUS.merged_file
    File? AF_table_preMinGQ_PCRPLUS = cat_AF_table_PCRPLUS.merged_file

    File? PCRMINUS_trio_tarball = GatherTrioData_PCRMINUS.tarball
    File? PCRPLUS_trio_tarball = GatherTrioData_PCRPLUS.tarball
    File? PCRMINUS_cleaned_trios_famfile = SplitFamfile_PCRMINUS.cleaned_trios_famfile
    File? PCRPLUS_cleaned_trios_famfile = SplitFamfile_PCRPLUS.cleaned_trios_famfile
  }

}

