version 1.0

import "Structs.wdl"
import "SVConcordancePacBioSample.wdl" as concordance
import "MakeGqRecalibratorTrainingSetFromPacBio.wdl" as makegq
import "Utils.wdl" as utils
import "TasksMakeCohortVcf.wdl" as mini_tasks

workflow MakeGqRecalibratorTrainingSetFromMultisampleVcf {

  input {
    # Cleaned GATK-formatted vcf
    # SVConcordance should be run first if the training set is a proper subset of the cohort
    # This can be either a single whole-genome vcf or multiple vcf shards.
    # Assumes all vcfs have indexes, i.e. at {VCF_PATH}.tbi
    Array[File] vcfs

    File training_sample_ids  # Sample IDs with PacBio or array data
    String? output_prefix
    File ploidy_table

    Array[String] pacbio_sample_ids  # Corresponding to files below (must be a subset of training_sample_ids)
    File truth_vcf
    Array[File] vapor_files

    # Optional: array intensity ratio files
    Array[File]? irs_sample_batches
    Array[File]? irs_test_reports
    File? irs_contigs_fai

    Boolean write_detail_report = true
    Int? vapor_max_cnv_size = 5000
    Float? vapor_min_precision = 0.999
    Int? vapor_pos_read_threshold = 2
    Int? irs_min_cnv_size = 10000
    Float? irs_good_pvalue_threshold = 0.000001
    Int? irs_min_probes = 5

    Float pesr_interval_overlap_strict = 0.1
    Float pesr_size_similarity_strict = 0.5
    Int pesr_breakend_window_strict = 5000
    File? clustering_config_strict
    File? stratification_config_strict

    Float pesr_interval_overlap_loose = 0
    Float pesr_size_similarity_loose = 0
    Int pesr_breakend_window_loose = 5000
    File? clustering_config_loose
    File? stratification_config_loose

    # These arrays give the names and intervals for reference contexts for stratification (same lengths)
    # Names must correspond to those in the stratification config files
    Array[String]? track_names
    Array[File]? track_bed_files

    File reference_dict

    String sv_utils_docker
    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    String linux_docker
  }

  String output_prefix_ =
    if defined(output_prefix) then
      select_first([output_prefix])
    else
        basename(vcfs[0], ".vcf.gz")

  scatter (i in range(length(vcfs))) {
    call utils.SubsetVcfBySamplesList as SubsetTrainingSamples {
      input:
        vcf = vcfs[i],
        vcf_idx = vcfs[i] + ".tbi",
        list_of_samples = training_sample_ids,
        outfile_name = "~{output_prefix_}.training_samples.shard_~{i}.vcf.gz",
        remove_samples = false,
        remove_private_sites = true,
        keep_af = true,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call mini_tasks.ConcatVcfs as ConcatTrainingSampleVcfs {
    input:
      vcfs=SubsetTrainingSamples.vcf_subset,
      vcfs_idx=SubsetTrainingSamples.vcf_subset_index,
      naive=true,
      outfile_prefix="~{output_prefix_}.concat_training_sample_vcfs",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call utils.WriteLines as WritePacBioSampleIds {
    input:
      lines = pacbio_sample_ids,
      output_filename = "~{output_prefix_}.pacbio_sample_ids.list",
      linux_docker = linux_docker
  }

  call makegq.GetVariantListsFromVaporAndIRS {
    input:
      vcf=ConcatTrainingSampleVcfs.concat_vcf,
      output_prefix=output_prefix_,
      vapor_sample_ids=pacbio_sample_ids,
      vapor_files=vapor_files,
      irs_sample_batches=select_first([irs_sample_batches, []]),
      irs_test_reports=select_first([irs_test_reports, []]),
      irs_contigs_fai=irs_contigs_fai,
      vapor_max_cnv_size=vapor_max_cnv_size,
      vapor_min_precision=vapor_min_precision,
      vapor_pos_read_threshold=vapor_pos_read_threshold,
      irs_min_cnv_size=irs_min_cnv_size,
      irs_good_pvalue_threshold=irs_good_pvalue_threshold,
      irs_min_probes=irs_min_probes,
      sv_utils_docker=sv_utils_docker
  }

  scatter (i in range(length(vcfs))) {
    call utils.SubsetVcfBySamplesList {
      input:
        vcf = vcfs[i],
        vcf_idx = vcfs[i] + ".tbi",
        list_of_samples = WritePacBioSampleIds.out,
        outfile_name = "~{output_prefix_}.pacbio_samples.shard_~{i}.vcf.gz",
        remove_samples = false,
        remove_private_sites = true,
        keep_af = true,
        sv_base_mini_docker = sv_base_mini_docker
    }
  }

  call mini_tasks.ConcatVcfs as ConcatPacbioSampleVcfs {
    input:
      vcfs=SubsetVcfBySamplesList.vcf_subset,
      vcfs_idx=SubsetVcfBySamplesList.vcf_subset_index,
      naive=true,
      outfile_prefix="~{output_prefix_}.concat_pacbio_sample_vcfs",
      sv_base_mini_docker=sv_base_mini_docker
  }

  scatter (i in range(length(pacbio_sample_ids))) {
    call makegq.PrepSampleVcf {
      input:
        sample_id=pacbio_sample_ids[i],
        vcf=ConcatPacbioSampleVcfs.concat_vcf,
        vcf_index=ConcatPacbioSampleVcfs.concat_vcf_idx,
        output_prefix=output_prefix_,
        sv_pipeline_docker=sv_pipeline_docker
    }
    call utils.SubsetVcfBySample as SubsetPacBioSampleTruth {
      input:
        vcf = vcfs[i],
        vcf_idx = vcfs[i] + ".tbi",
        sample_id = pacbio_sample_ids[i],
        outfile_name = "~{output_prefix_}.pacbio_samples.shard_~{i}.vcf.gz",
        remove_samples = false,
        remove_private_sites = true,
        keep_af = true,
        sv_base_mini_docker = sv_base_mini_docker
    }
    call concordance.SVConcordancePacBioSample {
      input:
        sample_id=pacbio_sample_ids[i],
        sample_vcf=PrepSampleVcf.out,
        pacbio_sample_vcfs=SubsetPacBioSampleTruth,
        tool_names=["truth"],
        prefix="~{output_prefix_}.loose",
        ploidy_table=ploidy_table,
        reference_dict=reference_dict,
        pesr_interval_overlap=pesr_interval_overlap_loose,
        pesr_size_similarity=pesr_size_similarity_loose,
        pesr_breakend_window=pesr_breakend_window_loose,
        clustering_config=clustering_config_loose,
        stratification_config=stratification_config_loose,
        track_names=track_names,
        track_bed_files=track_bed_files,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_docker=sv_pipeline_docker,
        gatk_docker=gatk_docker,
        linux_docker=linux_docker
    }
    call makegq.RefineSampleLabels {
      input:
        sample_id=pacbio_sample_ids[i],
        main_vcf=PrepSampleVcf.out,
        vapor_json=GetVariantListsFromVaporAndIRS.vapor_json,
        tool_names=["truth"],
        loose_concordance_vcfs=SVConcordancePacBioSample.pacbio_concordance_vcfs,
        strict_concordance_vcfs=SVConcordancePacBioSample.pacbio_concordance_vcfs,
        output_prefix="~{output_prefix_}.gq_training_labels.~{pacbio_sample_ids[i]}",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  call makegq.MergeJsons {
    input:
      jsons=flatten([[GetVariantListsFromVaporAndIRS.irs_json], RefineSampleLabels.out_json]),
      output_prefix="~{output_prefix_}.gq_training_labels",
      sv_pipeline_docker=sv_pipeline_docker
  }

  call mini_tasks.ConcatHeaderedTextFiles {
    input:
      text_files=RefineSampleLabels.out_table,
      gzipped=false,
      output_filename="~{output_prefix_}.gq_training_labels.tsv",
      linux_docker=linux_docker
  }

  output {
    File gq_recalibrator_training_json = MergeJsons.out
    File pacbio_support_summary_table = ConcatHeaderedTextFiles.out

    File training_sample_vcf = ConcatTrainingSampleVcfs.concat_vcf
    File training_sample_vcf_index = ConcatTrainingSampleVcfs.concat_vcf_idx
    File pacbio_sample_vcf = ConcatPacbioSampleVcfs.concat_vcf
    File pacbio_sample_vcf_index = ConcatPacbioSampleVcfs.concat_vcf_idx

    File vapor_output_json = GetVariantListsFromVaporAndIRS.vapor_json
    File irs_output_json = GetVariantListsFromVaporAndIRS.irs_json
    Array[File] pacbio_concordance_vcf_tars = SVConcordancePacBioSample.pacbio_concordance_vcfs_tar
  }
}

