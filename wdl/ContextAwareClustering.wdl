version 1.0

import "CombineSRBothsidePass.wdl" as CombineSRBothsidePass
import "FormatVcfForGatk.wdl" as GatkFormatting
import "TasksClusterBatch.wdl" as ClusterTasks
import "TasksGenotypeBatch.wdl" as GenotypeTasks
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow ContextAwareClustering {
  input {
	File vcf
	String output_prefix
    String variant_prefix

  	# Stratification parameters for each group
  	# Tab-delimited, with columns [NAME, SVTYPE, MIN_SIZE, MAX_SIZE, CONTEXT]
	# First line must be a header with column names (no "#")
	# Lines starting with "#" are ignored
	File stratify_config

	# Clustering parameters for each group in stratify_config
	# Tab-delimited, must have columns [NAME, RECIP_OVERLAP, SIZE_SIMILARITY, BREAKEND_WINDOW, SAMPLE_OVERLAP]
	# Note that one additional group named "not_matched" is required for variants not falling into the defined groups
	# No comments or header allowed
	File clustering_parameters_tsv

	Array[File] sr_bothside_pass_files
	Array[File] sr_background_fail_files

	File ploidy_table

	File reference_fasta
	File reference_fasta_fai
	File reference_dict

	Boolean use_hail = false
	String? gcs_project

	Float? java_mem_fraction

	String gatk_docker
	String sv_base_mini_docker
	String sv_pipeline_docker

	RuntimeAttr? runtime_attr_stratify
    RuntimeAttr? runtime_attr_svcluster
  }

  call SVStratifyTask {
	input:
      vcf=vcf,
      output_prefix="~{output_prefix}.stratify",
      reference_dict=reference_dict,
      java_mem_fraction=java_mem_fraction,
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_stratify
  }

  Array[Array[String]] clustering_params_values = read_tsv(clustering_parameters_tsv)
  #Scatter per group
  scatter (i in range(length(clustering_params_values))) {

	String group_name = clustering_params_values[i][0]
    Float reciprocal_overlap = float(clustering_params_values[i][1])
    Float size_similarity = float(clustering_params_values[i][2])
    Int breakend_window = int(clustering_params_values[i][3])
    Float sample_overlap = float(clustering_params_values[i][4])
	# TODO: disable sorting
	call ClusterTasks.SVCluster {
	  input:
        vcfs=SVStratifyTask.out_vcfs,
        vcf_grep_expression=".~{group_name}.",
        ploidy_table=ploidy_table,
        output_prefix="~{output_prefix}.svcluster.~{group_name}",
        fast_mode=false,
        pesr_sample_overlap=sample_overlap,
        pesr_interval_overlap=reciprocal_overlap,
        pesr_breakend_window=breakend_window,
        pesr_size_similarity=size_similarity,
        depth_sample_overlap=sample_overlap,
        depth_interval_overlap=reciprocal_overlap,
        depth_breakend_window=breakend_window,
        depth_size_similarity=size_similarity,
        mixed_sample_overlap=sample_overlap,
        mixed_interval_overlap=reciprocal_overlap,
        mixed_breakend_window=breakend_window,
        mixed_size_similarity=size_similarity,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{variant_prefix}_~{group_name}",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
	}

	#Subset bothside_pass & background_fail to chromosome of interest
	call SubsetVariantList as SubsetBothsidePass {
	  input:
		vid_list=CombineSRBothsidePass.out,
		vid_col=2,
		vcf=SVCluster.out,
		outfile_name="~{cohort_name}.combine_batches.sr_bothside_pass.~{contig}.subset.list",
		sv_base_mini_docker=sv_base_mini_docker,
		runtime_attr_override=runtime_override_subset_bothside_pass
	}
	call SubsetVariantList as SubsetBackgroundFail {
	  input:
		vid_list=CombineBackgroundFail.outfile,
		vid_col=1,
		vcf=SVCluster.out,
		outfile_name="~{cohort_name}.combine_batches.sr_background_fail.~{contig}.subset.list",
		sv_base_mini_docker=sv_base_mini_docker,
		runtime_attr_override=runtime_override_subset_background_fail
	}

	#Update SR background fail & bothside pass files (1)
	call MiniTasks.UpdateSrList as UpdateBackgroundFail {
	  input:
		vcf=SVCluster.out,
		original_list=SubsetBothsidePass.filtered_vid_list,
		outfile="~{cohort_name}.combine_batches.sr_bothside_pass.~{contig}.list",
		sv_pipeline_docker=sv_pipeline_docker,
		runtime_attr_override=runtime_override_update_sr_list
	}
	call MiniTasks.UpdateSrList as UpdateBothsidePass {
	  input:
		vcf=SVCluster.out,
		original_list=SubsetBackgroundFail.filtered_vid_list,
		outfile="~{cohort_name}.combine_batches.sr_background_fail.~{contig}.list",
		sv_pipeline_docker=sv_pipeline_docker,
		runtime_attr_override=runtime_override_update_sr_list
	}
  }

  # Merge resolved vcfs for QC
  if (merge_vcfs) {
	call MiniTasks.ConcatVcfs {
	  input:
		vcfs=SVCluster.out,
		vcfs_idx=SVCluster.out_index,
		naive=true,
		outfile_prefix="~{cohort_name}.combine_batches",
		sv_base_mini_docker=sv_base_mini_docker,
		runtime_attr_override=runtime_override_concat
	}
  }

  #Final outputs
  output {
	Array[File] combined_vcfs = SVCluster.out
	Array[File] combined_vcf_indexes = SVCluster.out_index
	Array[File] cluster_bothside_pass_lists = UpdateBothsidePass.updated_list
	Array[File] cluster_background_fail_lists = UpdateBackgroundFail.updated_list
	File? combine_batches_merged_vcf = ConcatVcfs.concat_vcf
	File? combine_batches_merged_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

task SVStratifyTask {
  input {
	File vcf
	String output_prefix
	File reference_dict
	String? contig
	String? additional_args

	Float? java_mem_fraction
	String gatk_docker
	RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
	truth_vcf: {
				 localization_optional: true
			   }
	eval_vcf:  {
				 localization_optional: true
			   }
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

  output {
	Array[File] out_vcfs = glob("out/~{output_prefix}.*.vcf.gz")
  	Array[File] out_vcf_indexes = glob("out/~{output_prefix}.*.vcf.gz.tbi")
  }
  command <<<
	set -euo pipefail

	function getJavaMem() {
	# get JVM memory in MiB by getting total memory from /proc/meminfo
	# and multiplying by java_mem_fraction
	cat /proc/meminfo \
	  | awk -v MEM_FIELD="$1" '{
	  	f[substr($1, 1, length($1)-1)] = $2
	  } END {
	  	printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
	  }'
	}
	JVM_MAX_MEM=$(getJavaMem MemTotal)
	echo "JVM memory: $JVM_MAX_MEM"

	mkdir -p out
	gatk --java-options "-Xmx${JVM_MAX_MEM}" SVStratify \
	--sequence-dictionary ~{reference_dict} \
	  ~{"-L " + contig} \
	  -V ~{vcf} \
	  --stratify-config /Users/markw/Work/talkowski/sv-pipe-testing/mw_sv_stratify/sv_stratify_config.rm_sr_sd.tsv \
	  -O out/ \
	  --output-prefix ~{output_prefix} \
	  --context-intervals /Users/markw/Work/talkowski/sv-pipe-testing/mw_sv_stratify/hg38.SegDup.sorted.merged.bed \
	  --context-intervals /Users/markw/Work/talkowski/sv-pipe-testing/mw_sv_stratify/hg38.SimpRep.sorted.merged.bed \
	  --context-intervals /Users/markw/Work/talkowski/sv-pipe-testing/mw_sv_stratify/hg38.RM.sorted.merged.bed \
	  --context-name SD \
	  --context-name SR \
	  --context-name RM
  >>>
  runtime {
	cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
	memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
	disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
	bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
	docker: gatk_docker
	preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
	maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

