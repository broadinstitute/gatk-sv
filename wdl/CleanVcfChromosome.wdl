version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "FormatVcfForGatk.wdl" as fvcf
import "HailMerge.wdl" as HailMerge

workflow CleanVcfChromosome {
	input {
		File vcf
		String contig
		File background_list
		File ped_file
		File allosome_fai
		String prefix
		Int max_shards_per_chrom_step1
		File bothsides_pass_list
		Int min_records_per_shard_step1
		Int samples_per_step2_shard
		File? outlier_samples_list
		Int? max_samples_per_shard_step3

		File HERVK_reference
		File LINE1_reference

		File ploidy_table
		String chr_x
		String chr_y

		File? svtk_to_gatk_script  # For debugging

		Boolean use_hail
		String? gcs_project

		String gatk_docker
		String linux_docker
		String sv_base_mini_docker
		String sv_pipeline_docker

		# overrides for local tasks
		RuntimeAttr? runtime_override_clean_vcf_1a
		RuntimeAttr? runtime_override_clean_vcf_1b
		RuntimeAttr? runtime_override_clean_vcf_2
		RuntimeAttr? runtime_override_clean_vcf_3
		RuntimeAttr? runtime_override_clean_vcf_4
		RuntimeAttr? runtime_override_clean_vcf_5
		RuntimeAttr? runtime_override_stitch_fragmented_cnvs
		RuntimeAttr? runtime_override_final_cleanup
		RuntimeAttr? runtime_override_rescue_me_dels
		RuntimeAttr? runtime_attr_add_high_fp_rate_filters

		RuntimeAttr? runtime_override_preconcat_step1
		RuntimeAttr? runtime_override_hail_merge_step1
		RuntimeAttr? runtime_override_fix_header_step1

		RuntimeAttr? runtime_override_preconcat_drc
		RuntimeAttr? runtime_override_hail_merge_drc
		RuntimeAttr? runtime_override_fix_header_drc

		# overrides for MiniTasks
		RuntimeAttr? runtime_override_split_vcf_to_clean
		RuntimeAttr? runtime_override_split_include_list
		RuntimeAttr? runtime_override_combine_clean_vcf_2
		RuntimeAttr? runtime_override_drop_redundant_cnvs
		RuntimeAttr? runtime_override_combine_step_1_vcfs
		RuntimeAttr? runtime_override_sort_drop_redundant_cnvs
		RuntimeAttr? runtime_attr_format
	}

	call fvcf.FormatVcf as FormatVcfToClean {
		input:
			vcf=vcf,
			ploidy_table=ploidy_table,
			output_prefix="~{prefix}.formatted",
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_attr_format
	}

	call MiniTasks.SplitVcf as SplitVcfToClean {
		input:
			vcf=FormatVcfToClean.out,
			contig=contig,
			prefix="~{prefix}.shard_",
			n_shards=max_shards_per_chrom_step1,
			min_vars_per_shard=min_records_per_shard_step1,
			sv_base_mini_docker=sv_base_mini_docker,
			runtime_attr_override=runtime_override_split_vcf_to_clean
	}

	scatter ( i in range(length(SplitVcfToClean.vcf_shards)) ) {
		call CleanVcf1a {
			input:
				vcf=SplitVcfToClean.vcf_shards[i],
				prefix="~{prefix}.clean_vcf_1a.shard_~{i}",
				background_fail_list=background_list,
				bothsides_pass_list=bothsides_pass_list,
				ped_file=ped_file,
				allosome_fai=allosome_fai,
				chr_x=chr_x,
				chr_y=chr_y,
				gatk_docker=gatk_docker,
				runtime_attr_override=runtime_override_clean_vcf_1a
		}
	}

	if (use_hail) {
		call HailMerge.HailMerge as CombineStep1VcfsHail {
			input:
				vcfs=CleanVcf1a.intermediate_vcf,
				prefix="~{prefix}.combine_step_1_vcfs",
				gcs_project=gcs_project,
				sv_base_mini_docker=sv_base_mini_docker,
				sv_pipeline_docker=sv_pipeline_docker,
				runtime_override_preconcat=runtime_override_preconcat_step1,
				runtime_override_hail_merge=runtime_override_hail_merge_step1,
				runtime_override_fix_header=runtime_override_fix_header_step1
		}
	}
	if (!use_hail) {
		call MiniTasks.ConcatVcfs as CombineStep1Vcfs {
			input:
				vcfs=CleanVcf1a.intermediate_vcf,
				vcfs_idx=CleanVcf1a.intermediate_vcf_idx,
				naive=true,
				generate_index=false,
				outfile_prefix="~{prefix}.combine_step_1_vcfs",
				sv_base_mini_docker=sv_base_mini_docker,
				runtime_attr_override=runtime_override_combine_step_1_vcfs
		}
	}

	call CleanVcf1b {
		input:
			vcf=select_first([CombineStep1Vcfs.concat_vcf, CombineStep1VcfsHail.merged_vcf]),
			prefix="~{prefix}.clean_vcf_1b",
			gatk_docker=gatk_docker,
			runtime_attr_override=runtime_override_clean_vcf_1b
	}

	call MiniTasks.SplitUncompressed as SplitIncludeList {
		input:
			whole_file=CleanVcf1a.include_list[0],
			lines_per_shard=samples_per_step2_shard,
			shard_prefix="~{prefix}.split_include_list.",
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_override_split_include_list
	}

	scatter ( i in range(length(SplitIncludeList.shards)) ){
		call CleanVcf2 {
			input:
				vcf=CleanVcf1b.out,
				prefix="~{prefix}.clean_vcf_2.shard_~{i}",
				include_list=SplitIncludeList.shards[i],
				gatk_docker=gatk_docker,
				runtime_attr_override=runtime_override_clean_vcf_2
			}
	}

	call MiniTasks.CatUncompressedFiles as CombineCleanVcf2 {
		input:
			shards=CleanVcf2.out,
			outfile_name="~{prefix}.combine_clean_vcf_2.txt",
			sv_base_mini_docker=sv_base_mini_docker,
			runtime_attr_override=runtime_override_combine_clean_vcf_2
	}

	call CleanVcf3 {
		input:
			rd_cn_revise=CombineCleanVcf2.outfile,
			max_samples_shard = max_samples_per_shard_step3,
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_override_clean_vcf_3
	}

	scatter ( i in range(length(CleanVcf3.shards)) ){
		call CleanVcf4 {
			input:
				vcf=CleanVcf1b.out,
				prefix="~{prefix}.clean_vcf_4.shard_~{i}",
				outlier_samples_list=outlier_samples_list,
				rd_cn_revise=CleanVcf3.shards[i],
				gatk_docker=gatk_docker,
				runtime_attr_override=runtime_override_clean_vcf_4
		}
	}

	if (use_hail) {
		call HailMerge.HailMerge as CombineStep4VcfsHail {
			input:
				vcfs=CleanVcf4.out,
				prefix="~{prefix}.combine_revised_4",
				gcs_project=gcs_project,
				sv_base_mini_docker=sv_base_mini_docker,
				sv_pipeline_docker=sv_pipeline_docker,
				runtime_override_preconcat=runtime_override_preconcat_step1,
				runtime_override_hail_merge=runtime_override_hail_merge_step1,
				runtime_override_fix_header=runtime_override_fix_header_step1
		}
	}
	if (!use_hail) {
		call MiniTasks.ConcatVcfs as CombineStep4Vcfs {
			input:
				vcfs=CleanVcf4.out,
				vcfs_idx=CleanVcf4.out_idx,
				naive=true,
				generate_index=true,
				outfile_prefix="~{prefix}.combine_revised_4",
				sv_base_mini_docker=sv_base_mini_docker,
				runtime_attr_override=runtime_override_combine_step_1_vcfs
		}
	}

	call CleanVcf5 {
		input:
			vcf=select_first([CombineStep4Vcfs.concat_vcf, CombineStep4VcfsHail.merged_vcf]),
			prefix="~{prefix}.clean_vcf_5",
			gatk_docker=gatk_docker,
			runtime_attr_override=runtime_override_clean_vcf_5
	}

	call DropRedundantCnvs {
		input:
			vcf=CleanVcf5.out,
			prefix="~{prefix}.drop_redundant_cnvs",
			contig=contig,
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_override_drop_redundant_cnvs
	}

	if (use_hail) {
		call HailMerge.HailMerge as SortDropRedundantCnvsHail {
			input:
				vcfs=[DropRedundantCnvs.out],
				prefix="~{prefix}.drop_redundant_cnvs.sorted",
				gcs_project=gcs_project,
				reset_cnv_gts=true,
				sv_base_mini_docker=sv_base_mini_docker,
				sv_pipeline_docker=sv_pipeline_docker,
				runtime_override_preconcat=runtime_override_preconcat_drc,
				runtime_override_hail_merge=runtime_override_hail_merge_drc,
				runtime_override_fix_header=runtime_override_fix_header_drc
		}
	}
	if (!use_hail) {
		call MiniTasks.SortVcf as SortDropRedundantCnvs {
			input:
				vcf=DropRedundantCnvs.out,
				outfile_prefix="~{prefix}.drop_redundant_cnvs.sorted",
				sv_base_mini_docker=sv_base_mini_docker,
				runtime_attr_override=runtime_override_sort_drop_redundant_cnvs
		}
	}

	call StitchFragmentedCnvs {
		input:
			vcf=select_first([SortDropRedundantCnvs.out, SortDropRedundantCnvsHail.merged_vcf]),
			prefix="~{prefix}.stitch_fragmented_cnvs",
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_override_stitch_fragmented_cnvs
	}

	call RescueMobileElementDeletions {
		input:
		vcf = StitchFragmentedCnvs.stitched_vcf_shard,
		prefix = "~{prefix}.rescue_me_dels",
		LINE1 = LINE1_reference,
		HERVK = HERVK_reference,
		sv_pipeline_docker = sv_pipeline_docker,
		runtime_attr_override = runtime_override_rescue_me_dels
  }

  call AddHighFDRFilters {
    input:
      vcf=RescueMobileElementDeletions.out,
      prefix="~{prefix}.high_fdr_filtered",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_add_high_fp_rate_filters
  }

  call FinalCleanup {
    input:
      vcf=AddHighFDRFilters.out,
      contig=contig,
      prefix="~{prefix}.final_cleanup",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_final_cleanup
  }

  call fvcf.FormatVcf as FormatVcfToOutput {
    input:
      vcf=FinalCleanup.final_cleaned_shard,
      ploidy_table=ploidy_table,
      args="--scale-down-gq",
      output_prefix="~{prefix}.final_format",
      script=svtk_to_gatk_script,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_format
  }
	
	output {
		File out = FormatVcfToOutput.out
		File out_idx = FormatVcfToOutput.out_index
	}
}


task CleanVcf1a {
	input {
		File vcf
		String prefix
		File background_fail_list
		File bothsides_pass_list
		File ped_file
		File allosome_fai
		String chr_x
		String chr_y
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size([vcf, background_fail_list, bothsides_pass_list], "GB")
	RuntimeAttr runtime_default = object {
																	mem_gb: 3.75,
																	disk_gb: ceil(10.0 + input_size * 2),
																	cpu_cores: 1,
																	preemptible_tries: 3,
																	max_retries: 1,
																	boot_disk_gb: 10
																}						
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: gatk_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
	String output_vcf = "~{prefix}.vcf.gz"
	String output_samples_list = "~{prefix}.includelist.txt"

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVCleanPt1a \
			-V ~{vcf} \
			-O ~{output_vcf} \
			--fail-list ~{background_fail_list} \
			--pass-list ~{bothsides_pass_list} \
			--chr-X ~{chr_x} \
			--chr-Y ~{chr_y} \
			--output-samples-list ~{output_samples_list}
	>>>

	output {
		File include_list="~{output_samples_list}"
		File intermediate_vcf="~{output_vcf}"
		File intermediate_vcf_idx="~{output_vcf}.tbi"
	}
}

task CleanVcf1b {
	input {
		File vcf
		String prefix
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size([vcf], "GB")
	RuntimeAttr runtime_default = object {
																	mem_gb: 3.75,
																	disk_gb: ceil(10.0 + input_size * 2),
																	cpu_cores: 1,
																	preemptible_tries: 3,
																	max_retries: 1,
																	boot_disk_gb: 10
																}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: gatk_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	String output_vcf = "~{prefix}.vcf.gz"
	Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVCleanPt1b \
			-V ~{vcf} \
			-O ~{output_vcf}
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}

task CleanVcf2 {
	input {
		File vcf
		String prefix
		File include_list
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size([vcf, include_list], "GB")
	Float base_disk_gb = 10.0
	Float input_disk_scale = 3.0
	RuntimeAttr runtime_default = object {
		mem_gb: 2.0,
		disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
		cpu_cores: 1,
		preemptible_tries: 3,
		max_retries: 1,
		boot_disk_gb: 10
	}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: gatk_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	String output_revised_list = "~{prefix}.txt"
	Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVCleanPt2 \
			-V ~{vcf} \
			--sample-list ~{include_list} \
			--output-revised-list ~{output_revised_list}
	>>>

	output {
		File out="~{output_revised_list}"
	}
}


task CleanVcf3 {
	input {
		File rd_cn_revise
		Int? max_samples_shard
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}
	
	Int max_samples_shard_ = select_first([max_samples_shard, 7000])
	Float input_size = size(rd_cn_revise, "GB")
	RuntimeAttr runtime_default = object {
		mem_gb: 3.75,
		disk_gb: ceil(10.0 + input_size * 2.0),
		cpu_cores: 1,
		preemptible_tries: 3,
		max_retries: 1,
		boot_disk_gb: 10
	}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: sv_pipeline_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	command <<<
		set -euo pipefail
		python /opt/sv-pipeline/04_variant_resolution/scripts/clean_vcf_part3.py ~{rd_cn_revise} -s ~{max_samples_shard_}
		# Ensure there is at least one shard
		touch shards/out.0_0.txt
	>>>

	output {
		 Array[File] shards = glob("shards/*")
	}
}


task CleanVcf4 {
	input {
		File vcf
		String prefix
		File rd_cn_revise
		File? outlier_samples_list
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size([vcf, rd_cn_revise], "GB")
	RuntimeAttr runtime_default = object {
																	mem_gb: 2.0,
																	disk_gb: 50,
																	cpu_cores: 1,
																	preemptible_tries: 3,
																	max_retries: 1,
																	boot_disk_gb: 10
																}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: gatk_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	String output_vcf = "~{prefix}.vcf.gz"
	Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVCleanPt4 \
			-V ~{vcf} \
			-O ~{output_vcf} \
			--revised-cn-list ~{rd_cn_revise} \
			~{if defined(outlier_samples_list) then "--outliers-list ~{outlier_samples_list}" else "" }
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}


task CleanVcf5 {
	input {
		File vcf
		String prefix
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size([vcf], "GB")
	RuntimeAttr runtime_default = object {
																	mem_gb: 2.0,
																	disk_gb: 50,
																	cpu_cores: 1,
																	preemptible_tries: 3,
																	max_retries: 1,
																	boot_disk_gb: 10
																}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: gatk_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	String output_vcf = "~{prefix}.vcf.gz"
	Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVCleanPt5 \
			-V ~{vcf} \
			-O ~{output_vcf}
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}


task RescueMobileElementDeletions {
	input {
		File vcf
		String prefix
		File LINE1
		File HERVK
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size(vcf, "GiB")
	RuntimeAttr runtime_default = object {
		mem_gb: 3.75 + input_size * 1.5,
		disk_gb: ceil(100.0 + input_size * 3.0),
		cpu_cores: 1,
		preemptible_tries: 3,
		max_retries: 1,
		boot_disk_gb: 10
	}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: sv_pipeline_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	command <<<
		set -euo pipefail

		python <<CODE
import os
import pysam
fin=pysam.VariantFile("~{vcf}")
fo=pysam.VariantFile("~{prefix}.bnd_del.vcf.gz", 'w', header = fin.header)
for record in fin:
		if record.info['SVTYPE'] in ['BND'] and record.info['STRANDS']=="+-" and record.chrom == record.info['CHR2'] and record.info['END2'] - record.start < 10000:
				record.info['SVLEN'] = record.info['END2'] - record.start
				fo.write(record)
fin.close()
fo.close()
CODE

		tabix -p vcf ~{prefix}.bnd_del.vcf.gz

		svtk vcf2bed ~{prefix}.bnd_del.vcf.gz -i ALL --include-filters ~{prefix}.bnd_del.bed
		bgzip ~{prefix}.bnd_del.bed

		bedtools coverage -wo -a ~{prefix}.bnd_del.bed.gz -b ~{LINE1} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_LINE1/' > manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv
		bedtools coverage -wo -a ~{prefix}.bnd_del.bed.gz -b ~{HERVK} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_HERVK/' >> manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv

		python <<CODE
import pysam
def SVID_MEI_DEL_readin(MEI_DEL_reset):
		out={}
		fin=open(MEI_DEL_reset)
		for line in fin:
				pin=line.strip().split()
				if not pin[0] in out.keys():
						out[pin[0]] = pin[3]
		fin.close()
		return out

hash_MEI_DEL_reset = SVID_MEI_DEL_readin("manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv")
fin=pysam.VariantFile("~{vcf}")
fo=pysam.VariantFile("~{prefix}.vcf.gz", 'w', header = fin.header)
for record in fin:
		if record.id in hash_MEI_DEL_reset.keys():
				del record.filter['UNRESOLVED']
				record.info['SVTYPE'] = 'DEL'
				record.info['SVLEN'] = record.info['END2'] - record.start
				record.stop = record.info['END2']
				record.info.pop("CHR2")
				record.info.pop("END2")
				record.info.pop("UNRESOLVED_TYPE")
				if hash_MEI_DEL_reset[record.id] == 'overlap_LINE1':
						record.alts = ('<DEL:ME:LINE1>',)
				if hash_MEI_DEL_reset[record.id] == 'overlap_HERVK':
						record.alts = ('<DEL:ME:HERVK>',)
		fo.write(record)
fin.close()
fo.close()
CODE
	>>>

	output {
		File out = "~{prefix}.vcf.gz"
	}
}


# Remove CNVs that are redundant with CPX events or other CNVs
task DropRedundantCnvs {
	input {
		File vcf
		String prefix
		String contig
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size(vcf, "GiB")
	# disk is cheap, read/write speed is proportional to disk size, and disk IO is a significant time factor:
	# in tests on large VCFs, memory usage is ~1.0 * input VCF size
	# the biggest disk usage is at the end of the task, with input + output VCF on disk
	Int cpu_cores = 2 # speed up compression / decompression of VCFs
	RuntimeAttr runtime_default = object {
		mem_gb: 3.75 + input_size * 1.5,
		disk_gb: ceil(100.0 + input_size * 2.0),
		cpu_cores: cpu_cores,
		preemptible_tries: 3,
		max_retries: 1,
		boot_disk_gb: 10
	}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: sv_pipeline_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	command <<<
		set -euo pipefail
		/opt/sv-pipeline/04_variant_resolution/scripts/resolve_cpx_cnv_redundancies.py \
			~{vcf} ~{prefix}.vcf.gz --temp-dir ./tmp
	>>>

	output {
		File out = "~{prefix}.vcf.gz"
	}
}


# Stitch fragmented RD-only calls found in 100% of the same samples
task StitchFragmentedCnvs {
	input {
		File vcf
		String prefix
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	Float input_size = size(vcf, "GB")
	RuntimeAttr runtime_default = object {
																	mem_gb: 7.5,
																	disk_gb: ceil(10.0 + input_size * 2),
																	cpu_cores: 1,
																	preemptible_tries: 3,
																	max_retries: 1,
																	boot_disk_gb: 10
																}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	Float mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
	Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

	runtime {
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: sv_pipeline_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	command <<<
		set -euo pipefail
		echo "First pass..."
		java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 ~{vcf} \
			| bgzip \
			> tmp.vcf.gz
		rm ~{vcf}
		echo "Second pass..."
		java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 tmp.vcf.gz \
			| bgzip \
			> ~{prefix}.vcf.gz
	>>>

	output {
		File stitched_vcf_shard = "~{prefix}.vcf.gz"
	}
}

# Add FILTER status for pockets of variants with high FP rate: wham-only DELs and Scramble-only SVAs with HIGH_SR_BACKGROUND
task AddHighFDRFilters {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + input_size * 3.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python <<CODE
import pysam
with pysam.VariantFile("~{vcf}", 'r') as fin:
  header = fin.header
  header.add_line("##FILTER=<ID=HIGH_ALGORITHM_FDR,Description=\"Categories of variants with low precision including Wham-only deletions and certain Scramble SVAs\">")
  with pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=header) as fo:
    for record in fin:
        if (record.info['ALGORITHMS'] == ('wham',) and record.info['SVTYPE'] == 'DEL') or \
          (record.info['ALGORITHMS'] == ('scramble',) and record.info['HIGH_SR_BACKGROUND'] and record.alts == ('<INS:ME:SVA>',)):
            record.filter.add('HIGH_ALGORITHM_FDR')
        fo.write(record)
CODE
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}



# Final VCF cleanup
task FinalCleanup {
	input {
		File vcf
		String contig
		String prefix
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	# generally assume working disk size is ~2 * inputs, and outputs are ~2 *inputs, and inputs are not removed
	# generally assume working memory is ~3 * inputs
	Float input_size = size(vcf, "GB")
	Float base_disk_gb = 10.0
	Float base_mem_gb = 2.0
	Float input_mem_scale = 3.0
	Float input_disk_scale = 5.0
	RuntimeAttr runtime_default = object {
		mem_gb: base_mem_gb + input_size * input_mem_scale,
		disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
		cpu_cores: 1,
		preemptible_tries: 3,
		max_retries: 1,
		boot_disk_gb: 10
	}
	RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
	runtime {
		memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
		disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
		cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
		preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
		maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
		docker: sv_pipeline_docker
		bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
	}

	command <<<
		set -eu -o pipefail
		
		/opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
			--chrom ~{contig} \
			--prefix ~{prefix} \
			~{vcf} stdout \
			| bcftools annotate --no-version -e 'SVTYPE=="CNV" && SVLEN<5000' -x INFO/MEMBERS -Oz -o ~{prefix}.vcf.gz
		tabix ~{prefix}.vcf.gz
	>>>

	output {
		File final_cleaned_shard = "~{prefix}.vcf.gz"
		File final_cleaned_shard_idx = "~{prefix}.vcf.gz.tbi"
	}
}