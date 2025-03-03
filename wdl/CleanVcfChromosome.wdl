version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "FormatVcfForGatk.wdl" as fvcf
import "HailMerge.wdl" as HailMerge

workflow CleanVcfChromosome {
	input {
		File vcf
		String contig
		String chr_x
		String chr_y
		String prefix

		File background_list
		File bothsides_pass_list
		File? outlier_samples_list
		File ped_file
		File ploidy_table
		File allosome_fai
		
		File HERVK_reference
		File LINE1_reference

		Boolean use_hail
		String? gcs_project

		String gatk_docker
		String linux_docker
		String sv_base_mini_docker
		String sv_pipeline_docker
    File? svtk_to_gatk_script  # For debugging
    File? make_clean_gq_script

		# overrides for local tasks
		RuntimeAttr? runtime_attr_preprocess
		RuntimeAttr? runtime_attr_revise_overlapping_cnvs
		RuntimeAttr? runtime_attr_revise_large_cnvs
		RuntimeAttr? runtime_attr_revise_abnormal_allosomes
		RuntimeAttr? runtime_attr_revise_multiallelics
		RuntimeAttr? runtime_attr_postprocess
		RuntimeAttr? runtime_override_stitch_fragmented_cnvs
		RuntimeAttr? runtime_override_final_cleanup
		RuntimeAttr? runtime_override_rescue_me_dels
		RuntimeAttr? runtime_attr_add_high_fp_rate_filters

		RuntimeAttr? runtime_override_preconcat_drc
		RuntimeAttr? runtime_override_hail_merge_drc
		RuntimeAttr? runtime_override_fix_header_drc

		RuntimeAttr? runtime_override_drop_redundant_cnvs
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

	call CleanVcfPreprocess {
		input:
			vcf=FormatVcfToClean.out,
			chr_x=chr_x,
			chr_y=chr_y,
			background_list=background_list,
			bothsides_pass_list=bothsides_pass_list,
			prefix="~{prefix}.preprocess",
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_attr_preprocess
	}

	call CleanVcfReviseOverlappingCnvs {
		input:
			vcf=CleanVcfPreprocess.out,
			prefix="~{prefix}.revise_overlapping_cnvs",
			gatk_docker=gatk_docker,
			runtime_attr_override=runtime_attr_revise_overlapping_cnvs
	}

	call CleanVcfReviseLargeCnvs {
		input:
			vcf=CleanVcfReviseOverlappingCnvs.out,
			outlier_samples_list=outlier_samples_list,
			prefix="~{prefix}.revise_large_cnvs",
			gatk_docker=gatk_docker,
			runtime_attr_override=runtime_attr_revise_large_cnvs
	}

	call CleanVcfReviseAbnormalAllosomes {
		input:
			vcf=CleanVcfReviseLargeCnvs.out,
			prefix="~{prefix}.revise_abnormal_allosomes",
			gatk_docker=gatk_docker,
			runtime_attr_override=runtime_attr_revise_abnormal_allosomes
	}

	call CleanVcfReviseMultiallelicCnvs {
		input:
			vcf=CleanVcfReviseAbnormalAllosomes.out,
			prefix="~{prefix}.revise_multiallelic_cnvs",
			gatk_docker=gatk_docker,
			runtime_attr_override=runtime_attr_revise_multiallelics
	}

	call CleanVcfPostprocess {
		input:
			vcf=CleanVcfReviseMultiallelicCnvs.out,
			prefix="~{prefix}.postprocess",
			sv_pipeline_docker=sv_pipeline_docker,
			runtime_attr_override=runtime_attr_postprocess
	}

	call DropRedundantCnvs {
		input:
			vcf=CleanVcfPostprocess.out,
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
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_format
  }
	
	output {
		File out = FormatVcfToOutput.out
		File out_idx = FormatVcfToOutput.out_index
	}
}

task CleanVcfPreprocess {
	input {
		File vcf
		String chr_x
		String chr_y
		File background_list
		File bothsides_pass_list
		String prefix
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr runtime_default = object {
		mem_gb: 3.75,
		disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
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

	Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
	String output_vcf = "~{prefix}.vcf.gz"

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi

		python /opt/sv-pipeline/04_variant_resolution/scripts/replace_ev_numeric_code_with_string.py \
			~{vcf} \
			processed.vcf.gz

		zgrep '^##' processed.vcf.gz > header.txt

		cat <<EOF >> header.txt
		##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved">
		##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description="Variant has high number of SR splits in background samples">
		##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description="Variant has read-level support for both sides of breakpoint">
		##INFO=<ID=REVISED_EVENT,Number=0,Type=Flag,Description="Variant has been revised due to a copy number mismatch">
		EOF

		zgrep '^#CHROM' processed.vcf.gz >> header.txt

		bcftools view processed.vcf.gz | bcftools reheader -h header.txt | bgzip -c > processed.reheader.vcf.gz

		rm processed.vcf.gz header.txt
		
		python /opt/sv-pipeline/04_variant_resolution/scripts/cleanvcf_preprocess.py \
			-V processed.reheader.vcf.gz \
			-O ~{output_vcf} \
			--chrX ~{chr_x} \
			--chrY ~{chr_y} \
			--fail-list ~{background_list} \
			--pass-list ~{bothsides_pass_list}

		tabix -p vcf ~{output_vcf}
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}

task CleanVcfReviseOverlappingCnvs {
	input {
		File vcf
		String prefix
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr runtime_default = object {
		mem_gb: 3.75,
		disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
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

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVReviseOverlappingCnvs \
			-V ~{vcf} \
			-O ~{output_vcf}
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}

task CleanVcfReviseMultiallelicCnvs {
	input {
		File vcf
		File? outlier_samples_list
		String prefix
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr runtime_default = object {
		mem_gb: 3.75,
		disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
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

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVReviseMutliallelicCnvs \
			-V ~{vcf} \
			-O ~{output_vcf} \
			~{if defined(outlier_samples_list) then "--outlier-samples ~{outlier_samples_list}" else "" }
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}

task CleanVcfReviseAbnormalAllosomes {
	input {
		File vcf
		String prefix
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr runtime_default = object {
		mem_gb: 3.75,
		disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
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

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVReviseAbnormalAllosomes \
			-V ~{vcf} \
			-O ~{output_vcf}
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}

task CleanVcfReviseOverlappingMultiallelics {
	input {
		File vcf
		String prefix
		String gatk_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr runtime_default = object {
		mem_gb: 3.75,
		disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
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

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi
		
		gatk --java-options "-Xmx~{java_mem_mb}m" SVReviseOverlappingMultiallelics \
			-V ~{vcf} \
			-O ~{output_vcf}
	>>>

	output {
		File out="~{output_vcf}"
		File out_idx="~{output_vcf}.tbi"
	}
}

task CleanVcfPostprocess {
	input {
		File vcf
		String prefix
		String sv_pipeline_docker
		RuntimeAttr? runtime_attr_override
	}

	RuntimeAttr runtime_default = object {
		mem_gb: 3.75,
		disk_gb: ceil(10.0 + size(vcf, "GB") * 2),
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

	Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
	String output_vcf = "~{prefix}.vcf.gz"

	command <<<
		set -euo pipefail

		if [ ! -f "~{vcf}.tbi" ]; then
			tabix -p vcf ~{vcf}
		fi

		python /opt/sv-pipeline/04_variant_resolution/scripts/cleanvcf_postprocess.py \
			-V ~{vcf} \
			-O processed.vcf.gz

		bcftools annotate -x INFO/MULTIALLELIC,INFO/UNRESOLVED,INFO/EVENT,INFO/REVISED_EVENT,INFO/MULTI_CNV,INFO/varGQ processed.vcf.gz -o processed.annotated.vcf.gz -O z

		bcftools view -h processed.annotated.vcf.gz | grep "^##" | \
			grep -v -E "CIPOS|CIEND|RMSSTD|source|bcftools|GATKCommandLine|##FORMAT=<ID=EV>|##ALT=<ID=UNR>|##INFO=<ID=(MULTIALLELIC|UNRESOLVED|EVENT|REVISED_EVENT|MULTI_CNV|varGQ)" > temp_header.txt
		echo '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">' >> temp_header.txt
		echo '##ALT=<ID=CNV,Description="Copy Number Polymorphism">' >> temp_header.txt

		bcftools view -h processed.annotated.vcf.gz | grep "^#CHROM" > chrom_header.txt

		cat temp_header.txt chrom_header.txt > header.txt
		
		bcftools reheader -h header.txt processed.annotated.vcf.gz -o ~{output_vcf}
		
		tabix -p vcf ~{output_vcf}
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