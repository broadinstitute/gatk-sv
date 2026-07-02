version 1.0

workflow SVShellPostProcessIntermediates {
	input {
		String sample_id            # Sample identifier

		File contigs_fai          	# Path to the contigs file
		File ploidy_table		  	# Path to the ploidy table file
		Int min_size               	# Minimum size for standardization

		File TEMP_intermediates     # Path to the tar.gz of intermediates

		File? dragen_sv             # Optional Dragen SV VCF
		File? dragen_cnv            # Optional Dragen CNV VCF

		String sv_pipeline_docker   # Docker image path for GATK-SV
	}

	call ExtractIntermediates {
		input:
			sample_id = sample_id,
			intermediates_tar = TEMP_intermediates,
			sv_pipeline_docker = sv_pipeline_docker
	}

	if (defined(dragen_sv)) {
		call StandardizeVcf as StandardizeDragenSv {
			input:
				sample_id = sample_id,
				vcf_path = select_first([dragen_sv]),
				caller = "dragen",
				contigs_fai = contigs_fai,
				min_size = min_size,
				sv_pipeline_docker = sv_pipeline_docker
		}

		call FormatVcfForGatk as FormatDragenSv {
			input:
				sample_id = sample_id,
				vcf_path = StandardizeDragenSv.standardized_vcf,
				ploidy_table = ploidy_table,
				sv_pipeline_docker = sv_pipeline_docker
		}
	}

	if (defined(dragen_cnv)) {
		call StandardizeVcf as StandardizeDragenCnv {
			input:
				sample_id = sample_id,
				vcf_path = select_first([dragen_cnv]),
				caller = "dragen",
				contigs_fai = contigs_fai,
				min_size = min_size,
				sv_pipeline_docker = sv_pipeline_docker
		}

		call FormatVcfForGatk as FormatDragenCnv {
			input:
				sample_id = sample_id,
				vcf_path = StandardizeDragenCnv.standardized_vcf,
				ploidy_table = ploidy_table,
				sv_pipeline_docker = sv_pipeline_docker
		}
	}

	if (ExtractIntermediates.has_scramble_vcf) {
		call StandardizeVcf as StandardizeScramble {
			input:
				sample_id = sample_id,
				vcf_path = ExtractIntermediates.scramble_vcf,
				caller = "scramble",
				contigs_fai = contigs_fai,
				min_size = min_size,
				sv_pipeline_docker = sv_pipeline_docker
		}

		call FormatVcfForGatk as FormatScramble {
			input:
				sample_id = sample_id,
				vcf_path = StandardizeScramble.standardized_vcf,
				ploidy_table = ploidy_table,
				sv_pipeline_docker = sv_pipeline_docker
		}
	}

	if (ExtractIntermediates.has_wham_vcf) {
		call StandardizeVcf as StandardizeWham {
			input:
				sample_id = sample_id,
				vcf_path = ExtractIntermediates.wham_vcf,
				caller = "wham",
				contigs_fai = contigs_fai,
				min_size = min_size,
				sv_pipeline_docker = sv_pipeline_docker
		}

		call FormatVcfForGatk as FormatWham {
			input:
				sample_id = sample_id,
				vcf_path = StandardizeWham.standardized_vcf,
				ploidy_table = ploidy_table,
				sv_pipeline_docker = sv_pipeline_docker
		}
	}

	output {
		# Raw files extracted from the intermediates tar.
		File scramble_tsv = ExtractIntermediates.scramble_tsv
		File scramble_vcf = ExtractIntermediates.scramble_vcf
		File scramble_vcf_index = ExtractIntermediates.scramble_vcf_index
		File wham_vcf = ExtractIntermediates.wham_vcf
		File wham_vcf_index = ExtractIntermediates.wham_vcf_index

		# Formatted VCFs (only produced when the corresponding raw VCF was present).
		File? dragen_sv_vcf_formatted = FormatDragenSv.formatted_vcf
		File? dragen_sv_vcf_formatted_index = FormatDragenSv.formatted_vcf_index
		File? dragen_cnv_vcf_formatted = FormatDragenCnv.formatted_vcf
		File? dragen_cnv_vcf_formatted_index = FormatDragenCnv.formatted_vcf_index
		File? scramble_vcf_formatted = FormatScramble.formatted_vcf
		File? scramble_vcf_formatted_index = FormatScramble.formatted_vcf_index
		File? wham_vcf_formatted = FormatWham.formatted_vcf
		File? wham_vcf_formatted_index = FormatWham.formatted_vcf_index
	}
}

task StandardizeVcf {
	input {
		String sample_id
		File vcf_path
        String caller
		File contigs_fai
		Int min_size
		String sv_pipeline_docker
	}

	command <<<
		set -eu -o pipefail
		
		svtk standardize \
			--sample-names ~{sample_id} \
			--contigs ~{contigs_fai} \
			--min-size ~{min_size} \
			~{vcf_path} \
			~{sample_id}.std.vcf.gz \
			~{caller}

		tabix -p vcf ~{sample_id}.std.vcf.gz
	>>>

	output {
		File standardized_vcf = "~{sample_id}.std.vcf.gz"
		File standardized_vcf_index = "~{sample_id}.std.vcf.gz.tbi"
	}

	runtime {
		cpu: 1
		memory: "2 GiB"
		disks: "local-disk 2 HDD"
		docker: sv_pipeline_docker
	}
}

task FormatVcfForGatk {
	input {
		String sample_id
		File vcf_path
		File ploidy_table
		String sv_pipeline_docker
	}

	command <<<
		set -eu -o pipefail

		python /opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py \
			--vcf ~{vcf_path} \
			--out ~{sample_id}.fmt.vcf.gz \
			--ploidy-table ~{ploidy_table} \
			--fix-end

		tabix -p vcf ~{sample_id}.fmt.vcf.gz
	>>>

	output {
		File formatted_vcf = "~{sample_id}.fmt.vcf.gz"
		File formatted_vcf_index = "~{sample_id}.fmt.vcf.gz.tbi"
	}

	runtime {
		cpu: 1
		memory: "4 GiB"
		disks: "local-disk 5 SSD"
		docker: sv_pipeline_docker
	}
}

task ExtractIntermediates {
	input {
		String sample_id
		File intermediates_tar
		String sv_pipeline_docker
	}

	command <<<
		set -eu -o pipefail

		mkdir -p extracted
		tar -xzf ~{intermediates_tar} -C extracted

		# Locate the output_GatherSampleEvidence_* directory (name suffix varies per run).
		run_dir=$(find extracted -type d -name 'output_GatherSampleEvidence_*' | head -n 1 || true)

		# Copy each expected output if present, otherwise create an empty placeholder
		# so the WDL 1.0 File outputs always resolve. A companion ".present" file
		# records whether the real artifact was found, so the workflow can skip
		# downstream tasks for callers whose VCFs are missing.
		for fname in \
			"~{sample_id}.scramble.tsv.gz" \
			"~{sample_id}.scramble.vcf.gz" \
			"~{sample_id}.scramble.vcf.gz.tbi" \
			"~{sample_id}.wham.vcf.gz" \
			"~{sample_id}.wham.vcf.gz.tbi"; do
			src=""
			if [ -n "$run_dir" ] && [ -f "$run_dir/$fname" ]; then
				src="$run_dir/$fname"
			fi
			if [ -z "$src" ]; then
				# Fall back to a broader search in case the layout differs slightly.
				src=$(find extracted -type f -name "$fname" | head -n 1 || true)
			fi
			if [ -n "$src" ] && [ -f "$src" ]; then
				cp "$src" "$fname"
				echo "true" > "$fname.present"
			else
				: > "$fname"
				echo "false" > "$fname.present"
			fi
		done
	>>>

	output {
		File scramble_tsv       = "~{sample_id}.scramble.tsv.gz"
		File scramble_vcf       = "~{sample_id}.scramble.vcf.gz"
		File scramble_vcf_index = "~{sample_id}.scramble.vcf.gz.tbi"
		File wham_vcf           = "~{sample_id}.wham.vcf.gz"
		File wham_vcf_index     = "~{sample_id}.wham.vcf.gz.tbi"

		Boolean has_scramble_vcf = read_boolean("~{sample_id}.scramble.vcf.gz.present")
		Boolean has_wham_vcf     = read_boolean("~{sample_id}.wham.vcf.gz.present")
	}

	runtime {
		cpu: 1
		memory: "2 GiB"
		disks: "local-disk 20 HDD"
		docker: sv_pipeline_docker
	}
}