version 1.0

workflow PreprocessVcfForMakeGq {
	input {
		String sample_id            # Sample identifier
		File vcf_path             	# Path to the input VCF file
        String caller               # Caller to standardize VCF for

		File contigs_fai          	# Path to the contigs file
		File ploidy_table		  	# Path to the ploidy table file
		Int min_size               	# Minimum size for standardization
		
		String sv_pipeline_docker   # Docker image path for GATK-SV
	}

	call StandardizeVcf {
		input:
			sample_id = sample_id,
			vcf_path = vcf_path,
            caller = caller,
			contigs_fai = contigs_fai,
			min_size = min_size,
			sv_pipeline_docker = sv_pipeline_docker
	}

	call FormatVcfForGatk {
		input:
			sample_id = sample_id,
			vcf_path = StandardizeVcf.standardized_vcf,
			ploidy_table = ploidy_table,
			sv_pipeline_docker  = sv_pipeline_docker
	}

	output {
		File dragen_vcf_std = FormatVcfForGatk.formatted_vcf
		File dragen_vcf_std_idx = FormatVcfForGatk.formatted_vcf_index
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