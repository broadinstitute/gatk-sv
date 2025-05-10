version 1.0

workflow PreprocessVcfForVapor {
	input {
		String sample_id            # Sample identifier
		File vcf_path             	# Path to the input VCF file
		String caller				# Caller to standardize VCF for

		File contigs_fai          	# Path to the contigs file
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

	call Vcf2Bed {
		input:
			sample_id = sample_id,
			vcf_path = StandardizeVcf.standardized_vcf,
			sv_pipeline_docker = sv_pipeline_docker
	}

	output {
		File dragen_sr_bed = Vcf2Bed.vcf2bed_vapor
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
			~{sample_id}.std_dragen.vcf.gz \
			~{caller}
	>>>

	output {
		File standardized_vcf = "~{sample_id}.std_dragen.vcf.gz"
	}

	runtime {
		cpu: 1
		memory: "2 GiB"
		disks: "local-disk 2 HDD"
		docker: sv_pipeline_docker
	}
}

task Vcf2Bed {
	input {
		String sample_id
		File vcf_path
		String sv_pipeline_docker
	}

	command <<<
		set -eu -o pipefail

		svtk vcf2bed --info SVTYPE --info SVLEN ~{vcf_path} - | awk '$7 != "BND"' > ~{sample_id}.bed
	>>>

	output {
		File vcf2bed_vapor = "~{sample_id}.bed"
	}

	runtime {
		cpu: 1
		memory: "2 GiB"
		disks: "local-disk 2 HDD"
		docker: sv_pipeline_docker
	}
}
