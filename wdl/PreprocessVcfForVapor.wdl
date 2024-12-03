version 1.0

workflow PreprocessVcfForVapor {
	input {
		String sample_id            # Sample identifier
		String vcf_path             # Path to the input VCF file
		String sv_pipeline_docker   # Docker image path for GATK-SV
		String svtk_docker          # Docker image path for svtk
		File contigs_file          	# Path to the contigs file
		Int min_size               	# Minimum size for standardization
	}

	call StandardizeVcf {
		input:
			sample_id = sample_id,
			vcf_path = vcf_path,
			contigs_file = contigs_file,
			min_size = min_size,
			svtk_docker = svtk_docker
	}

	call SortVcf {
		input:
			sample_id = sample_id,
			vcf_path = StandardizeVcf.standardized_vcf,
			sv_pipeline_docker = sv_pipeline_docker
	}

	call Vcf2Bed {
		input:
			sample_id = sample_id,
			vcf_path = SortVcf.sorted_vcf,
			sv_pipeline_docker = sv_pipeline_docker
	}

	output {
		File vcf2bed_vapor = Vcf2Bed.vcf2bed_vapor
	}
}

task StandardizeVcf {
	input {
		String sample_id            # Sample identifier
		String vcf_path             # Path to the input VCF file
		File contigs_file           # Path to the contigs file
		Int min_size                # Minimum size for standardization
		String svtk_docker          # Docker image path for svtk
	}

	command <<<
		set -eu -o pipefail

		svtk standardize \
			--sample-names ~{sample_id} \
			--contigs ~{contigs_file} \
			--min-size ~{min_size} \
			~{vcf_path} \
			~{sample_id}.std_dragen.vcf.gz \
			dragen
	>>>

	output {
		File standardized_vcf = "~{sample_id}.std_dragen.vcf.gz"
	}

	runtime {
		cpu: 1
		memory: "2 GiB"
		disks: "local-disk 2 HDD"
		docker: svtk_docker
	}
}

task SortVcf {
	input {
		String sample_id            # Sample identifier
		File vcf_path               # Path to the standardized VCF file
		String sv_pipeline_docker   # Docker image path for GATK-SV
	}

	command <<<
		set -eu -o pipefail

		bcftools sort -O z -o ~{sample_id}.std_dragen_sorted.vcf.gz ~{vcf_path}
	>>>

	output {
		File sorted_vcf = "~{sample_id}.std_dragen_sorted.vcf.gz"
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
		String sample_id            # Sample identifier
		File vcf_path             	# Path to the sorted VCF file
		String sv_pipeline_docker   # Docker image path for GATK-SV
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
