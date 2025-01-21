version 1.0

workflow CombineVcfsForMakeGq {
	input {
		Array[String] sample_ids            # Array of sample identifiers
		Array[File] vcf_paths               # Array of input VCF file paths (same order as sample_ids)

		File ploidy_table                   # Path to the ploidy table file
		File contigs_fai                    # Path to the contigs FASTA index file
		File ref_fasta                      # Path to the reference FASTA file
		File ref_fasta_fai                  # Path to the reference FASTA index file
		File ref_dict                       # Path to the reference dictionary
		Int min_size               					# Minimum size for standardization

		String sv_pipeline_docker           # Docker image path for GATK-SV
		String gatk_docker                  # Docker image path for GATK
	}

	scatter (idx in range(length(sample_ids))) {
		call StandardizeVcf {
			input:
				sample_id 						= sample_ids[idx],
				vcf_path       				= vcf_paths[idx],
				contigs_fai   				= contigs_fai,
				min_size      				= min_size,
				sv_pipeline_docker  	= sv_pipeline_docker
		}

		call SortVcf {
			input:
				sample_id          = sample_ids[idx],
				vcf_path           = StandardizeVcf.standardized_vcf,
				sv_pipeline_docker = sv_pipeline_docker
		}

		call FormatVcfForGatk {
			input:
				sample_id       = sample_ids[idx],
				vcf_path        = SortVcf.sorted_vcf,
				ploidy_table    = ploidy_table,
				sv_pipeline_docker  = sv_pipeline_docker
		}
	}

	Array[File] formatted_vcfs = FormatVcfForGatk.formatted_vcf
	Array[File] formatted_vcf_indices = FormatVcfForGatk.formatted_vcf_index
	
	call SVCluster {
		input:
			formatted_vcfs  			= formatted_vcfs,
			formatted_vcf_indices = formatted_vcf_indices,
			ploidy_table   				= ploidy_table,
			ref_fasta      				= ref_fasta,
			ref_fasta_fai  				= ref_fasta_fai,
			ref_dict       				= ref_dict,
			variant_prefix				= "clustered",
			gatk_docker    				= gatk_docker
	}

	output {
		File clustered_vcf = SVCluster.clustered_vcf
		File clustered_vcf_index = SVCluster.clustered_vcf_index
	}
}

task StandardizeVcf {
	input {
		String sample_id
		File vcf_path
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
			dragen
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

task SortVcf {
	input {
		String sample_id
		File vcf_path
		String sv_pipeline_docker
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
			--out ~{sample_id}.formatted.vcf.gz \
			--ploidy-table ~{ploidy_table} \
			--fix-end

		tabix -p vcf ~{sample_id}.formatted.vcf.gz
	>>>

	output {
		File formatted_vcf = "~{sample_id}.formatted.vcf.gz"
		File formatted_vcf_index = "~{sample_id}.formatted.vcf.gz.tbi"
	}

	runtime {
		cpu: 1
		memory: "4 GiB"
		disks: "local-disk 5 SSD"
		docker: sv_pipeline_docker
	}
}

task SVCluster {
	input {
		Array[File] formatted_vcfs
		Array[File] formatted_vcf_indices
		File ploidy_table
		File ref_fasta
		File ref_fasta_fai
		File ref_dict
		String variant_prefix
		String gatk_docker
	}

	command <<<
		set -eu -o pipefail

		awk '{print "-V "$0}' ~{write_lines(formatted_vcfs)} > arguments.txt

		gatk SVCluster \
			--arguments_file arguments.txt \
			--output clustered.vcf.gz \
			--ploidy-table ~{ploidy_table} \
			--reference ~{ref_fasta} \
			--variant-prefix ~{variant_prefix} \
			--depth-interval-overlap 1 --depth-breakend-window 0 \
			--mixed-interval-overlap 1 --mixed-breakend-window 0 \
			--pesr-interval-overlap 1 --pesr-breakend-window 0 \
	>>>

	output {
		File clustered_vcf = "clustered.vcf.gz"
		File clustered_vcf_index = "clustered.vcf.gz.tbi"
	}

	runtime {
		cpu: 4
		memory: "16 GiB"
		disks: "local-disk 20 SSD"
		docker: gatk_docker
	}
}
