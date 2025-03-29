version 1.0

workflow CombineVcfsForMakeGq {
	input {
		Array[File] vcf_paths               # Array of input VCF file paths (same order as sample_ids)
		Array[File] vcf_paths_idx           # Array of input VCF index file paths (same order as sample_ids)

		File ploidy_table                   # Path to the ploidy table file
		File ref_fasta                      # Path to the reference FASTA file
		File ref_fasta_fai                  # Path to the reference FASTA index file
		File ref_dict                       # Path to the reference dictionary

		String gatk_docker                  # Docker image path for GATK
	}
	
	call SVCluster {
		input:
			formatted_vcfs  			= vcf_paths,
			formatted_vcf_indices = vcf_paths_idx,
			ploidy_table   				= ploidy_table,
			ref_fasta      				= ref_fasta,
			ref_fasta_fai  				= ref_fasta_fai,
			ref_dict       				= ref_dict,
			gatk_docker    				= gatk_docker
	}

	output {
		File clustered_vcf = SVCluster.clustered_vcf
		File clustered_vcf_index = SVCluster.clustered_vcf_index
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
			--pesr-interval-overlap 1 --pesr-breakend-window 0 \
			--depth-interval-overlap 1 --depth-breakend-window 0 \
			--mixed-interval-overlap 1 --mixed-breakend-window 0
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