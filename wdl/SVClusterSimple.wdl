version 1.0

workflow SVClusterSimple {
	input {
		File vcf               # Input VCF
		File vcf_idx           # Input VCF index file
		Float defrag_fration	 # Defragment padding fraction

		File ploidy_table      # Path to the ploidy table file
		File ref_fasta         # Path to the reference FASTA file
		File ref_fasta_fai     # Path to the reference FASTA index file
		File ref_dict          # Path to the reference dictionary

		String gatk_docker     # Docker image path for GATK
	}
	
	call SVCluster {
		input:
			vcf  						= vcf,
			vcf_idx 				= vcf_idx,
			defrag_fration	= defrag_fration,
			ploidy_table  	= ploidy_table,
			ref_fasta     	= ref_fasta,
			ref_fasta_fai 	= ref_fasta_fai,
			ref_dict      	= ref_dict,
			gatk_docker   	= gatk_docker
	}

	output {
		File defrag_vcf	 		= SVCluster.defrag_vcf
		File defrag_vcf_idx = SVCluster.defrag_vcf_idx
	}
}

task SVCluster {
	input {
		File vcf
		File vcf_idx
		Float defrag_fration
		File ploidy_table
		File ref_fasta
		File ref_fasta_fai
		File ref_dict
		String gatk_docker
	}

	command <<<
		set -eu -o pipefail

		gatk SVCluster \
			-V ~{vcf} \
			--output defragmented.vcf.gz \
			--ploidy-table ~{ploidy_table} \
			--reference ~{ref_fasta} \
			--algorithm DEFRAGMENT_CNV \
			--defrag-padding-fraction ~{defrag_fration} \
			--defrag-sample-overlap 1.0
	>>>

	output {
		File defrag_vcf = "defragmented.vcf.gz"
		File defrag_vcf_idx = "defragmented.vcf.gz.tbi"
	}

	runtime {
		cpu: 1
		memory: "4 GiB"
		disks: "local-disk 4 HDD"
		docker: gatk_docker
	}
}