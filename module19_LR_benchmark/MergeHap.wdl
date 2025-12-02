version 1.0

import "Structs.wdl"

workflow MergeHap {
	input {
		Array[File] vcf_pat_vcfs
		Array[File] vcf_pat_idxes
		Array[File] vcf_mat_vcfs
		Array[File] vcf_mat_idxes
		Arrayp[String] samples
		String collapse_params
		File ref_fa
		File ref_fai
		String truvari_docker
		RuntimeAttr? runtime_attr_override

	}

	scatter(i in range(length(vcf_pat_vcfs))){
		call TruvariCollapse{
			input:
				vcf_pat = vcf_pat_vcfs[i],
				vcf_pat_idx = vcf_pat_idxes[i],
				vcf_mat = vcf_mat_vcfs[i],
				vcf_mat_idx = vcf_mat_idxes[i],
				sample = samples[i],
				collapse_params = collapse_params,
				ref_fa = ref_fa,
				ref_fai = ref_fai,
				docker = truvari_docker,
				runtime_attr_override = runtime_attr_override
		}
	}

	output{
		Array[File] merged_vcf = TruvariCollapse.diploid_vcf
		Array[File] merged_idx = TruvariCollapse.diploid_vcf_idx
	}
}


task TruvariCollapse {
	input {
		File vcf_pat
		File vcf_pat_idx
		File vcf_mat
		File vcf_mat_idx
		String sample
		String collapse_params
		File ref_fa
		File ref_fai
		String docker
		RuntimeAttr? runtime_attr_override
	}

	command <<<
		set -euo pipefail

		bcftools concat \
			-a \
			~{vcf_pat} \
			~{vcf_mat} \
		| bcftools sort \
			-Oz -o combined.vcf.gz
		
		tabix combined.vcf.gz

		truvari collapse \
			--reference ~{ref_fa} \
			-i combined.vcf.gz \
			-o ~{sample}.merged.vcf.gz \
			-c ~{sample}.collapsed.vcf.gz \
			--hap \
			~{collapse_params}

		bcftools sort \
			~{sample}.merged.vcf.gz \
			-Oz -o ~{sample}.merged.sorted.vcf.gz
		
		tabix ~{sample}.merged.sorted.vcf.gz
	>>>

	output {
		File diploid_vcf = "~{sample}.merged.sorted.vcf.gz"
		File diploid_vcf_idx = "~{sample}.merged.sorted.vcf.gz.tbi"
	}

	RuntimeAttr default_attr = object {
		cpu_cores: 1,
		mem_gb: 4,
		disk_gb: ceil(size(vcf_pat, "GB") + size(vcf_mat, "GB")) * 3 + 10,
		boot_disk_gb: 10,
		preemptible_tries: 1,
		max_retries: 0
	}
	RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
	runtime {
		cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
		memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
		disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
		bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
		docker: docker
		preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
		maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
	}
}
