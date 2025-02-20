version 1.0

workflow SVConcordanceSimple {
  input {
    File eval_vcf
    File truth_vcf
    String output_prefix

    Float depth_interval_overlap

    File reference_dict

    String gatk_docker
  }

  call SVConcordanceTask {
    input:
      eval_vcf=eval_vcf,
      truth_vcf=truth_vcf,
      output_prefix=output_prefix,
      reference_dict=reference_dict,
      additional_args="--depth-interval-overlap ~{depth_interval_overlap}",
      gatk_docker=gatk_docker
  }

  output {
    File conc_vcf = SVConcordanceTask.out
    File conc_vcf_idx = SVConcordanceTask.out_idx
  }
}

task SVConcordanceTask {
  input {
    File truth_vcf
    File eval_vcf
    String output_prefix
    
    File reference_dict

    String? additional_args
    String gatk_docker
  }

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_idx = "~{output_prefix}.vcf.gz.tbi"
  }

  command <<<
    set -euo pipefail

    gatk SVConcordance \
      --eval ~{eval_vcf} \
      --truth ~{truth_vcf} \
      --sequence-dictionary ~{reference_dict} \
      ~{additional_args} \
      -O ~{output_prefix}.vcf.gz
      
  >>>

  runtime {
		cpu: 1
		memory: "2 GiB"
		disks: "local-disk 2 HDD"
		docker: gatk_docker
	}
}