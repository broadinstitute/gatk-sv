version 1.0

workflow GTConcordanceWorkflow {
  input {
    File pav_vcf
    File kg_vcf
    File pg_vcf
    File pav_vcf_tbi
    File kg_vcf_tbi
    File pg_vcf_tbi
    Array[String] contigs
    File rscript  # GT_concordant_2.calcu_GT_concordance.R
    String sv_base_mini_docker
    String sv_pipeline_docker
  }


  scatter (contig in contigs) {

    call SplitVCFs as split_input_vcf {
      input:
        input_vcf = pav_vcf,
        input_vcf_tbi = pav_vcf_tbi,
        contig = contig,
        docker_image  = sv_base_mini_docker
    }


   call SplitVCFs as split_kg_vcf {
      input:
        input_vcf = kg_vcf,
        input_vcf_tbi = kg_vcf_tbi,
        contig = contig,
        docker_image  = sv_base_mini_docker
    }

   call SplitVCFs as split_pg_vcf {
      input:
        input_vcf = pg_vcf,
        input_vcf_tbi = pg_vcf_tbi,
        contig = contig,
        docker_image  = sv_base_mini_docker
    }

    call RunConcordance {
      input:
        pav_contig = split_input_vcf.contig_vcf,
        kg_contig = split_kg_vcf.contig_vcf,
        pg_contig = split_pg_vcf.contig_vcf,
        rscript   = rscript,
        docker_image = sv_pipeline_docker
    }
  }

  output {
    Array[File] kg_outputs = RunConcordance.kg_output
    Array[File] pg_outputs = RunConcordance.pg_output
  }
}

task SplitVCFs {
  input {
    File input_vcf
    File input_vcf_tbi
    String contig
    String docker_image
  }

  String prefix = basename(input_vcf, ".vcf.gz") 
  
  command <<<
    bcftools view -r ~{contig} -Oz -o ~{prefix}.~{contig}.vcf.gz ~{input_vcf}
   >>>

  output {
    File contig_vcf = "~{prefix}.~{contig}.vcf.gz"

  }

  runtime {
    docker: docker_image
    memory: "4G"
    cpu: 1
  }
}


task RunConcordance {
  input {
    File pav_contig
    File kg_contig
    File pg_contig
    File rscript
    String docker_image
  }

  String kg_output = sub(basename(kg_contig), ".vcf.gz", ".SVID_concor")
  String pg_output = sub(basename(pg_contig), ".vcf.gz", ".SVID_concor")

  command <<<
    Rscript ~{rscript} \
      --pav_input ~{pav_contig} \
      --KG_input ~{kg_contig} \
      --PG_input ~{pg_contig} \
      --KG_output ~{kg_output} \
      --PG_output ~{pg_output}
  >>>

  output {
    File kg_output = "~{kg_output}"
    File pg_output = "~{pg_output}"
  }

  runtime {
    docker: docker_image
    memory: "8G"
    cpu: 2
  }
}
