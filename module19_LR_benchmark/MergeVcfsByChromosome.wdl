version 1.0

workflow MergeVcfsByChromosome {
  input {
    Array[File] input_vcfs            # bgzipped VCFs with .tbi
    String chrom         # e.g. ["1", "2", ..., "22", "X"]
  }

    scatter (vcf in input_vcfs) {
      call ExtractChromosomeVcf {
        input:
          input_vcf = vcf,
          chromosome = chrom
      }
    }

    call MergeVcfs {
      input:
        input_vcfs = ExtractChromosomeVcf.output_vcf,
        output_name = "${chrom}.vcf.gz"
    }

  output {
    File merged_vcf = MergeVcfs.output_merged_vcf
    File merged_vcf_index = MergeVcfs.output_merged_vcf_index
  }
}

# Task 1: Extract a chromosome from a VCF
task ExtractChromosomeVcf {
  input {
    File input_vcf
    String chromosome
  }

  command <<<
    set -e
    bcftools view -r ~{chromosome} ~{input_vcf} -Oz -o chr~{chromosome}.vcf.gz
    tabix -p vcf chr~{chromosome}.vcf.gz
  >>>

  output {
    File output_vcf = "chr~{chromosome}.vcf.gz"
  }

  runtime {
    docker: "biocontainers/bcftools:v1.17-1-deb-py3"
    cpu: 1
    memory: "2G"
  }
}

# Task 2: Merge multiple VCFs
task MergeVcfs {
  input {
    Array[File] input_vcfs
    String output_name
  }

  command <<<
    set -e
    bcftools merge ~{sep=' ' input_vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}
  >>>

  output {
    File output_merged_vcf = output_name
    File output_merged_vcf_index = "${output_name}.tbi"
  }

  runtime {
    docker: "biocontainers/bcftools:v1.17-1-deb-py3"
    cpu: 2
    memory: "4G"
  }
}

# Task 3: Concatenate per-chromosome VCFs
task ConcatVcfs {
  input {
    Array[File] input_vcfs
    String output_name
  }

  command <<<
    set -e
    bcftools concat ~{sep=' ' input_vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}
  >>>

  output {
    File output_vcf = output_name
    File output_vcf_index = "${output_name}.tbi"
  }

  runtime {
    docker: "biocontainers/bcftools:v1.17-1-deb-py3"
    cpu: 2
    memory: "4G"
  }
}