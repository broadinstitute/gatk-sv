version 1.0

## Given a VCF (+ index) and a list of sample IDs (one per line), extract
## each sample's informative variants -- sites where that sample's GT is
## neither missing nor hom-ref -- and write one single-sample VCF per
## sample.

workflow ExtractInformativeGenotypesPerSample {
  input {
    File vcf
    File vcf_idx
    File sample_list
    String docker = "staphb/bcftools:1.19"
  }

  Array[String] samples = read_lines(sample_list)

  scatter (sample in samples) {
    call ExtractSampleVariants {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        sample = sample,
        docker = docker,
    }
  }

  output {
    Array[File] sample_vcfs = ExtractSampleVariants.sample_vcf
    Array[File] sample_vcf_indices = ExtractSampleVariants.sample_vcf_idx
  }
}

task ExtractSampleVariants {
  input {
    File vcf
    File vcf_idx
    String sample
    String docker
  }

  Int disk_gb = ceil(size(vcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    ln -s ~{vcf} input.vcf.gz
    ln -s ~{vcf_idx} input.vcf.gz.tbi

    bcftools view -s ~{sample} -i 'GT="alt"' input.vcf.gz -Oz -o ~{sample}.informative.vcf.gz
    bcftools index -t ~{sample}.informative.vcf.gz
  >>>

  output {
    File sample_vcf = "~{sample}.informative.vcf.gz"
    File sample_vcf_idx = "~{sample}.informative.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}
