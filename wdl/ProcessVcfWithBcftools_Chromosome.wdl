
workflow BcftoolsProcessing {
    input {
        Array[String] combined_vcfs
        Int cpus = 3
    }

    scatter (vcf_path in combined_vcfs) {
        call ProcessVCF {
            input:
                vcf_gcs_path = vcf_path,
                cpus = cpus
        }
    }

    output {
        Array[File] processed_vcfs = ProcessVCF.output_vcf
    }
}

task ProcessVCF {
    input {
        String vcf_gcs_path
        Int cpus
    }

    # Extract chromosome name from the GCS path
    String chrom = sub(vcf_gcs_path, ".*chr([0-9XYM]+).*", "$1")

    command {
        gsutil cat ~{vcf_gcs_path} | \
        bcftools view -i '(SVTYPE=="DEL" || SVTYPE=="DUP") && SVLEN>=1000000' -Ou - | \
        bcftools +fill-tags - -O z -- -t AF | \
        bcftools view -i 'AF>0.005' -O z -G --threads ~{cpus} - > processed_~{chrom}.vcf.gz
    }

    output {
        File output_vcf = "processed_~{chrom}.vcf.gz"
    }

    runtime {
        cpu: ~{cpus}
        memory: "3GB"
        disks: "20GB"
        docker: "broadinstitute/genomes-in-the-cloud:2.4.3-1552931386"
    }
}
