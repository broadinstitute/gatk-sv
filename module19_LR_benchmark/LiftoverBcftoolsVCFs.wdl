version 1.0

import "Structs.wdl"

workflow BcftoolsLiftoverVCFs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File chain_file

    File fasta_t2t
    File fasta_hg38
    String liftover_ref_version

    String bcftools_liftover_docker
    RuntimeAttr? runtime_attr_liftover
  }

  scatter(i in range(length(vcfs))) {
    call LiftoverBcftools {
        input:
            vcf = vcfs[i],
            vcf_idx = vcf_idxs[i],

            fasta_t2t = fasta_t2t,
            fasta_hg38 = fasta_hg38,
            chain_file = chain_file,
            liftover_ref_version = liftover_ref_version,

            docker_image = bcftools_liftover_docker,
            runtime_attr_override = runtime_attr_liftover
    }
  }

    output {
        Array[File] lifted_sorted_vcfs = LiftoverBcftools.lifted_vcf
        Array[File] lifted_sorted_tbis = LiftoverBcftools.lifted_vcf_tbi
    }
}



task LiftoverBcftools {
    input {
        File vcf
        File vcf_idx
        File fasta_t2t
        File fasta_hg38
        File chain_file
        String liftover_ref_version
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(bed, ".bed")

    command <<<
        set -euo pipefail

        bcftools +liftover -Ou ~{vcf} -- -s ~{fasta_t2t} -f ~{fasta_hg38} \
        -c ~{chain_file} | bcftools sort -Oz -o ~{prefix}.~{liftover_ref_version}.liftover.vcf.gz -W=tbi

    >>>

    output {
        File lifted_vcf = "~{prefix}.~{liftover_ref_version}.liftover.vcf.gz"
        File lifted_vcf_tbi = "~{prefix}.~{liftover_ref_version}.liftover.vcf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1,
        mem_gb: 10,
        disk_gb: 15 + ceil(size(bed, "GiB") *5),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: docker_image
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}


