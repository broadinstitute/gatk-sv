version 1.0

import "Structs.wdl"

workflow UpdateVcfGenotypes {
  input {
    File vcf_A
    File vcf_A_index
    File vcf_B
    File vcf_B_index
    File patch_script          # patch_vcf.py
    String output_prefix
    Array[String] contigs
    String docker_image
  }

  scatter (ctg in contigs) {
    call ExtractContigVCF {
      input:
        vcf = vcf_A,
        vcf_index = vcf_A_index,
        contig = ctg,
        docker_image = docker_image
    }

    call PatchGenotypes {
      input:
        vcf_A_contig = ExtractContigVCF.out_vcf,
        vcf_B = vcf_B,
        vcf_B_index = vcf_B_index,
        patch_script = patch_script,
        contig = ctg,
        docker_image = docker_image
    }
  }

  call ConcatPatchedVCFs{
    input:
      vcfs = PatchGenotypes.out_vcf,
      output_prefix = output_prefix, 
      docker_image = docker_image
  }

  output {
    File updated_vcf = ConcatPatchedVCFs.out_vcf
    File updated_vcf_idx = ConcatPatchedVCFs.out_vcf_index
  }
}

task ExtractContigVCF {
  input {
    File vcf
    File vcf_index
    String contig
    String docker_image
  }

  command <<<
    set -euo pipefail

    bcftools view -r ~{contig} -Oz -o ~{contig}.A.vcf.gz ~{vcf}
    tabix -p vcf ~{contig}.A.vcf.gz
  >>>

  output {
    File out_vcf = "~{contig}.A.vcf.gz"
    File out_vcf_index = "~{contig}.A.vcf.gz.tbi"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 20 HDD"
  }
}

task PatchGenotypes {
  input {
    File vcf_A_contig
    File vcf_B
    File vcf_B_index
    File patch_script
    String contig
    String docker_image
  }

  command <<<
    set -euo pipefail

    python3 ~{patch_script} \
      ~{vcf_A_contig} \
      ~{vcf_B} \
      ~{contig}.patched.vcf.gz

    tabix -p vcf ~{contig}.patched.vcf.gz
  >>>

  output {
    File out_vcf = "~{contig}.patched.vcf.gz"
    File out_vcf_index = "~{contig}.patched.vcf.gz.tbi"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "6 GiB"
    disks: "local-disk 30 HDD"
  }
}

task ConcatPatchedVCFs {
  input {
    Array[File] vcfs
    String output_prefix
    String docker_image
  }

  command <<<
    set -euo pipefail

    # Ensure deterministic order
    printf "%s\n" ~{sep='\n' vcfs} > vcf.list

    bcftools concat \
      -f vcf.list \
      -Oz \
      -o ~{output_prefix}.vcf.gz

    tabix -p vcf ~{output_prefix}.vcf.gz
  >>>

  output {
    File out_vcf = "~{output_prefix}.vcf.gz"
    File out_vcf_index = "~{output_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "8 GiB"
    disks: "local-disk 50 HDD"
  }
}

