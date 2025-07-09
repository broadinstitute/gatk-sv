version 1.0

import "Structs.wdl"

workflow MergeVcfsByChromosome {
  input {
    String chrom
    Array[File] input_vcfs
    Array[File] input_vcfs_idx
    Array[String] sample_list
    Boolean convert_to_biallelic = false
    String sv_base_mini_docker
    String sv_pipeline_base_docker
  }


    scatter (idx in range(length(input_vcfs))) {

      call ExtractChromosomeVcf {
        input:
          input_vcf = input_vcfs[idx],
          input_vcf_idx = input_vcfs_idx[idx],
          chromosome = chrom,
          sv_base_mini_docker = sv_base_mini_docker
      }
    }

    call MergeVcfs {
      input:
        input_vcfs = ExtractChromosomeVcf.output_vcf,
        input_vcfs_idx = ExtractChromosomeVcf.output_vcf_idx,
        output_name = "${chrom}.vcf.gz",
        sv_base_mini_docker = sv_base_mini_docker

    }

    if(convert_to_biallelic){
      call ConvertMultiToBiAllelic {
        input:
          vcf = MergeVcfs.output_merged_vcf,
          vcf_idx = MergeVcfs.output_merged_vcf_idx,
          sv_pipeline_base_docker = sv_pipeline_base_docker
        }
      }




  output {
    File merged_vcf = select_first([ConvertMultiToBiAllelic.bi_allelic_vcf , MergeVcfs.output_merged_vcf])
    File merged_vcf_idx = select_first([ConvertMultiToBiAllelic.bi_allelic_vcf_idx, MergeVcfs.output_merged_vcf_idx])
  }
}

# Task 1: Extract a chromosome from a VCF
task ExtractChromosomeVcf {
  input {
    File input_vcf
    File input_vcf_idx
    String chromosome
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(input_vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    set -e
    bcftools view -r ~{chromosome} ~{input_vcf} -Oz -o ~{chromosome}.vcf.gz
    tabix -p vcf ~{chromosome}.vcf.gz
  >>>

  output {
    File output_vcf = "~{chromosome}.vcf.gz"
    File output_vcf_idx = "~{chromosome}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Task 2: Merge multiple VCFs
task MergeVcfs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_idx
    String output_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10,
    disk_gb: ceil(10 + size(input_vcfs, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -e
    bcftools merge --merge none ~{sep=' ' input_vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}
  >>>

  output {
    File output_merged_vcf = output_name
    File output_merged_vcf_idx = "${output_name}.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Task 3: Concatenate per-chromosome VCFs
task ConcatVcfs {
  input {
    Array[File] input_vcfs
    Array[File] input_vcfs_idx
    String output_name
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10,
    disk_gb: ceil(10 + size(input_vcfs, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    set -e
    bcftools concat ~{sep=' ' input_vcfs} -Oz -o ~{output_name}
    tabix -p vcf ~{output_name}
  >>>

  output {
    File output_vcf = output_name
    File output_vcf_idx = "${output_name}.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Task 4: index VCFs
task IndexVcf {
  input {
    File vcf                # input VCF (.vcf.gz)
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: ceil(size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    set -e
      tabix -p vcf ~{vcf}
  >>>

  output {
    File indexed_vcf_idx = "~{vcf}.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Task 5: convert multi-allelic vcf to bi-allelic vcf
task ConvertMultiToBiAllelic {
  input {
    File vcf                # input VCF (.vcf.gz)
    File vcf_idx
    String sv_pipeline_base_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf, "GB") * 3),
    disk_gb: 20 + ceil(size(vcf, "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String prefix = basename(vcf, ".vcf.gz")
  command <<<
    set -e

    python3 <<CODE

    import pysam
    import sys

    def split_multiallelic_to_biallelic(input_vcf, output_vcf):
      vcf_in = pysam.VariantFile(input_vcf)
      vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)
        for rec in vcf_in:
            if len(rec.alts) <= 1:
                vcf_out.write(rec)
                continue
            for alt_index, alt_allele in enumerate(rec.alts):
                new_rec = rec
                rec.alts = alts=(alt_allele,)
                for sample in rec.samples:
                    orig_gt = rec.samples[sample]["GT"]
                    if orig_gt is None:
                        new_rec.samples[sample]["GT"] = (None, None)
                        continue
                    # remap genotype to bi-allelic
                    new_gt = []
                    for allele in orig_gt:
                        if allele == 0:
                            new_gt.append(0)
                        elif allele == alt_index + 1:
                            new_gt.append(1)
                        elif allele is None:
                            new_gt.append(None)
                        else:
                            new_gt.append(0)  # non-matching ALT treated as REF
                    new_rec.samples[sample]["GT"] = tuple(new_gt)
                vcf_out.write(new_rec)
        vcf_in.close()
        vcf_out.close()
    split_multiallelic_to_biallelic("~{vcf}", "~{prefix}.biallelic.vcf.gz")

    CODE

    bcftools sort "~{prefix}.biallelic.vcf.gz" -Oz -o "~{prefix}.biallelic.sorted.vcf.gz" 
    tabix -p vcf "~{prefix}.biallelic.sorted.vcf.gz" 

  >>>

  output {
    File bi_allelic_vcf = "~{prefix}.biallelic.sorted.vcf.gz" 
    File bi_allelic_vcf_idx = "~{prefix}.biallelic.sorted.vcf.gz.tbi" 

  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}





