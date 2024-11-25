version 1.0

import "Structs.wdl"
import "Utils.wdl" as utils

workflow PreprocessPESR {
  input {
    Array[String] samples         # Sample ID
    Array[File]? dragen_vcfs       # Dragen VCF
    Array[File]? manta_vcfs        # Manta VCF
    Array[File]? melt_vcfs         # Melt VCF
    Array[File]? scramble_vcfs     # Scramble VCF
    Array[File]? wham_vcfs         # Wham VCF
    File contigs          # .fai file of included contigs
    Int min_svsize        # Minimum SV length to include
    String batch
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr
  }

  Array[String] algorithms = ["dragen", "manta", "melt", "scramble", "wham"]
  Array[Array[File]?] vcfs = [dragen_vcfs, manta_vcfs, melt_vcfs, scramble_vcfs, wham_vcfs]

  scatter (i in range(length(algorithms))) {
    if (defined(vcfs[i]) && (length(select_first([vcfs[i]])) > 0)) {
      String algorithm = algorithms[i]
      call StandardizeVCFs {
        input:
          raw_vcfs=select_first([vcfs[i]]),
          caller = algorithm,
          samples=samples,
          contigs=contigs,
          min_svsize=min_svsize,
          prefix="~{batch}.~{algorithm}_std",
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override=runtime_attr
      }
    }
  }
  output {
    File? std_dragen_vcf_tar = StandardizeVCFs.out[0]
    File? std_manta_vcf_tar = StandardizeVCFs.out[1]
    File? std_melt_vcf_tar = StandardizeVCFs.out[2]
    File? std_scramble_vcf_tar = StandardizeVCFs.out[3]
    File? std_wham_vcf_tar = StandardizeVCFs.out[4]
  }
}

task StandardizeVCFs {
  input {
    Array[File] raw_vcfs
    Array[String] samples
    String caller
    File contigs
    Int min_svsize
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int num_samples = length(samples)

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 3
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    vcfs=(~{sep=" " raw_vcfs})
    sample_ids=(~{sep=" " samples})
    mkdir out
    for (( i=0; i<~{num_samples}; i++ ));
    do
      vcf=${vcfs[$i]}
      sample_id=${sample_ids[$i]}
      svtk standardize --sample-names ${sample_id} --prefix ~{caller}_${sample_id} --contigs ~{contigs} --min-size ~{min_svsize} $vcf tmp.vcf ~{caller}
      sample_no=`printf %03d $i`
      bcftools sort tmp.vcf -Oz -o out/std_${sample_no}.~{caller}.${sample_id}.vcf.gz
    done
    tar czf ~{prefix}.tar.gz -C out/ .
  >>>

  output {
    File out = "~{prefix}.tar.gz"
  }
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
