version 1.0

import "Structs.wdl"

workflow BAFTestChromosome {
  input {
    File vcf
    Array[String] samples
    String chrom
    File baf_metrics
    String batch
    String algorithm
    Int split_size
    Int? suffix_len
    File ref_dict

    String linux_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_baftest
    RuntimeAttr? runtime_attr_split_baf_vcf
    RuntimeAttr? runtime_attr_merge_baf
  }

  File baf_metrics_idx = baf_metrics + ".tbi"
  
  call SplitBafVcf {
    input:
      vcf = vcf,
      batch = batch,
      algorithm = algorithm,
      chrom = chrom,
      split_size = split_size,
      suffix_len = select_first([suffix_len, 4]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_split_baf_vcf
  }

  scatter (split in SplitBafVcf.split_beds) {
    call BAFTest {
      input:
        bed = split,
        baf_metrics = baf_metrics,
        baf_metrics_idx = baf_metrics_idx,
        ref_dict = ref_dict,
        samples = samples,
        prefix = basename(split),
        batch = batch,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_baftest
    }
  }

  call MergeBAFSplits {
    input:
      stats = BAFTest.stats,
      prefix = "${batch}.${algorithm}.${chrom}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_baf
  }

  output {
    File stats = MergeBAFSplits.merged_stats
  }
}

task BAFTest {
  input {
    File bed
    File baf_metrics
    File baf_metrics_idx
    File ref_dict
    Array[String] samples
    String prefix
    String batch
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    baf_metrics: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    File stats = "${prefix}.metrics"
  }
  command <<<

    set -euo pipefail

    echo -e "sample\tgroup\tbatch" > batch.key
    awk -v batch=~{batch} -v OFS="\t" '{print $1, $1, batch}' ~{write_lines(samples)} >> batch.key

    if [ -s ~{bed} ]; then
      set +o pipefail
      start=$(cut -f2 ~{bed} | sort -k1,1n | head -n1)
      end=$(cut -f3 ~{bed} | sort -k1,1n | tail -n1)
      chrom=$(cut -f1 ~{bed} | head -n1)
      set -o pipefail

      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --skip-header \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{baf_metrics} \
        -L "${chrom}:${start}-${end}" \
        -O local.BAF.txt.gz
    else
      touch local.BAF.txt
      bgzip local.BAF.txt
    fi

    tabix -s1 -b2 -e2 local.BAF.txt.gz
    svtk baf-test ~{bed} local.BAF.txt.gz --batch batch.key > ~{prefix}.metrics
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task MergeBAFSplits {
  input {
    Array[File] stats
    String prefix
    String linux_docker
    Int disk_gb_baseline = 10
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = disk_gb_baseline + ceil(2 * size(stats, "GiB"))
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File merged_stats = "${prefix}.stats"
  }
  command <<<

    set -eu
    echo -n "chrom start end name samples svtype delstat snp_ratio " > ~{prefix}.stats;
    echo -n "del_loglik dupstat KS_stat KS_pval total_case_snps " >> ~{prefix}.stats;
    echo -n "total_snps n_nonROH_cases n_samples mean_control_snps " >> ~{prefix}.stats;
    echo "n_nonROH_controls n_controls" >> ~{prefix}.stats;
    sed -i -e 's/ /\t/g' ~{prefix}.stats;
    while read split; do
      cat $split;
    done < ~{write_lines(stats)} >> ~{prefix}.stats
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: linux_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task SplitBafVcf {
  input {
    File vcf
    String batch
    String algorithm
    String chrom
    Int split_size
    Int suffix_len
    String sv_pipeline_docker
    Int disk_gb_baseline = 10
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = disk_gb_baseline + ceil(size(vcf, "GiB"))
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[File] split_beds = glob("${batch}.${algorithm}.split.*")
  }
  command <<<

    set -euo pipefail 
    tabix -p vcf ~{vcf}
    #TODO : split -a parameter should be scaled properly (using suffix_len does not always work)
    tabix -h ~{vcf} ~{chrom} \
      | svtk vcf2bed --no-header stdin stdout > all.bed
    if fgrep -q -e "DEL" -e "DUP" < all.bed; then
      fgrep -e "DEL" -e "DUP" < all.bed \
        | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
        | awk '($3-$2>=10000 && $3-$2<10000000)' \
        | split -a ~{suffix_len} -d -l 300 - ~{batch}.~{algorithm}.split.gt10kb.
      fgrep -e "DEL" -e "DUP" < all.bed \
        | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
        | awk '($3-$2<10000)' \
        | sort -k1,1V -k2,2n \
        | split -a ~{suffix_len} -d -l ~{split_size} - ~{batch}.~{algorithm}.split.
    fi
  
  >>>
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

