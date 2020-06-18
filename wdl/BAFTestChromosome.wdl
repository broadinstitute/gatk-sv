##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_baftest_autosome/15/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

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
    Int tabix_retries

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
        samples = samples,
        prefix = basename(split),
        batch = batch,
        tabix_retries = tabix_retries,
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
    Array[String] samples
    String prefix
    String batch
    Int tabix_retries
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

  output {
    File stats = "${prefix}.metrics"
  }
  command <<<

    set -euxo pipefail
    echo -e "sample\tgroup\tbatch" > batch.key
    awk -v batch=~{batch} -v OFS="\t" '{print $1, $1, batch}' ~{write_lines(samples)} >> batch.key
    set +o pipefail
    start=$(cut -f2 ~{bed} | sort -k1,1n | head -n1)
    end=$(cut -f3 ~{bed} | sort -k1,1n | tail -n1)
    chrom=$(cut -f1 ~{bed} | head -n1)
    set -o pipefail

    # Temporary workaround for corrupted tabix downloads
    x=0
    while [ $x -lt ~{tabix_retries} ]
    do
      # Download twice
      GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ~{baf_metrics} "$chrom":"$start"-"$end" | bgzip -c > local_baf_1.bed.gz
      GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix ~{baf_metrics} "$chrom":"$start"-"$end" | bgzip -c > local_baf_2.bed.gz
      # Done if the downloads were identical, otherwise retry
      cmp --silent local_baf_1.bed.gz local_baf_2.bed.gz && break       
      x=$(( $x + 1))
    done
    echo "BAF-tabix retry:" $x
    cmp --silent local_baf_1.bed.gz local_baf_2.bed.gz || exit 1

    mv local_baf_1.bed.gz local_baf.bed.gz
    tabix -b2 local_baf.bed.gz
    svtk baf-test ~{bed} local_baf.bed.gz --batch batch.key > ~{prefix}.metrics
  
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

task MergeBAFSplits {
  input {
    Array[File] stats
    String prefix
    String linux_docker
    RuntimeAttr? runtime_attr_override
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
    RuntimeAttr? runtime_attr_override
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

  output {
    Array[File] split_beds = glob("${batch}.${algorithm}.split.*")
  }
  command <<<

    set -euo pipefail 
    tabix -p vcf ~{vcf};
    #TODO : split -a parameter should be scaled properly (using suffix_len does not always work)
    tabix -h ~{vcf} ~{chrom} \
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2>=10000 && $3-$2<10000000)' \
      | split -a ~{suffix_len} -d -l 300 - ~{batch}.~{algorithm}.split.gt10kb.
    tabix -h ~{vcf} ~{chrom} \
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2<10000)' \
      | sort -k1,1V -k2,2n \
      | split -a ~{suffix_len} -d -l ~{split_size} - ~{batch}.~{algorithm}.split.
  
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

