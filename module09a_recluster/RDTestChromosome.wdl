version 1.0

import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics

workflow RDTestChromosome {
  input {
    File vcf
    String algorithm
    File ped_file
    String chrom
    String batch
    File coveragefile
    Int split_size
    File medianfile
    String flags
    Int? suffix_len
    File male_samples
    File female_samples
    File male_only_variant_ids
    File samples
    Boolean allosome
    File ref_dict

    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String linux_docker
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_split_rd_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
  }
  
  File coveragefile_idx = coveragefile + ".tbi"

  call SplitRDVcf {
    input:
      vcf = vcf,
      batch = batch,
      algorithm = algorithm,
      chrom = chrom,
      split_size = split_size,
      suffix_len = select_first([suffix_len, 4]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_split_rd_vcf
  }

  scatter (split in SplitRDVcf.split_beds) {
    if (allosome) {
      call RDTest as RDTestFemale {
        input:
          bed = split,
          coveragefile = coveragefile,
          coveragefile_idx = coveragefile_idx,
          medianfile = medianfile,
          ped_file = ped_file,
          include_list = female_samples,
          prefix = basename(split),
          flags = flags,
          ref_dict = ref_dict,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest
      }

      call RDTest as RDTestMale {
        input:
          bed = split,
          coveragefile = coveragefile,
          coveragefile_idx = coveragefile_idx,
          medianfile = medianfile,
          ped_file = ped_file,
          include_list = male_samples,
          prefix = basename(split),
          flags = flags,
          ref_dict = ref_dict,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest
      }

      call tasksbatchmetrics.MergeAllosomes as MergeAllosomes {
        input:
          male_test = RDTestMale.stats,
          female_test = RDTestFemale.stats,
          male_only_ids_list = male_only_variant_ids,
          chrom = chrom,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_merge_allo,
          male_only_expr = "females.CNVID.isin(male_only_ids)"
      }
    }

    if (!allosome) {
      call RDTest as RDTestAutosome {
        input:
          bed = split,
          coveragefile = coveragefile,
          coveragefile_idx = coveragefile_idx,
          medianfile = medianfile,
          include_list = samples,
          ped_file = ped_file,
          prefix = basename(split),
          flags = flags,
          ref_dict = ref_dict,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest
      }
    }
  }

  Array[File?] stats = if allosome then MergeAllosomes.merged_test else RDTestAutosome.stats

  call MergeRDSplits {
    input:
      stats = select_all(stats),
      prefix = "${batch}.${algorithm}.${chrom}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  output {
    File out_stats = MergeRDSplits.merged_stats
  }
}

task RDTest {
  input {
    File bed
    File coveragefile
    File coveragefile_idx
    File medianfile
    File ped_file
    File include_list
    File ref_dict
    String prefix
    String flags
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    coveragefile: {
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
    File stats = "~{prefix}.metrics"
  }
  command <<<

    set -euo pipefail

    if [ -s ~{bed} ]; then
      set +o pipefail
      start=$(( $(cut -f2 ~{bed} | sort -k1,1n | head -n1) ))
      end=$(( $(cut -f3 ~{bed} | sort -k1,1n | tail -n1) ))
      chrom=$(cut -f1 ~{bed} | head -n1)
      set -o pipefail

      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{coveragefile} \
        -L "${chrom}:${start}-${end}" \
        -O local.RD.txt.gz
    else
      touch local.RD.txt
      bgzip local.RD.txt
    fi

    tabix -p bed local.RD.txt.gz

    Rscript /opt/RdTest/RdTest.R \
      -b ~{bed} \
      -n ~{prefix} \
      -c local.RD.txt.gz \
      -m ~{medianfile} \
      -f ~{ped_file} \
      -w ~{include_list} \
      ~{flags}
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_rdtest_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task SplitRDVcf {
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
    tabix -h ~{vcf} ~{chrom} \
      | svtk vcf2bed --no-header stdin stdout > all.bed
    if fgrep -q -e "DEL" -e "DUP" < all.bed; then
      fgrep -e "DEL" -e "DUP" < all.bed \
        | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
        | awk '($3-$2>=10000)' \
        > ~{batch}.~{algorithm}.split.gt10kb
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

task MergeRDSplits {
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
    echo 'chr	Start	End	CNVID	SampleIDs	Type	Median_Power	P	2ndMaxP	Model	Median_Rank	Median_Separation' > ~{prefix}.stats
    while read split; do
      sed 1d $split;
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
