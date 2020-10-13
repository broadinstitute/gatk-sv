##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_srtest_allosome/12/wdl
## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_srtest_autosome/13/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Tasks02.wdl" as tasks02

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
    File samples
    Boolean allosome
    Int tabix_retries

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
          tabix_retries = tabix_retries,
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
          tabix_retries = tabix_retries,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest
      }

      call tasks02.MergeAllosomes as MergeAllosomes {
        input:
          male_test = RDTestMale.stats,
          female_test = RDTestFemale.stats,
          chrom = chrom,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_merge_allo,
          male_only_expr = "females.P.astype(str) == 'No_samples_for_analysis'"
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
          tabix_retries = tabix_retries,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest
      }
    }
  }

  Array[File?] stats = if allosome then MergeAllosomes.merged_test else RDTestAutosome.stats

  call tasks02.MergeStats as MergeStats {
    input:
      stats = select_all(stats),
      prefix = "${batch}.${algorithm}.${chrom}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  output {
    File stats = MergeStats.merged_stats
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
    String prefix
    String flags
    Int tabix_retries
    String sv_pipeline_rdtest_docker
    Int disk_gb_baseline = 10
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    coveragefile: {
      localization_optional: true
    }
  }

  Int disk_gb = disk_gb_baseline + ceil(
                                    size([bed, medianfile, ped_file, include_list], "GiB")
                                    + 2 * size([coveragefile, coveragefile_idx], "GiB")
                                 )
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
    File stats = "~{prefix}.metrics"
  }
  command <<<

    set -eu
    start=$(cut -f2 ~{bed} | sort -k1,1n | head -n1);
    end=$(cut -f3 ~{bed} | sort -k1,1n | tail -n1);
    chrom=$(cut -f1 ~{bed} | head -n1);
    set -o pipefail

    # Temporary workaround for corrupted tabix downloads
    x=0
    while [ $x -lt ~{tabix_retries} ]
    do
      # Download twice
      GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -h ~{coveragefile} "$chrom":"$start"-"$end" | sed 's/Chr/chr/g' | sed 's/Start/start/g' | sed 's/End/end/' | bgzip -c > local_coverage_1.bed.gz
      GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` tabix -h ~{coveragefile} "$chrom":"$start"-"$end" | sed 's/Chr/chr/g' | sed 's/Start/start/g' | sed 's/End/end/' | bgzip -c > local_coverage_2.bed.gz
      # Done if the downloads were identical, otherwise retry
      cmp --silent local_coverage_1.bed.gz local_coverage_2.bed.gz && break       
      x=$(( $x + 1))
    done
    echo "RD-tabix retry:" $x
    cmp --silent local_coverage_1.bed.gz local_coverage_2.bed.gz || exit 1
    
    mv local_coverage_1.bed.gz local_coverage.bed.gz
    tabix -p bed local_coverage.bed.gz;

    Rscript /opt/RdTest/RdTest.R \
      -b ~{bed} \
      -n ~{prefix} \
      -c local_coverage.bed.gz \
      -m ~{medianfile} \
      -f ~{ped_file} \
      -w ~{include_list} \
      ~{flags}
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
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
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2>=10000)' \
      > ~{batch}.~{algorithm}.split.gt10kb;
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

