##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_srtest_allosome/12/wdl
## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_srtest_autosome/13/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Tasks02.wdl" as tasks02

workflow RDTestChromosome {
  input {
    File bed
    File ped_file
    String chrom
    String batch
    File coveragefile
    Int split_size
    File medianfile
    String flags
    Int? suffix_len
    Boolean allosome
    Array[String] samples
    String output_gs_folder

    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String linux_docker
    String sv_base_docker
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_split_rd_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
  }
  
  File coveragefile_idx = coveragefile + ".tbi"

  call GetSampleLists {
    input:
      ped_file = ped_file,
      samples = samples,
      sv_base_docker = sv_base_docker
  }

  call SplitRDbed {
    input:
      bed = bed,
      batch = batch,
      chrom = chrom,
      split_size = split_size,
      suffix_len = select_first([suffix_len, 4]),
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_split_rd_vcf
  }

  scatter (split in SplitRDbed.split_beds) {
    if (allosome) {
      call RDTest as RDTestFemale {
        input:
          bed = split,
          coveragefile = coveragefile,
          coveragefile_idx = coveragefile_idx,
          medianfile = medianfile,
          ped_file = ped_file,
          whitelist = GetSampleLists.female_samples,
          prefix = basename(split),
          flags = flags,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest,
          output_gs_folder = output_gs_folder
      }

      call RDTest as RDTestMale {
        input:
          bed = split,
          coveragefile = coveragefile,
          coveragefile_idx = coveragefile_idx,
          medianfile = medianfile,
          ped_file = ped_file,
          whitelist = GetSampleLists.male_samples,
          prefix = basename(split),
          flags = flags,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest,
          output_gs_folder = output_gs_folder
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
          whitelist = GetSampleLists.samples_file,
          ped_file = ped_file,
          prefix = basename(split),
          flags = flags,
          sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
          runtime_attr_override = runtime_attr_rdtest,
          output_gs_folder = output_gs_folder
      }
    }
  }

  Array[File?] stats = if allosome then MergeAllosomes.merged_test else RDTestAutosome.stats

  call tasks02.MergeStats as MergeStats {
    input:
      stats = select_all(stats),
      prefix = "${batch}.${chrom}",
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_merge_stats
  }

  output {
    File stats = MergeStats.merged_stats
  }
}


task GetSampleLists {
  input {
    File ped_file
    Array[String] samples
    String sv_base_docker
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

  File samples_list = write_lines(samples)

  output {
    File male_samples = "male.list"
    File female_samples = "female.list"
    File samples_file = "samples.list"
  }
  command <<<

    set -eu
    awk -v sex=1 '($5==sex) {print $2}' ~{ped_file} > ped_males.list
    awk -v sex=2 '($5==sex) {print $2}' ~{ped_file} > ped_females.list
    cat ~{samples_list} > samples.list

    python3 <<CODE
    with open("ped_males.list",'r') as ped_m, open("ped_females.list",'r') as ped_f:
      male_samples = set([x.strip() for x in ped_m.readlines() if x.strip()])
      female_samples = set([x.strip() for x in ped_f.readlines() if x.strip()])
      with open("male.list", 'w') as samples_m, open("female.list",'w') as samples_f, open("samples.list",'r') as samples:
        for line in samples:
          if line.strip():
            if (line.strip() in male_samples):
              samples_m.write(line)
            if (line.strip() in female_samples):
              samples_f.write(line)
    CODE
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}



task RDTest {
  input {
    File bed
    File coveragefile
    File coveragefile_idx
    File medianfile
    File ped_file
    File whitelist
    String prefix
    String flags
    String sv_pipeline_rdtest_docker
    String output_gs_folder
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    coveragefile: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 10, 
    disk_gb: 20,
    boot_disk_gb: 20,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File stats = "./tmp/~{prefix}.metrics"
  }
  command <<<

    set -eu
    start=$(cut -f2 ~{bed} | sort -k1,1n | head -n1);
    end=$(cut -f3 ~{bed} | sort -k1,1n | tail -n1);
    chrom=$(cut -f1 ~{bed} | head -n1);
    set -o pipefail
    GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token` \
      tabix -h ~{coveragefile} "$chrom":"$start"-"$end" \
      | sed 's/Chr/chr/g' \
      | sed 's/Start/start/g' \
      | sed 's/End/end/' \
      | bgzip -c > local_coverage.bed.gz
    tabix -p bed local_coverage.bed.gz;
    mkdir ./tmp;
    Rscript /opt/RdTest/RdTest.R \
      -b ~{bed} \
      -n ~{prefix} \
      -c local_coverage.bed.gz \
      -m ~{medianfile} \
      -f ~{ped_file} \
      -w ~{whitelist} \
      ~{flags} \
      --plot TRUE \
      -o  ./tmp
      gsutil cp ./tmp/* ~{output_gs_folder}

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
    String chrom
    Int split_size
    Int suffix_len
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 10, 
    disk_gb: 20,
    boot_disk_gb: 20,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[File] split_beds = glob("${batch}.split.*")
  }
  command <<<

    set -euo pipefail
    tabix -p vcf ~{vcf};
    tabix -h ~{vcf} ~{chrom} \
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2>=10000)' \
      > ~{batch}.split.gt10kb;
    tabix -h ~{vcf} ~{chrom} \
      | svtk vcf2bed --no-header stdin stdout \
      | fgrep -e "DEL" -e "DUP" \
      | awk -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' \
      | awk '($3-$2<10000)' \
      | sort -k1,1V -k2,2n \
      | split -a ~{suffix_len} -d -l ~{split_size} - ~{batch}.split.
  
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

task SplitRDbed {
  input {
    File bed
    String batch
    String chrom
    Int split_size
    Int suffix_len
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 10, 
    disk_gb: 20,
    boot_disk_gb: 20,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[File] split_beds = glob("${batch}.split.*")
  }
  command <<<

    set -euo pipefail
    split -a ~{suffix_len} -d -l ~{split_size} ~{bed} ~{batch}.split.
  
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

