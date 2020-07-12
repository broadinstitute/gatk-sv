##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/02_aggregate_MMDLW/3/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "PETest.wdl" as pet
import "RDTest.wdl" as rdt
import "SRTest.wdl" as srt
import "BAFTest.wdl" as baft
import "Tasks02.wdl" as tasks02
import "Utils.wdl" as util

workflow Module02 {
  input {
    String batch

    File? depth_vcf
    File? melt_vcf
    File? delly_vcf
    File? wham_vcf
    File? manta_vcf

    File baf_metrics
    File discfile
    File coveragefile
    File splitfile
    File medianfile

    Int BAF_split_size
    Int RD_split_size
    Int PE_split_size
    Int SR_split_size
    Int common_cnv_size_cutoff

    File rmsk
    File segdups
    File ped_file
    File autosome_contigs
    File allosome_contigs
    File ref_dict

    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_base_docker
    String linux_docker

    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_sample_list
    RuntimeAttr? runtime_attr_baf_samples
    RuntimeAttr? runtime_attr_aggregate_tests
    RuntimeAttr? runtime_attr_aggregate_callers
    RuntimeAttr? runtime_attr_petest
    RuntimeAttr? runtime_attr_srtest
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_baftest
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_split_rd_vcf
    RuntimeAttr? runtime_attr_split_baf_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_baf
    RuntimeAttr? runtime_attr_merge_stats
  }

  Array[String] algorithms = ["depth", "melt", "delly", "wham", "manta"]
  Array[File?] vcfs = [depth_vcf, melt_vcf, delly_vcf, wham_vcf, manta_vcf]

  call util.GetSampleIdsFromVcf {
    input:
      vcf = select_first(vcfs),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  call GetSampleLists {
    input:
      ped_file = ped_file,
      samples = GetSampleIdsFromVcf.out_array,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_sample_list
  }

  scatter (i in range(length(algorithms))) {
    
    if (defined(vcfs[i])) {

      String algorithm = algorithms[i]
      File vcf = select_first([vcfs[i]])

      if (algorithm != "melt") {
        call rdt.RDTest as RDTest {
          input:
            coveragefile = coveragefile,
            medianfile = medianfile,
            ped_file = ped_file,
            vcf = vcf,
            autosome_contigs = autosome_contigs,
            split_size = RD_split_size,
            flags = "",
            algorithm = algorithm,
            allosome_contigs = allosome_contigs,
            ref_dict = ref_dict,
            batch = batch,
            samples = GetSampleLists.samples_file,
            male_samples = GetSampleLists.male_samples,
            female_samples = GetSampleLists.female_samples,
            sv_pipeline_docker = sv_pipeline_docker,
            sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
            linux_docker = linux_docker,
            runtime_attr_rdtest = runtime_attr_rdtest,
            runtime_attr_split_rd_vcf = runtime_attr_split_rd_vcf,
            runtime_attr_merge_allo = runtime_attr_merge_allo,
            runtime_attr_merge_stats = runtime_attr_merge_stats
        }

        call baft.BAFTest as BAFTest {
          input:
            baf_metrics = baf_metrics,
            vcf = vcf,
            autosome_contigs = autosome_contigs,
            ref_dict = ref_dict,
            split_size = BAF_split_size,
            algorithm = algorithm,
            batch = batch,
            samples = GetSampleIdsFromVcf.out_array,
            linux_docker = linux_docker,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_baftest = runtime_attr_baftest,
            runtime_attr_split_baf_vcf = runtime_attr_split_baf_vcf,
            runtime_attr_merge_baf = runtime_attr_merge_baf,
            runtime_attr_merge_stats = runtime_attr_merge_stats
        }
      }

      if (algorithm != "depth") {
        call srt.SRTest as SRTest {
          input:
            splitfile = splitfile,
            medianfile = medianfile,
            ped_file = ped_file,
            vcf = vcf,
            autosome_contigs = autosome_contigs,
            ref_dict = ref_dict,
            split_size = SR_split_size,
            algorithm = algorithm,
            allosome_contigs = allosome_contigs,
            batch = batch,
            samples = GetSampleLists.samples_file,
            male_samples = GetSampleLists.male_samples,
            female_samples = GetSampleLists.female_samples,
            run_common = true,
            common_cnv_size_cutoff = common_cnv_size_cutoff,
            sv_base_mini_docker = sv_base_mini_docker,
            linux_docker = linux_docker,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_srtest = runtime_attr_srtest,
            runtime_attr_split_vcf = runtime_attr_split_vcf,
            runtime_attr_merge_allo = runtime_attr_merge_allo,
            runtime_attr_merge_stats = runtime_attr_merge_stats
        }
      }

      if (algorithm != "depth" && algorithm != "melt") {
        call pet.PETest as PETest {
          input:
            discfile = discfile,
            medianfile = medianfile,
            ped_file = ped_file,
            vcf = vcf,
            autosome_contigs = autosome_contigs,
            ref_dict = ref_dict,
            split_size = PE_split_size,
            algorithm = algorithm,
            allosome_contigs = allosome_contigs,
            batch = batch,
            samples = GetSampleLists.samples_file,
            male_samples = GetSampleLists.male_samples,
            female_samples = GetSampleLists.female_samples,
            common_cnv_size_cutoff = common_cnv_size_cutoff,
            sv_base_mini_docker = sv_base_mini_docker,
            linux_docker = linux_docker,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_petest = runtime_attr_petest,
            runtime_attr_split_vcf = runtime_attr_split_vcf,
            runtime_attr_merge_allo = runtime_attr_merge_allo,
            runtime_attr_merge_stats = runtime_attr_merge_stats
        }
      }

      call AggregateTests {
        input:
          vcf = vcf,
          petest = PETest.petest,
          srtest = SRTest.srtest,
          rdtest = RDTest.rdtest,
          baftest = BAFTest.baftest,
          segdups = segdups,
          rmsk = rmsk,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_aggregate_tests
      }

      call tasks02.GetCommonVCF {
        input:
          vcf = vcf,
          cnv_size_cutoff = common_cnv_size_cutoff,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_split_vcf
      }

      call AggregateTests as AggregateTestsCommon {
        input:
          vcf = GetCommonVCF.common_vcf,
          petest = PETest.petest_common,
          srtest = SRTest.srtest_common,
          segdups = segdups,
          rmsk = rmsk,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_aggregate_tests
      }
    }
  }

  call AggregateCallers {
    input:
      batch = batch,
      input_metrics = select_all(AggregateTests.metrics),
      common = false,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_aggregate_callers
  }

  call AggregateCallers as AggregateCallersCommon {
    input:
      batch = batch,
      input_metrics = select_all(AggregateTestsCommon.metrics),
      common = true,
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_aggregate_callers
  }

  output {
    File metrics = AggregateCallers.metrics
    File metrics_common = AggregateCallersCommon.metrics
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

task AggregateTests {
  input {
    File vcf
    File? rdtest
    File? baftest
    File? petest
    File? srtest
    File segdups
    File rmsk
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 7.5,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File metrics = "aggregated.metrics"
  }
  command <<<

    /opt/sv-pipeline/02_evidence_assessment/02e_metric_aggregation/scripts/aggregate.py \
      -v ~{vcf} \
      ~{if defined(rdtest) then "-r ~{rdtest}" else "" } \
      ~{if defined(baftest) then "-b ~{baftest}" else "" } \
      ~{if defined(petest) then "-p ~{petest}" else "" } \
      ~{if defined(srtest) then "-s ~{srtest}" else "" } \
      --segdups ~{segdups} \
      --rmsk ~{rmsk} \
      aggregated.metrics
  
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

task AggregateCallers {
  input {
    String batch
    Array[File] input_metrics
    Boolean common
    String sv_pipeline_base_docker
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

  String output_file = if common then "${batch}.common.metrics" else "${batch}.metrics"

  output {
    File metrics = "~{output_file}"
  }
  command <<<

    set -eu
    python3 <<CODE
    import pandas as pd
    metrics = ["~{sep='", "' input_metrics}"]
    dfs=[]
    for df in metrics:
      dfs.append(pd.read_table(df))
    df = pd.concat(dfs)
    df.to_csv("~{output_file}", index=False, sep='\t')
    CODE
        
  >>>
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
