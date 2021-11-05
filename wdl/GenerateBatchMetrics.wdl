version 1.0

import "PETest.wdl" as pet
import "SRTest.wdl" as srt
import "TasksGenerateBatchMetrics.wdl" as tasksbatchmetrics
import "Utils.wdl" as util
import "GenerateBatchMetricsMetrics.wdl" as metrics

workflow GenerateBatchMetrics {
  input {
    String batch

    File original_metrics

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
    String? chr_x

    # Module metrics parameters
    # Run module metrics workflow at the end - on by default
    Boolean? run_module_metrics
    File? primary_contigs_list  # required if run_module_metrics = true

    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    String sv_base_mini_docker
    String sv_base_docker
    String sv_pipeline_base_docker
    String linux_docker

    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_subset_ped
    RuntimeAttr? runtime_attr_sample_list
    RuntimeAttr? runtime_attr_aggregate_tests
    RuntimeAttr? runtime_attr_aggregate_callers
    RuntimeAttr? runtime_attr_petest
    RuntimeAttr? runtime_attr_srtest
    RuntimeAttr? runtime_attr_split_vcf
    RuntimeAttr? runtime_attr_merge_allo
    RuntimeAttr? runtime_attr_merge_stats
    RuntimeAttr? runtime_attr_get_male_only
    RuntimeAttr? runtime_attr_update_metrics
  }

  Array[String] algorithms = ["depth", "melt", "delly", "wham", "manta"]
  Array[File?] vcfs = [depth_vcf, melt_vcf, delly_vcf, wham_vcf, manta_vcf]

  call util.GetSampleIdsFromVcf {
    input:
      vcf = select_first(vcfs),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  call util.SubsetPedFile {
    input:
      ped_file = ped_file,
      sample_list = GetSampleIdsFromVcf.out_file,
      subset_name = batch,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_subset_ped
  }

  call GetSampleLists {
    input:
      ped_file = SubsetPedFile.ped_subset_file,
      samples = GetSampleIdsFromVcf.out_array,
      sv_base_docker = sv_base_docker,
      runtime_attr_override = runtime_attr_sample_list
  }

  scatter (i in range(length(algorithms))) {
    
    if (defined(vcfs[i])) {

      String algorithm = algorithms[i]
      File vcf = select_first([vcfs[i]])

      call GetMaleOnlyVariantIDs {
        input:
          vcf = vcf,
          female_samples = GetSampleLists.female_samples,
          male_samples = GetSampleLists.male_samples,
          contig = select_first([chr_x, "chrX"]),
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_get_male_only
      }

      if (algorithm != "depth") {
        call srt.SRTest as SRTest {
          input:
            splitfile = splitfile,
            medianfile = medianfile,
            ped_file = SubsetPedFile.ped_subset_file,
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
            male_only_variant_ids = GetMaleOnlyVariantIDs.male_only_variant_ids,
            run_common = false,
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
            ped_file = SubsetPedFile.ped_subset_file,
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
            male_only_variant_ids = GetMaleOnlyVariantIDs.male_only_variant_ids,
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

  call UpdateMetrics {
    input:
      old_metrics = original_metrics,
      new_metrics = AggregateCallers.metrics,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_update_metrics
  }

  output {
    File metrics = UpdateMetrics.updated_metrics
  }
}

task UpdateMetrics {
  input {
    File old_metrics
    File new_metrics
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
  String metrics_basename = basename(old_metrics, ".metrics")
  output {
    File updated_metrics = "~{metrics_basename}.updated.metrics"
  }
  command <<<

    set -euo pipefail
    python3 <<CODE
    update_cols = "PE_log_pval PE_called_median PE_bg_median PE_bg_frac SR_posA_log_pval SR_posB_log_pval SR_sum_log_pval SR_posA_called_median SR_posB_called_median SR_sum_called_median SR_posA_bg_median SR_posB_bg_median SR_sum_bg_median SR_posA_bg_frac SR_posB_bg_frac SR_sum_bg_frac SR_posA_pos SR_posB_pos PESR_log_pval PESR_called_median PESR_bg_median PESR_bg_frac".split()
    with open("~{old_metrics}", 'r') as old, open("~{new_metrics}", 'r') as new, open("~{metrics_basename}.updated.metrics", 'w') as out:
      header = None
      for line_old, line_new in zip(old, new):
        if line_old.startswith("name\tchrom\tsvtype"):
          header = line_old.strip().split('\t')
          out.write(line_old) # write header unchanged
          continue
        # non-header lines from here
        fields_old = dict(zip(header, line_old.split('\t')))
        fields_new = dict(zip(header, line_new.split('\t')))
        # check variant ID to ensure files are in same order
        if fields_old["name"] != fields_new["name"]:
          print(f"Variant ID mismatch! Old: {fields_old['name']}, new: {fields_new['name']}")
          exit(1)
        if fields_old["chrom"] == "chrX":
          for col in update_cols:
            fields_old[col] = fields_new[col]
          out.write("\t".join([fields_old[field].strip() for field in header]) + "\n")
        else:
          out.write(line_old) # if not chrX, do not update - write out unchanged
    CODE
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

task GetMaleOnlyVariantIDs {
  input {
    File vcf
    File female_samples
    File male_samples
    String contig
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
    File male_only_variant_ids = "male_only_variant_ids.txt"
  }
  command <<<
    bcftools view -t ~{contig} -S ~{male_samples} ~{vcf} | bcftools view --min-ac 1 | bcftools query -f '%ID\n' > variant_ids_in_males.txt
    bcftools view -t ~{contig} -S ~{female_samples} ~{vcf} | bcftools view --min-ac 1 | bcftools query -f '%ID\n' > variant_ids_in_females.txt
    awk 'NR==FNR{a[$0];next} !($0 in a)' variant_ids_in_females.txt variant_ids_in_males.txt > male_only_variant_ids.txt
    
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
