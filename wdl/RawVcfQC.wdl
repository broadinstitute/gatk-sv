## This WDL pipeline implements quality check on vcfs from each SV calling algorithms (https://github.com/arq5x/lumpy-sv)

## two output will be generated from this script: 
##  XXX.low describes individuals with low number of SVs on certain chromosome, indicating FC failure while running  algorithm on the sample
##  XXX.high describes individuals with high number of SVs on certain chromosome,indicating potential library fail on those samples

version 1.0

import "Structs.wdl"

workflow RawVcfQC {
  input {
    Array[File] vcfs
    Array[String] samples
    String prefix
    String caller
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_qc
    RuntimeAttr? runtime_attr_outlier
    RuntimeAttr? runtime_attr_counts
  }

  scatter (i in range(length(vcfs))) {
    call RunIndividualQC {
      input:
        vcf = vcfs[i],
        sample_id = samples[i],
        caller = caller,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_qc
    }
  }

  call PickOutliers {
    input:
      stat_files = RunIndividualQC.stat,
      prefix = prefix,
      caller = caller,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_outlier
  }

  call MergeVariantCounts {
    input:
      stat_files = RunIndividualQC.stat,
      prefix = prefix,
      caller = caller,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_counts
  }

  output {
    File low = PickOutliers.low
    File high = PickOutliers.high
    File variant_counts = MergeVariantCounts.variant_counts
  }
}

task RunIndividualQC {
  input {
    File vcf
    String sample_id
    String caller
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 1, 
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File stat = "${caller}.${sample_id}.QC.stat"
  }
  command <<<
    set -euo pipefail

    python /opt/sv-pipeline/pre_SVCalling_and_QC/raw_vcf_qc/calcu_num_SVs.by_type_chromo.py ~{vcf} ~{caller}.~{sample_id}.QC.stat ~{sample_id}
    
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

task PickOutliers {
  input {
    Array[File] stat_files
    String prefix
    String caller
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(50 + 2 * size(stat_files, "GiB"))

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 7.5,
    disk_gb: disk_gb,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File low = "${prefix}.${caller}.QC.outlier.low"
    File high = "${prefix}.${caller}.QC.outlier.high"
  }
  command <<<
    set -eu -o pipefail
    # concatenate all stat_files into one input file
    xargs cat < ~{write_lines(stat_files)} > ~{caller}.QC.input
    # pick outliers from input stats
    /opt/sv-pipeline/pre_SVCalling_and_QC/raw_vcf_qc/calc_num_svs_pick_outlier.py ~{caller}.QC.input ~{prefix}.~{caller}.QC.outlier -z
    
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

task MergeVariantCounts {
  input {
    Array[File] stat_files
    String prefix
    String caller
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int disk_gb = ceil(50 + 2 * size(stat_files, "GiB"))

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
    File variant_counts = "${prefix}.${caller}.variant_counts.tsv"
  }
  command <<<
    set -euo pipefail
    
    # Create a list of stat files for Python to read
    echo "~{sep='\n' stat_files}" > stat_files.list
    
    python <<CODE
import os

import pandas as pd

dfs = []
with open('stat_files.list', 'r') as file_list:
    for file_path in file_list:
        file_path = file_path.strip()
        if not file_path or not os.path.exists(file_path):
            continue

        df = pd.read_csv(file_path, sep='\t')
        if df.empty:
            continue

        df = df.rename(columns={'SAMPLE': 'sample_id', 'NUM': 'count'})
        df['metric'] = (
            '~{caller}_'
            + df['SVTYPE'].astype(str)
            + '_'
            + df['#CHROM'].astype(str)
        )
        dfs.append(df[['sample_id', 'metric', 'count']])

if dfs:
    combined = pd.concat(dfs, ignore_index=True)
    pivoted = combined.pivot(index='sample_id', columns='metric', values='count')
    pivoted = pivoted.fillna(0).astype(int).reset_index()
else:
    pivoted = pd.DataFrame(columns=['sample_id'])

pivoted.to_csv('~{prefix}.~{caller}.variant_counts.tsv', sep='\t', index=False)
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

