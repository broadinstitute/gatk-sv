## Author: Ryan L. Collins <rlcollins@g.harvard.edu>
## 
## Workflow to run 00 matrix QC for PE, SR, RD, and BAF in Talkowski SV pipeline

version 1.0

import "Structs.wdl"

workflow MatrixQC {
  input {
    File BAF_file
    Int distance
    String batch
    File PE_file
    File BAF_idx
    File genome_file
    File ref_dict
    File RD_file
    File RD_idx
    File PE_idx
    File SR_idx
    File SR_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_pesrbaf
    RuntimeAttr? runtime_attr_rd
    RuntimeAttr? runtime_attr_plot
  }

  call PESRBAF_QC as BAF_QC {
    input:
      matrix_file = BAF_file,
      matrix_index = BAF_idx,
      genome_file = genome_file,
      ref_dict = ref_dict,
      prefix = "${batch}.BAF",
      ev = "BAF",
      batch = batch,
      distance = distance,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_pesrbaf
  }

  call RD_QC {
    input:
      matrix_file = RD_file,
      matrix_index = RD_idx,
      genome_file = genome_file,
      ref_dict = ref_dict,
      prefix = "${batch}.RD",
      ev = "RD",
      batch = batch,
      distance = distance,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rd
  }

  call PESRBAF_QC as PE_QC {
    input:
      matrix_file = PE_file,
      matrix_index = PE_idx,
      genome_file = genome_file,
      ref_dict = ref_dict,
      prefix = "${batch}.PE",
      ev = "PE",
      batch = batch,
      distance = distance,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_pesrbaf
  }

  call PESRBAF_QC as SR_QC {
    input:
      matrix_file = SR_file,
      matrix_index = SR_idx,
      genome_file = genome_file,
      ref_dict = ref_dict,
      prefix = "${batch}.SR",
      ev = "SR",
      batch = batch,
      distance = distance,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_pesrbaf
  }

  call PlotMatrixQC {
    input:
      bafmatrix = BAF_QC.stats,
      pematrix = PE_QC.stats,
      srmatrix = SR_QC.stats,
      rdmatrix = RD_QC.stats,
      batch = batch,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_plot
  }

  output {
    File PE_stats = PE_QC.stats
    File RD_stats = RD_QC.stats
    File SR_stats = SR_QC.stats
    File BAF_stats = BAF_QC.stats
    File QC_plot=PlotMatrixQC.plot
  }
}

task PESRBAF_QC {
  input {
    File matrix_file
    File matrix_index
    File genome_file
    File ref_dict
    String prefix
    Int distance
    String batch
    String ev
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    matrix_file: {
      localization_optional: true
    }
    matrix_index: {
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
    File stats = "~{batch}.~{ev}.QC_matrix.txt"
  }

  String print_ev_output = "local.~{ev}.txt.gz"

  command <<<

    set -euo pipefail
    fgrep -v "#" ~{genome_file} | awk -v distance=~{distance} -v OFS="\t" '{ print $1, $2-distance, $2 }' > regions.bed

    if [ -s regions.bed ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{matrix_file} \
        -L regions.bed \
        -O ~{print_ev_output}
    else
      touch ~{print_ev_output}
      bgzip ~{print_ev_output}
    fi

    tabix -f -s 1 -b 2 -e 2 ~{print_ev_output}

    /opt/sv-pipeline/00_preprocessing/misc_scripts/nonRD_matrix_QC.sh \
      -d ~{distance} \
      ~{print_ev_output} \
      ~{genome_file} \
      ~{batch}.~{ev}.QC_stats.txt
    cut -f1 ~{genome_file} > contigs.list
    python /opt/sv-pipeline/00_preprocessing/misc_scripts/qcstat2matrix.py ~{batch}.~{ev}.QC_stats.txt ~{batch} ~{ev} contigs.list
  
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

task RD_QC {
  input {
    File matrix_file
    File matrix_index
    File genome_file
    File ref_dict
    String prefix
    String ev
    String batch
    Int distance
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    matrix_file: {
      localization_optional: true
    }
    matrix_index: {
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
    File stats = "~{batch}.~{ev}.QC_matrix.txt"
  }

  String print_ev_output = "local.~{ev}.txt.gz"

  command <<<

    set -euo pipefail
    fgrep -v "#" ~{genome_file} | awk -v distance=~{distance} -v OFS="\t" '{ print $1, $2-distance, $2 }' > regions.bed

    if [ -s regions.bed ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{matrix_file} \
        -L regions.bed \
        -O ~{print_ev_output}
    else
      touch local.RD.txt
      bgzip local.RD.txt
    fi

    tabix -f -p bed ~{print_ev_output}

    /opt/sv-pipeline/00_preprocessing/misc_scripts/RD_matrix_QC.sh \
      -d ~{distance} \
      ~{print_ev_output} \
      ~{genome_file} \
      ~{batch}.~{ev}.QC_stats.txt
    cut -f1 ~{genome_file} > contigs.list
    python /opt/sv-pipeline/00_preprocessing/misc_scripts/qcstat2matrix.py ~{batch}.~{ev}.QC_stats.txt ~{batch} ~{ev} contigs.list
  
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

# Combine and plot
task PlotMatrixQC {
  input {
    File bafmatrix
    File pematrix
    File srmatrix
    File rdmatrix
    String batch
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

  command <<<
    cp ~{bafmatrix} . 
    cp ~{pematrix} .
    cp ~{srmatrix} .
    cp ~{rdmatrix} .
    Rscript /opt/sv-pipeline/00_preprocessing/misc_scripts/plot_00_matrix_FC_QC.R ~{batch}
  >>>
  
  output {
    File plot = "~{batch}.00_matrix_FC_QC.png"
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

