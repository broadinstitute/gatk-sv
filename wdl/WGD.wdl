version 1.0

import "Structs.wdl"

workflow WGD {
  input {
    String batch
    File wgd_scoring_mask
    File bincov_matrix
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_build
    RuntimeAttr? runtime_attr_score
  }

  call BuildWGDMatrix {
    input:
      bincov_matrix = bincov_matrix,
      wgd_scoring_mask = wgd_scoring_mask,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      batch = batch
  }

  call WGDScore {
    input:
      wgd_scoring_mask = wgd_scoring_mask,
      WGD_matrix = BuildWGDMatrix.WGD_matrix,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      batch = batch
  }

  output {
    File WGD_dist = WGDScore.WGD_dist
    File WGD_matrix = BuildWGDMatrix.WGD_matrix
    File WGD_scores = WGDScore.WGD_scores
  }
}

task BuildWGDMatrix {
  input {
    File bincov_matrix
    File wgd_scoring_mask
    String batch
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File WGD_matrix = "${batch}_WGD_scoring_matrix_output.bed.gz"
  }
  command <<<

    set -eu
    zcat ~{bincov_matrix} | head -n 1 > header.txt
    sed -i 's/#//g' header.txt
    
    set -o pipefail
    zcat ~{bincov_matrix} \
    | bedtools intersect -f 0.49 -wa -u \
      -a - \
      -b ~{wgd_scoring_mask} \
    | sort -Vk1,1 -k2,2n -k3,3n \
    > ~{batch}_WGD_scoring_matrix.bed
    
    cat header.txt ~{batch}_WGD_scoring_matrix.bed \
    | bgzip -c \
    > ~{batch}_WGD_scoring_matrix_output.bed.gz

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_qc_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

task WGDScore {
  input {
    File wgd_scoring_mask
    File WGD_matrix
    String batch
    String sv_pipeline_qc_docker
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
    File WGD_scores = "${batch}_WGD_scores.txt.gz"
    File WGD_dist = "${batch}_WGD_score_distributions.pdf"
  }
  command <<<

    set -euo pipefail
    Rscript /opt/WGD/bin/scoreDosageBiases.R -z -O . ~{WGD_matrix} ~{wgd_scoring_mask}
    mv WGD_scores.txt.gz ~{batch}_WGD_scores.txt.gz
    mv WGD_score_distributions.pdf ~{batch}_WGD_score_distributions.pdf
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_qc_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

