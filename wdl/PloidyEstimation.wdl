version 1.0

import "Structs.wdl"

workflow Ploidy {
  input {
    File bincov_matrix
    String batch
    String sv_base_mini_docker
    String sv_pipeline_qc_docker
    RuntimeAttr? runtime_attr_score
    RuntimeAttr? runtime_attr_build
  }

  Int bin_size = 1000000

  call BuildPloidyMatrix {
    input:
      bincov_matrix = bincov_matrix,
      batch = batch,
      bin_size = bin_size,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_build
  }

  call PloidyScore {
    input:
      ploidy_matrix = BuildPloidyMatrix.ploidy_matrix,
      batch = batch,
      sv_pipeline_qc_docker = sv_pipeline_qc_docker,
      runtime_attr_override = runtime_attr_score
  }

  output {
    File ploidy_matrix = BuildPloidyMatrix.ploidy_matrix
    File ploidy_plots = PloidyScore.ploidy_plots
  }
}

task BuildPloidyMatrix {
  input {
    File bincov_matrix
    Int bin_size
    String batch
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File ploidy_matrix = "${batch}_ploidy_matrix.bed.gz"
  }

  command <<<
    set -euo pipefail
    zcat ~{bincov_matrix} \
     | awk ' \
       function printRow() \
         {printf "%s\t%d\t%d",chr,start,stop; \
          for(i=4;i<=nf;++i) {printf "\t%d",vals[i]; vals[i]=0}; \
          print ""} \
       BEGIN {binSize=~{bin_size}} \
       NR==1 {print substr($0,2)} \
       NR==2 {chr=$1; start=$2; stop=start+binSize; nf=NF; for(i=4;i<=nf;++i) {vals[i]=$i}} \
       NR>2  {if($1!=chr){printRow(); chr=$1; start=$2; stop=start+binSize} \
              else if($2>=stop) {printRow(); while($2>=stop) {start=stop; stop=start+binSize}} \
              for(i=4;i<=nf;++i) {vals[i]+=$i}} \
       END   {if(nf!=0)printRow()}' \
     | bgzip > ~{batch}_ploidy_matrix.bed.gz
  >>>

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

task PloidyScore {
  input {
    File ploidy_matrix
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
    File ploidy_plots = "${batch}_ploidy_plots.tar.gz"
  }
  
  command <<<
    set -euo pipefail
    
    mkdir ploidy_est
    Rscript /opt/WGD/bin/estimatePloidy.R -z -O ./ploidy_est ~{ploidy_matrix}

    sleep 10
    
    python /opt/sv-pipeline/02_evidence_assessment/estimated_CN_denoising.py \
      --binwise-copy-number ./ploidy_est/binwise_estimated_copy_numbers.bed.gz \
      --estimated-copy-number ./ploidy_est/estimated_copy_numbers.txt.gz \
      --output-stats cn_denoising_stats.tsv \
      --output-pdf cn_denoising_plots.pdf
    
    cp cn_denoising_stats.tsv ./ploidy_est/
    cp cn_denoising_plots.pdf ./ploidy_est/
    
    tar -zcf ./ploidy_est.tar.gz ./ploidy_est
    mv ploidy_est.tar.gz ~{batch}_ploidy_plots.tar.gz
  >>>
  
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_qc_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

