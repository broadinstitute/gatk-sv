version 1.0

import "Structs.wdl"

# Plots CNV depth profiles across batches

workflow VisualizeCnvs {
  input{
    # Note vcf will be faster
    File vcf_or_bed  # bed columns: chrom,start,end,name,svtype,samples

    String prefix
    Array[File] median_files
    Array[File] rd_files
    File ped_file
    Int min_size
    Boolean? restrict_samples
    File? cutoffs
    String flags
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_rdtest
    RuntimeAttr? runtime_attr_rdtestprep
  }

  scatter (file in rd_files) {
    File rd_file_indexes = file + ".tbi"
  }

  call RdTestPrep {
    input:
      vcf_or_bed=vcf_or_bed,
      median_files=median_files,
      rd_files = rd_files,
      rd_file_indexes=rd_file_indexes,
      min_size=min_size,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rdtestprep
  }

  call RdTestPlot {
    input:
      median_file=RdTestPrep.median_file,
      cnvs_bed=RdTestPrep.cnvs_bed,
      samples_in_bed=RdTestPrep.samples_in_bed,
      ped_file=ped_file,
      rd_subsets_tar=RdTestPrep.rd_subsets_tar,
      prefix=prefix,
      cutoffs=cutoffs,
      min_size=min_size,
      restrict_samples=select_first([restrict_samples, true]),
      flags=flags,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rdtest
  }
  output{
    File rdtest_plots = RdTestPlot.plots
    File rdtest_median_geno = RdTestPlot.median_geno
    File rdtest_geno = RdTestPlot.geno
  }
}

task RdTestPrep {
  input{
    File vcf_or_bed
    Array[File] median_files
    Array[File] rd_files
    Array[File] rd_file_indexes
    Int min_size
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(100 + size(vcf_or_bed, "GB") * 10 + size(rd_files, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail

    if [[ ~{vcf_or_bed} == *.vcf.gz ]]; then
      # Subset to DEL/DUP above min size and covert to bed format
      bcftools view -i '(SVTYPE=="DEL" || SVTYPE=="DUP") && SVLEN>=~{min_size}' ~{vcf_or_bed} \
        | svtk vcf2bed stdin raw.bed
      # Swap columns 5/6 for RdTest
      awk -F '\t' -v OFS="\t" '{print $1,$2,$3,$4,$6,$5}' raw.bed > cnvs.bed
      rm raw.bed
      # Get list of all sample IDs
      bcftools query -l ~{vcf_or_bed} > samples.txt
    elif [[ ~{vcf_or_bed} == *.bed || ~{vcf_or_bed} == *.bed.gz ]]; then
      if [[ ~{vcf_or_bed} == *.gz ]]; then
        DECOMPRESSED_BED="raw.bed"
        zcat ~{vcf_or_bed} > $DECOMPRESSED_BED
      else
        DECOMPRESSED_BED="~{vcf_or_bed}"
      fi
      # Subset to DEL/DUP above min size and swap columns 5/6 for RdTest
      awk -F '\t' -v OFS="\t" '{ if ($0!~"#" && $3-$2>=~{min_size} && ($5=="DEL" || $5=="DUP")) {print $1,$2,$3,$4,$6,$5} }' $DECOMPRESSED_BED > cnvs.bed
      # Get list of all sample IDs
      awk -F '\t' -v OFS="\t" '{ if ($0!~"#") {print $6} }' $DECOMPRESSED_BED \
        | sed 's/\,/\n/g' \
        | sort -u \
        > samples.txt
    else
      echo "Invalid extension for input calls. Must be .vcf.gz, .bed.gz, or .bed"
      exit 1
    fi

    paste ~{sep=" " median_files} > median_file.txt
    bgzip median_file.txt


    bedtools merge -i cnvs.bed | cut -f1-3 > merged.bed
    bgzip cnvs.bed

    i=0
    mkdir rd_subsets
    while read FILE; do
      OUT="rd_subsets/$i.bed.gz"
      tabix -h $FILE -R merged.bed | bgzip > $OUT
      tabix -p bed $OUT
      i=$((i+1))
    done<~{write_lines(rd_files)}

    tar -czvf rd_subsets.tar.gz rd_subsets/
  >>>

  output {
    File cnvs_bed = "cnvs.bed.gz"
    File median_file = "median_file.txt.gz"
    File samples_in_bed = "samples.txt.gz"
    File rd_subsets_tar = "rd_subsets.tar.gz"
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

task RdTestPlot {
  input{
    File median_file
    File samples_in_bed
    File rd_subsets_tar
    File cnvs_bed
    Int min_size
    File ped_file
    String prefix
    File? cutoffs
    Boolean restrict_samples
    String sv_pipeline_docker
    String flags
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(100 + size(cnvs_bed, "GB") * 10 + size(rd_subsets_tar, "GB") * 10),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail

    mv ~{median_file} median_file.txt.gz
    gunzip median_file.txt.gz

    mv ~{cnvs_bed} cnvs.bed.gz
    gunzip cnvs.bed.gz

    mv ~{samples_in_bed} samples.txt.gz
    gunzip samples.txt.gz

    tar -xzvf ~{rd_subsets_tar}


    Rscript /opt/RdTest/RdTestV2.R \
      -b cnvs.bed \
      -n ~{prefix} \
      -x rd_subsets \
      -m median_file.txt \
      -f ~{ped_file} \
      -p TRUE \
      ~{if (restrict_samples) then "-w samples.txt" else ""} \
      ~{if defined(cutoffs) then "-g TRUE" else ""} \
      ~{"-r " + cutoffs} \
      ~{flags}

    mkdir ~{prefix}_rd_plots
    mv *jpg ~{prefix}_rd_plots
    tar -czvf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots/
  >>>

  output {
    File plots = "~{prefix}_rd_plots.tar.gz"
    File median_geno = "~{prefix}.median_geno"
    File geno = "~{prefix}.median_geno"
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
