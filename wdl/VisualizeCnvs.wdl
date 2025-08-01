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
    String flags
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_rdtest
  }

  scatter (file in rd_files) {
    File rd_file_indexes = file + ".tbi"
  }

  call RdTestPlot {
    input:
      vcf_or_bed=vcf_or_bed,
      median_files=median_files,
      ped_file=ped_file,
      rd_files=rd_files,
      rd_file_indexes=rd_file_indexes,
      prefix=prefix,
      min_size=min_size,
      flags=flags,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rdtest
  }
  output{
    File rdtest_plots = RdTestPlot.plots
  }
}

task RdTestPlot {
  input{
    File vcf_or_bed
    Array[File] rd_files
    Array[File] rd_file_indexes
    Array[File] median_files
    Int min_size
    File ped_file
    String prefix
    String sv_pipeline_docker
    String flags
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
    bedtools merge -i cnvs.bed | cut -f1-3 > merged.bed

    i=0
    mkdir rd_subsets
    while read FILE; do
      OUT="rd_subsets/$i.bed.gz"
      tabix -h $FILE -R merged.bed | bgzip > $OUT
      tabix -p bed $OUT
      i=$((i+1))
    done<~{write_lines(rd_files)}

    Rscript /opt/RdTest/RdTestV2.R \
      -b cnvs.bed \
      -n ~{prefix} \
      -x rd_subsets \
      -m median_file.txt \
      -f ~{ped_file} \
      -p TRUE \
      -w samples.txt \
      ~{flags}

    mkdir ~{prefix}_rd_plots
    for file in *.jpg; do
      mv $file ~{prefix}_rd_plots/
    done
    tar -czvf ~{prefix}_rd_plots.tar.gz ~{prefix}_rd_plots/
  >>>

  output {
    File plots = "~{prefix}_rd_plots.tar.gz"
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
