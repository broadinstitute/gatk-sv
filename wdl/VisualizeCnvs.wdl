version 1.0

import "Structs.wdl"

# Plots CNV depth profiles across batches

workflow VisualizeCnvs {
  input {
    # Note vcf will be faster
    File vcf_or_bed  # bed columns: chrom,start,end,name,svtype,samples
    String prefix
    Array[File] median_files
    Array[File] rd_files
    Array[String] batches  # in same order as median_files, rd_files
    File samples_in_batches  # batch ID \t sample ID
    File ped_file
    Int min_size
    Int max_size
    Int records_per_shard
    String flags
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_rdtest
  }

  scatter (file in rd_files) {
    File rd_file_indexes = file + ".tbi"
  }

  call RdTestScatter {
    input:
      vcf_or_bed=vcf_or_bed,
      shard_prefix="~{prefix}.shard",
      min_size=min_size,
      lines_per_shard=records_per_shard,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rdtest
  }

  scatter (shard in RdTestScatter.shards) {
    call RdTestPlot {
      input:
        bed_shard=shard,
        median_files=write_lines(median_files),
        ped_file=ped_file,
        rd_files=write_lines(rd_files),
        rd_file_indexes=write_lines(rd_file_indexes),
        batches = write_lines(batches),
        samples_in_batches=samples_in_batches,
        prefix=prefix,
        flags=flags,
        max_size=max_size,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override = runtime_attr_rdtest
    }
  }

  call TarRdPlots {
    input:
      raw_plots=flatten(RdTestPlot.plots),
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override = runtime_attr_rdtest
  }

  output {
    File rdtest_plots = TarRdPlots.plots
  }
}

task RdTestScatter {
  input {
    File vcf_or_bed
    Int lines_per_shard
    Int min_size
    String shard_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(100 + size(vcf_or_bed, "GB") * 10),
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
      awk -F '\t' -v OFS="\t" '{ if ($0!~"#") {print $1,$2,$3,$4,$6,$5} }' raw.bed > cnvs.bed
      rm raw.bed
    elif [[ ~{vcf_or_bed} == *.bed || ~{vcf_or_bed} == *.bed.gz ]]; then
      if [[ ~{vcf_or_bed} == *.gz ]]; then
        DECOMPRESSED_BED="raw.bed"
        zcat ~{vcf_or_bed} > $DECOMPRESSED_BED
      else
        DECOMPRESSED_BED="~{vcf_or_bed}"
      fi
      # Subset to DEL/DUP above min size and swap columns 5/6 for RdTest
      awk -F '\t' -v OFS="\t" '{ if ($0!~"#" && $3-$2>=~{min_size} && ($5=="DEL" || $5=="DUP")) {print $1,$2,$3,$4,$6,$5} }' $DECOMPRESSED_BED > cnvs.bed
    else
      echo "Invalid extension for input calls. Must be .vcf.gz, .bed.gz, or .bed"
      exit 1
    fi

    split -d -a 6 -l ~{lines_per_shard} \
      --additional-suffix ".bed" \
      cnvs.bed \
      ~{shard_prefix}
  >>>

  output {
     Array[File] shards=glob("~{shard_prefix}*")
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
  input {
    File bed_shard
    File rd_files
    File rd_file_indexes
    File median_files
    File batches
    Int max_size
    File ped_file
    File samples_in_batches
    String prefix
    String sv_pipeline_docker
    String flags
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 7.5,
                               disk_gb: ceil(200 + size(bed_shard, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    # for each record in shard, find relevant samples & batches and subset the data files appropriately, then run RdTest
    while read RECORD; do
      python <<CODE
import random
import math

chrom, start, end, name, samples, svtype = "$RECORD".strip().split('\t')
with open("cnvs.bed", 'w') as out:
    # segment large CNVs into multiple records
    chunk_size = int(~{max_size})
    start, end = int(start), int(end)
    svlen = end - start
    chunks = math.ceil(svlen / chunk_size)
    if chunks > 1:
        for i in range(chunks):
            chunk_start = str(int(start + i*chunk_size))
            chunk_end = str(int(min(start + (i+1)*chunk_size, end)))
            chunk_name = name + "." + str(i)
            out.write("\t".join([chrom, chunk_start, chunk_end, chunk_name, samples, svtype]) + "\n")
    else:
        out.write("$RECORD" + "\n")

# get batches containing record
batches_record = set()
with open("~{samples_in_batches}", 'r') as samples_batches:
    samples = set(samples.split(","))
    for line in samples_batches:
        batch, sample = line.strip().split("\t")
        if sample in samples:
            batches_record.add(batch)
    # randomly choose 10 batches if there are more
    if len(batches_record) > 10:
        batches_record = random.sample(batches_record, 10)

with open("~{samples_in_batches}", 'r') as samples_batches, open("samples.txt", 'w') as samples_file:
    # get sample IDs for record and write to file: all samples in batches containing record
    for line in samples_batches:
        batch, sample = line.strip().split("\t")
        if batch in batches_record:
            samples_file.write(f"{sample}\n")

# get median files, rd files, and rd indexes for chosen batches
with open("median_files.txt", 'w') as med, open("rd_files.txt", 'w') as rdf, open("rd_indexes.txt", 'w') as rdi, \
    open("~{median_files}", 'r') as med_list, open("~{rd_files}", 'r') as rdf_list, open("~{rd_file_indexes}", 'r') as rdi_list, \
    open("~{batches}", 'r') as batches_list:
    for batch, med_file, rd_file, rd_index in zip(batches_list, med_list, rdf_list, rdi_list):
        if batch.strip() in batches_record:
            med.write(med_file)
            rdf.write(rd_file)
            rdi.write(rd_index)
CODE

      # localize and merge relevant batch median files
      # delete any local dirs if they exist from the previous record
      if [ -d "med_files" ]; then rm -r med_files; fi
      mkdir med_files
      cat median_files.txt | gsutil -m cp -I med_files/
      paste med_files/* > median_file.txt

      # localize and tabix region of interest from relevant batch RD files
      echo $RECORD > merged.bed  # unsegmented record
      i=0
      if [ -d "rd_files" ]; then rm -r rd_files; fi
      mkdir rd_files
      cat rd_files.txt | gsutil -m cp -I rd_files/
      cat rd_indexes.txt | gsutil -m cp -I rd_files/
      if [ -d "rd_subsets" ]; then rm -r rd_subsets; fi
      mkdir rd_subsets
      while read FILE; do
        local_file="rd_files/"$(basename $FILE)
        OUT="rd_subsets/$i.bed.gz"
        tabix -h $local_file -R merged.bed | bgzip > $OUT
        tabix -p bed $OUT
        i=$((i+1))
      done<rd_files.txt

      # run RdTest on segmented CNVs with appropriate sample/batch/data inputs
      Rscript /opt/RdTest/RdTestV2.R \
        -b cnvs.bed \
        -n ~{prefix} \
        -x rd_subsets \
        -m median_file.txt \
        -f ~{ped_file} \
        -p TRUE \
        -w samples.txt \
        ~{flags}

    done < ~{bed_shard}
  >>>

  output {
    Array[File] plots = glob("*jpg")
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


task TarRdPlots {
  input {
    Array[File] raw_plots
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(100 + size(raw_plots, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail

    mkdir ~{prefix}_rd_plots
    while read plot; do
      mv "$plot" ~{prefix}_rd_plots
    done < ~{write_lines(raw_plots)}
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
