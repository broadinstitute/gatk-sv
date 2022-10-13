version 1.0

import "Structs.wdl"
import "TasksBenchmark.wdl" as tasks10

workflow VaPoRBed {
  input {
    String prefix
    String bam_or_cram_file
    String bam_or_cram_index
    File bed_file
    String? sample_to_extract
    File ref_fasta
    File ref_fai
    File ref_dict
    Array[String] contigs
    String vapor_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_Vapor 
    RuntimeAttr? runtime_attr_bcf2vcf
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_SplitVcf
    RuntimeAttr? runtime_attr_ConcatBeds
    RuntimeAttr? runtime_attr_LocalizeCram
  }

  scatter ( contig in contigs ) {

    call PreprocessBedForVaPoR {
      input:
        contig = contig,
        sample_to_extract = sample_to_extract,
        bed_file = bed_file,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override=runtime_attr_SplitVcf
    }
    
    call RunVaPoRWithCram as RunVaPoR {
      input:
        prefix = prefix,
        contig = contig,
        bam_or_cram_file=bam_or_cram_file,
        bam_or_cram_index=bam_or_cram_index,
        bed = PreprocessBedForVaPoR.contig_bed,
        ref_fasta = ref_fasta,
        ref_fai = ref_fai,
        ref_dict = ref_dict,
        vapor_docker = vapor_docker,
        runtime_attr_override = runtime_attr_Vapor
    }
  }

  call tasks10.ConcatVaPoR {
    input:
      shard_bed_files=RunVaPoR.vapor,
      shard_plots = RunVaPoR.vapor_plot,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_ConcatBeds
  }

  output {
      File bed = ConcatVaPoR.merged_bed_file
      File plots = ConcatVaPoR.merged_bed_plot
    }
  }

# extract specific contig from BED, and sites for sample if provided, and add SVLEN to INS if header contains SVLEN column
task PreprocessBedForVaPoR {
  input {
    String contig
    String? sample_to_extract
    File bed_file  # first 5 columns must be chrom, start, end, name, svtype (or VaPoR description). if >5 columns, use header or assume samples is 6th. Need header & SVLEN column unless already appended to INS descriptions
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File contig_bed = "~{contig}.vapor.bed"
  }

  command <<<
    python3 <<CODE
    import gzip

    def is_gzipped(path):
      return path.endswith(".gz")


    open_fn = gzip.open if is_gzipped("~{bed_file}") else open
    open_mode = 'rt' if is_gzipped("~{bed_file}") else 'r'

    remove_types = {"BND", "CPX", "CNV"}

    default_columns = "chrom start end name svtype".split()
    default_num_columns = len(default_columns)
    columns = dict(zip(default_columns, range(default_num_columns)))
    first = True
    with open_fn("~{bed_file}", open_mode) as inp, open("~{contig}.vapor.bed", 'w') as out:
      for line in inp:
        fields = line.lstrip("#").rstrip('\n').split('\t')
        if first:
          if line.startswith("#"):
            for i, name in enumerate(fields[default_num_columns:]):
              columns[name] = default_num_columns + i  # get column names beyond first default ones from header if available
          else:
            if len(fields) >= default_num_columns:
              columns['samples'] = default_num_columns  # if no header but extra fields, assume samples is next column
        if fields[columns["chrom"]] != "~{contig}":
          continue  # extract only contig of interest. also drops header if exists
        if fields[columns["svtype"]] in remove_types:
          continue  # drop BND, CNV, CPX
        if "~{sample_to_extract}" != "" and "samples" in columns and "~{sample_to_extract}" not in fields[columns["samples"]]:
          continue  # extract events in sample of interest if provided
        svtype_write = fields[columns["svtype"]]
        if "INS" in svtype_write or "MEI" in svtype_write:
          if "SVLEN" in columns:
            svtype_write = f"INS_{fields[columns['SVLEN']]}"  # for INS, format svtype as INS_SVLEN for VaPoR if SVLEN column in input BED
        # write chrom, pos, end, SVID, svtype/description for vapor
        out.write("\t".join([fields[columns["chrom"]], fields[columns["start"]], fields[columns["end"]], fields[columns["name"]], svtype_write]) + "\n")
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

task RunVaPoR {
  input {
    String prefix
    String contig
    File bam_or_cram_file
    File bam_or_cram_index
    File bed
    File ref_fasta
    File ref_fai
    File ref_dict
    String vapor_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vapor = "~{prefix}.~{contig}.vapor.gz"
    File vapor_plot = "~{prefix}.~{contig}.tar.gz"
  }

  command <<<

    set -Eeuo pipefail

    mkdir ~{prefix}.~{contig}
  
    vapor bed \
    --sv-input ~{bed} \
    --output-path ~{prefix}.~{contig} \
    --output-file ~{prefix}.~{contig}.vapor \
    --reference ~{ref_fasta} \
    --PB-supp 0 \
    --pacbio-input ~{bam_or_cram_file}

    tar -czf ~{prefix}.~{contig}.tar.gz ~{prefix}.~{contig}
    bgzip  ~{prefix}.~{contig}.vapor
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: vapor_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task RunVaPoRWithCram {
  input {
    String prefix
    String contig
    String bam_or_cram_file
    String bam_or_cram_index
    File bed
    File ref_fasta
    File ref_fai
    File ref_dict
    String vapor_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 15, 
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  output {
    #File local_bam = "~{contig}.bam"
    File vapor = "~{prefix}.~{contig}.vapor.gz"
    File vapor_plot = "~{prefix}.~{contig}.tar.gz"
  }

  command <<<

    set -Eeuo pipefail

    #localize cram files
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`   
    samtools view -h -o ~{contig}.bam ~{bam_or_cram_file} ~{contig}
    samtools index ~{contig}.bam
  
    #run vapor
    mkdir ~{prefix}.~{contig}

    vapor bed \
    --sv-input ~{bed} \
    --output-path ~{prefix}.~{contig} \
    --output-file ~{prefix}.~{contig}.vapor \
    --reference ~{ref_fasta} \
    --PB-supp 0 \
    --pacbio-input ~{contig}.bam

    tar -czf ~{prefix}.~{contig}.tar.gz ~{prefix}.~{contig}
    bgzip  ~{prefix}.~{contig}.vapor
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: vapor_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
