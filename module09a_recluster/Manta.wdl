## This WDL pipeline implements SV calling with Illumina's Manta software
##
## Requirements/expectations :
## - Human whole-genome pair-end sequencing data in mapped BAM format

version 1.0

import "Structs.wdl"

workflow Manta {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File region_bed
    File? region_bed_index
    Float? jobs_per_cpu
    Int? mem_gb_per_job
    String manta_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    bam_or_cram_file: ".bam or .cram file to search for SVs. crams are preferable because they localise faster and use less disk."
    bam_or_cram_index: "[optional] Index for bam_or_cram_file. If not specified, index is assumed to be at bam_file_path + '.bai' or cram_file_path + '.crai'"
    sample_id: "sample name. Outputs will be sample_name + 'manta.vcf.gz' and sample_name + 'manta.vcf.gz.tbi'"
    reference_fasta: ".fasta file with reference used to align bam or cram file"
    reference_index: "[optional] If omitted, the WDL will look for an index by appending .fai to the .fasta file"
    region_bed: "[optional] gzipped bed file with included regions where manta should make SV calls."
    region_bed_index: "[optional]If omitted, the WDL will look for an index by appending .tbi to the region_bed file"
    jobs_per_cpu: "[optional] number of manta threads, i.e. num_jobs = round(num_cpu * jobs_per_cpu). If omitted, defaults to 1.3."
    mem_gb_per_job: "[optional] Memory to request for VM (in GB) = round(num_jobs * mem_gb_per_job). If omitted, defaults to 2."
  }
  
  call RunManta {
    input:
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      sample_id = sample_id,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      region_bed = region_bed,
      region_bed_index = region_bed_index,
      jobs_per_cpu = jobs_per_cpu,
      mem_gb_per_job = mem_gb_per_job,
      manta_docker = manta_docker,
      runtime_attr_override = runtime_attr_override
  }

  output {
    File vcf = RunManta.vcf
    File index = RunManta.index
  }
}

task RunManta {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    String sample_id
    File reference_fasta
    File? reference_index
    File region_bed
    File? region_bed_index
    Float? jobs_per_cpu
    Int? mem_gb_per_job
    String manta_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean is_bam = basename(bam_or_cram_file, ".bam") + ".bam" == basename(bam_or_cram_file)
  String bam_ext = if is_bam then ".bam" else ".cram"
  String index_ext = if is_bam then ".bai" else ".crai"
  File ref_index = select_first([reference_index, reference_fasta + ".fai"])
  File bed_index = select_first([region_bed_index, region_bed + ".tbi"])

  # select number of cpus and jobs
  Int num_cpu_use = if defined(runtime_attr_override) then select_first([select_first([runtime_attr_override]).cpu_cores, 8]) else 8
  # select number of jobs (threads) to run
  Float jobs_per_cpu_use = select_first([jobs_per_cpu, 1.3])
  Int num_jobs = round(num_cpu_use * jobs_per_cpu_use)

  # ensure there's sufficient memory.
  Float mem_size_gb = num_jobs * select_first([mem_gb_per_job, 0.4])
  # ALSO: manta will scale down number of jobs if less than 2GB per
  # job are reported, even if that memory is not needed. The memory
  # reported must be an integer
  # Note: unsure if this applies to v1.6.0
  Int mem_reported_to_manta = 2 * num_jobs
  
  # ensure there's sufficient disk space
  Float disk_overhead = 10.0
  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float bam_or_cram_index_size = size(bam_or_cram_index, "GiB")
  Float ref_size = size(reference_fasta, "GiB")
  Float ref_index_size = size(ref_index, "GiB")
  Float region_bed_size = size(region_bed, "GiB")
  Float region_bed_index_size = size(region_bed_index, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + bam_or_cram_index_size + ref_size + ref_index_size + region_bed_size + region_bed_index_size + disk_overhead)

  String expected_index_name = basename(bam_or_cram_file) + index_ext

  RuntimeAttr default_attr = object {
    cpu_cores: num_cpu_use,
    mem_gb: mem_size_gb, 
    disk_gb: vm_disk_size,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File vcf = "${sample_id}.manta.vcf.gz"
    File index = "${sample_id}.manta.vcf.gz.tbi"
  }
  command <<<

    set -Eeuo pipefail

    # if a preemptible instance restarts and runWorkflow.py already
    # exists, manta will throw an error
    if [ -f ./runWorkflow.py ]; then
      rm ./runWorkflow.py
    fi

    ln -s ~{bam_or_cram_file} sample~{bam_ext}
    ln -s ~{bam_or_cram_index} sample~{bam_ext}~{index_ext}


    # prepare the analysis job
    /usr/local/bin/manta/bin/configManta.py \
      --bam sample~{bam_ext} \
      --referenceFasta ~{reference_fasta} \
      --runDir . \
      --callRegions ~{region_bed}

    # always tell manta there are 2 GiB per job, otherwise it will
    # scale back the requested number of jobs, even if they won't
    # need that much memory
    ./runWorkflow.py \
      --mode local \
      --jobs ~{num_jobs} \
      --memGb $((~{num_jobs} * 2))

    # inversion conversion, then compression and index
    python2 /usr/local/bin/manta/libexec/convertInversion.py \
      /usr/local/bin/samtools \
      ~{reference_fasta} \
      results/variants/diploidSV.vcf.gz \
      | bcftools reheader -s <(echo "~{sample_id}") \
      > diploidSV.vcf

    bgzip -c diploidSV.vcf > ~{sample_id}.manta.vcf.gz
    tabix -p vcf ~{sample_id}.manta.vcf.gz
    
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: manta_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}

