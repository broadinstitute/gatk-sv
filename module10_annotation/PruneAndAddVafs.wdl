# Workflow to perform final sample pruning & compute all relevant AF statistics
# for a VCF from the Talkowski SV pipeline

version 1.0

import "TasksMakeCohortVcf.wdl" as MiniTasks
import "ChromosomeAlleleFrequencies.wdl" as calcAF

# Prune off samples in annotated VCF, add VAF annotation
workflow PruneAndAddVafs {
  
  input {

    File   vcf
    File   vcf_idx
    String prefix
    String contig

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? ped_file                # Used for M/F AF calculations
    File? par_bed
    File? allosomes_list
    File sample_list              # List of samples to be retained from the output vcf

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_combine_vcfs
    RuntimeAttr? runtime_attr_concat_vcfs
    RuntimeAttr? runtime_attr_extract_subset_samples_from_vcf
  }
  
  # Prune VCF
  call ExtractSubsetSamples {
    input:
      vcf        = vcf,
      vcf_idx    = vcf_idx,
      sample_list = sample_list,
      midfix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_extract_subset_samples_from_vcf
  }

  # Compute AC, AN, and AF per population & sex combination
  call calcAF.ChromosomeAlleleFrequencies as ChromosomeAlleleFrequencies {
    input:
      vcf                    = ExtractSubsetSamples.out_vcf,
      vcf_idx                = ExtractSubsetSamples.out_vcf_idx,
      contig                 = contig,
      prefix                 = prefix,
      sample_pop_assignments = sample_pop_assignments,
      ped_file               = ped_file,
      par_bed                = par_bed,
      allosomes_list         = allosomes_list,
      sv_base_mini_docker    = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_shard_vcf    = runtime_attr_shard_vcf,
      runtime_attr_compute_AFs  = runtime_attr_compute_AFs,
      runtime_attr_combine_vcfs = runtime_attr_combine_vcfs
  }

  output {
    File output_vcf     = ChromosomeAlleleFrequencies.vcf_wAFs
    File output_vcf_idx = ChromosomeAlleleFrequencies.vcf_wAFs_idx
  }
}

# Prune off samples from annotated VCF
task PruneVcf {
  
  input {
    File   vcf
    File   vcf_idx
    String contig
    String prefix
    
    File? prune_list

    String sv_base_mini_docker
    
    RuntimeAttr? runtime_attr_override
  }
  
  output {
    File pruned_vcf     = "${prefix}.${contig}.pruned.vcf.gz"
    File pruned_vcf_idx = "${prefix}.${contig}.pruned.vcf.gz.tbi"
  }

  command <<<

    set -euo pipefail
    
    # Tabix chromosome of interest
    tabix -h ~{vcf} ~{contig} | bgzip -c > ~{contig}.vcf.gz
    
    # Get column indexes corresponding to samples to drop, if any exist
    if ~{defined(prune_list)}; then
      dropidx=$( zcat ~{contig}.vcf.gz \
        | sed -n '1,500p' \
        | grep "^#CHROM" \
        | sed 's/\t/\n/g' \
        | awk -v OFS="\t" '{ print NR, $1 }' \
        | fgrep -wf ~{prune_list} \
        | cut -f1 | paste -s -d, )
      zcat ~{contig}.vcf.gz \
        | cut --complement -f "$dropidx" \
        | bgzip -c \
        > "~{prefix}.~{contig}.pruned.vcf.gz"
    else
      cp "~{contig}.vcf.gz" "~{prefix}.~{contig}.pruned.vcf.gz"
    fi
    
    tabix -f "~{prefix}.~{contig}.pruned.vcf.gz"
  
  >>>

  #########################
  RuntimeAttr default_attr = object {
    cpu_cores:          1, 
    mem_gb:             3.75, 
    disk_gb:            250,
    boot_disk_gb:       10,
    preemptible_tries:  3,
    max_retries:        1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])  
  runtime {
    cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
    preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
    docker:                 sv_base_mini_docker
  }
}

task ExtractSubsetSamples {
    input {
        File vcf
        File vcf_idx
        File sample_list
        String midfix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }


    Float input_size = size(vcf, "GB")
    Float base_disk_gb = 10.0
    RuntimeAttr runtime_default = object {
            mem_gb: 3,
            disk_gb: ceil(base_disk_gb + (input_size * 2.0)),
            cpu_cores: 1,
            preemptible_tries: 3,
            max_retries: 1,
            boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
            memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
            disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
            cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
            preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
            maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
            docker: sv_pipeline_docker
            bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, '.vcf.gz')
    command <<<
        set -eu -o pipefail

        bcftools view -S ~{sample_list} ~{vcf} \
        | bgzip > ~{prefix}.~{midfix}.vcf.gz

        tabix -p vcf ~{prefix}.~{midfix}.vcf.gz

    >>>

    output {
        File out_vcf = "~{prefix}.~{midfix}.vcf.gz"
        File out_vcf_idx = "~{prefix}.~{midfix}.vcf.gz.tbi"
    }
}

