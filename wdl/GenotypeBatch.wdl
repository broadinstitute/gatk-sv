version 1.0

import "TasksGenotypeBatch.wdl" as tasks_genotype_batch
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Structs.wdl"

workflow GenotypeBatch {
  input {
    String batch
    File vcf

    File training_intervals
    File median_coverage
    File rd_file
    File pe_file
    File sr_file
    File reference_dict
    File ploidy_table
    File depth_exclusion_intervals
    File pesr_exclusion_intervals
    File rf_cutoffs

    File contig_list

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_attr_train
    RuntimeAttr? runtime_attr_genotype
    RuntimeAttr? runtime_override_concat_shards
    RuntimeAttr? runtime_attr_regeno_coverage_medians
  }

  call TrainSVGenotyping {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      output_name = batch,
      training_intervals = training_intervals,
      median_coverage = median_coverage,
      rd_file = rd_file,
      rd_file_index = rd_file + ".tbi",
      pe_file = pe_file,
      pe_file_index = pe_file + ".tbi",
      sr_file = sr_file,
      sr_file_index = sr_file + ".tbi",
      reference_dict = reference_dict,
      ploidy_table = ploidy_table,
      depth_exclusion_intervals = depth_exclusion_intervals,
      depth_exclusion_intervals_index = depth_exclusion_intervals + ".tbi",
      pesr_exclusion_intervals = pesr_exclusion_intervals,
      pesr_exclusion_intervals_index = pesr_exclusion_intervals + ".tbi",
      rf_cutoffs = rf_cutoffs,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_train
  }


  scatter (contig in read_lines(contig_list)) {
    call GenotypeSVs {
      input:
        vcf = vcf,
        vcf_index = vcf + ".tbi",
        output_prefix = "~{batch}.genotype_batch.~{contig}",
        contig = contig,
        median_coverage = median_coverage,
        rd_file = rd_file,
        rd_file_index = rd_file + ".tbi",
        pe_file = pe_file,
        pe_file_index = pe_file + ".tbi",
        sr_file = sr_file,
        sr_file_index = sr_file + ".tbi",
        reference_dict = reference_dict,
        ploidy_table = ploidy_table,
        depth_exclusion_intervals = depth_exclusion_intervals,
        depth_exclusion_intervals_index = depth_exclusion_intervals + ".tbi",
        pesr_exclusion_intervals = pesr_exclusion_intervals,
        pesr_exclusion_intervals_index = pesr_exclusion_intervals + ".tbi",
        rd_table = TrainSVGenotyping.rd_table,
        pe_table = TrainSVGenotyping.pe_table,
        sr_table = TrainSVGenotyping.sr_table,
        gatk_docker = gatk_docker,
        runtime_attr_override = runtime_attr_genotype
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs = GenotypeSVs.out,
      vcfs_idx = GenotypeSVs.out_index,
      naive = true,
      outfile_prefix = batch + ".genotype_batch",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat_shards
  }

  call SeparateDepthPesr {
    input:
      vcf = ConcatVcfs.concat_vcf,
      vcf_index = ConcatVcfs.concat_vcf_idx,
      prefix = batch + ".genotype_batch",
      sv_base_mini_docker = sv_base_mini_docker
  }

  call GenerateRegenoCoverageMedians {
    input:
      vcf = SeparateDepthPesr.depth_vcf,
      prefix = "~{batch}.regeno_coverage_medians",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_regeno_coverage_medians
  }

  output {
    File genotyped_depth_vcf = SeparateDepthPesr.depth_vcf
    File genotyped_depth_vcf_index = SeparateDepthPesr.depth_vcf_index
    File genotyped_pesr_vcf = SeparateDepthPesr.pesr_vcf
    File genotyped_pesr_vcf_index = SeparateDepthPesr.pesr_vcf_index
    File genotyping_rd_table = TrainSVGenotyping.rd_table
    File genotyping_pe_table = TrainSVGenotyping.pe_table
    File genotyping_sr_table = TrainSVGenotyping.sr_table
    File regeno_coverage_medians = GenerateRegenoCoverageMedians.out
  }
}

task TrainSVGenotyping {
  input {
    File vcf
    File vcf_index
    File training_intervals
    File median_coverage
    File rd_file
    File rd_file_index
    File pe_file
    File pe_file_index
    File sr_file
    File sr_file_index
    File reference_dict
    File ploidy_table
    File pesr_exclusion_intervals
    File pesr_exclusion_intervals_index
    File depth_exclusion_intervals
    File depth_exclusion_intervals_index
    String output_name

    File rf_cutoffs

    String gatk_docker
    Float? java_mem_fraction
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 16,
                               disk_gb: ceil(50 + size([vcf, rd_file, pe_file, sr_file], "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    PEQ=$(awk -F '\t' '{if ( $5=="PEQ") print $2 }' ~{rf_cutoffs})
    SRQ=$(awk -F '\t' '{if ( $5=="SRQ") print $2 }' ~{rf_cutoffs})

    function getJavaMem() {
      cat /proc/meminfo | awk -v MEM_FIELD="$1" '{
        f[substr($1, 1, length($1)-1)] = $2
      } END {
        printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
      }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)

    gatk --java-options "-Xmx${JVM_MAX_MEM}" TrainSVGenotyping \
      -XL chrX -XL chrY \
      -V ~{vcf} \
      --training-intervals ~{training_intervals} \
      -O ~{output_name}.vcf.gz \
      --median-coverage ~{median_coverage} \
      --rd-file ~{rd_file} \
      --split-reads-file ~{sr_file} \
      --discordant-pairs-file ~{pe_file} \
      --sequence-dictionary ~{reference_dict} \
      --ploidy-table ~{ploidy_table} \
      --depth-exclusion-intervals ~{depth_exclusion_intervals} \
      --pesr-exclusion-intervals ~{pesr_exclusion_intervals} \
      --pe-quality ${PEQ} \
      --sr-quality ${SRQ} \
      --output-dir ./ \
      --output-name ~{output_name}
  >>>

  output {
    File rd_table = "~{output_name}.rd_geno_params.tsv"
    File pe_table = "~{output_name}.pe_geno_params.tsv"
    File sr_table = "~{output_name}.sr_geno_params.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task GenotypeSVs {
  input {
    File vcf
    File vcf_index
    String output_prefix
    String? contig
    File median_coverage
    File rd_file
    File rd_file_index
    File pe_file
    File pe_file_index
    File sr_file
    File sr_file_index
    File reference_dict
    File ploidy_table
    File pesr_exclusion_intervals
    File pesr_exclusion_intervals_index
    File depth_exclusion_intervals
    File depth_exclusion_intervals_index
    File rd_table
    File pe_table
    File sr_table

    String gatk_docker
    Float? java_mem_fraction
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    rd_file: {
                localization_optional: true
              }
    pe_file: {
                localization_optional: true
              }
    sr_file: {
                localization_optional: true
             }
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size([vcf, rd_file, pe_file, sr_file], "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    function getJavaMem() {
      cat /proc/meminfo | awk -v MEM_FIELD="$1" '{
        f[substr($1, 1, length($1)-1)] = $2
      } END {
        printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
      }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)

    gatk --java-options "-Xmx${JVM_MAX_MEM}" PrintSVEvidence \
      --sequence-dictionary ~{reference_dict} \
      --evidence-file ~{rd_file} \
      ~{"-L " + contig} \
      -O local.rd.txt.gz

    gatk --java-options "-Xmx${JVM_MAX_MEM}" PrintSVEvidence \
      --sequence-dictionary ~{reference_dict} \
      --evidence-file ~{pe_file} \
      ~{"-L " + contig} \
      -O local.pe.txt.gz

    gatk --java-options "-Xmx${JVM_MAX_MEM}" PrintSVEvidence \
      --sequence-dictionary ~{reference_dict} \
      --evidence-file ~{sr_file} \
      ~{"-L " + contig} \
      -O local.sr.txt.gz

    gatk --java-options "-Xmx${JVM_MAX_MEM}" GenotypeSVs \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      ~{"-L " + contig} \
      --median-coverage ~{median_coverage} \
      --rd-file local.rd.txt.gz \
      --discordant-pairs-file local.pe.txt.gz \
      --split-reads-file local.sr.txt.gz \
      --sequence-dictionary ~{reference_dict} \
      --ploidy-table ~{ploidy_table} \
      --pesr-exclusion-intervals ~{pesr_exclusion_intervals} \
      --depth-exclusion-intervals ~{depth_exclusion_intervals} \
      --rd-table ~{rd_table} \
      --pe-table ~{pe_table} \
      --sr-table ~{sr_table}
  >>>

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Separate depth-only and pesr-only records from in VCF
task SeparateDepthPesr {
  input {
    File vcf
    File vcf_index
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(vcf, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euxo pipefail
    bcftools view -i 'INFO/ALGORITHMS=="depth"' ~{vcf} -Oz -o ~{prefix}.depth.vcf.gz
    tabix ~{prefix}.depth.vcf.gz
    bcftools view -i 'INFO/ALGORITHMS!="depth"' ~{vcf} -Oz -o ~{prefix}.pesr.vcf.gz
    tabix ~{prefix}.pesr.vcf.gz
  >>>

  output {
    File depth_vcf = "~{prefix}.depth.vcf.gz"
    File depth_vcf_index = "~{prefix}.depth.vcf.gz.tbi"
    File pesr_vcf = "~{prefix}.pesr.vcf.gz"
    File pesr_vcf_index = "~{prefix}.pesr.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

  task GenerateRegenoCoverageMedians {
    input {
        File vcf
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
                                 cpu_cores: 1,
                                 mem_gb: 3.75,
                                 disk_gb: ceil(50 + size(vcf, "GB")),
                                 boot_disk_gb: 10,
                                 preemptible_tries: 3,
                                 max_retries: 1
                               }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    command <<<
      set -euxo pipefail
      python /opt/sv-pipeline/scripts/extract_format_table.py \
        --format-field "RD_MCR" \
        --id-column "cnvID" \
        --vcf ~{vcf} \
        --out ~{prefix}.tsv
      gzip ~{prefix}.tsv
    >>>

    output {
        File out = "~{prefix}.tsv.gz"
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