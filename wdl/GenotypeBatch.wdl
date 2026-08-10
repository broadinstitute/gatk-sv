version 1.0

import "TasksGenotypeBatch.wdl" as tasks_genotype_batch
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Structs.wdl"

workflow GenotypeBatch {
  input {
    String batch
    # TODO: this vcf is sites-only and all samples in rd_file/pe_file/sr_file will be genotyped at these sites,
    # but we should add an option to subset to just retained samples in case of sample filtering
    File vcf

    # Sites used to train the PE/SR genotyping models: normally this batch's own FilterBatch PESR
    # calls. Deliberately separate from `vcf`, which is the (potentially cohort-wide) set of sites
    # to genotype.
    #
    # These should not be the same file when `vcf` is cohort-wide. The SR frequency cutoff optimizer
    # learns recovery-fraction thresholds by contrasting variants that have depth/PE support against
    # those that do not. Cohort-wide sites absent from this batch contribute background-only
    # fractions to both classes, which erases the separation the optimizer needs. On one AoU batch
    # that drove the pass:fail ratio from ~1:10 to ~1:26 and collapsed all four cutoffs to zero,
    # silently disabling SR background filtering.
    #
    # Only PE and SR training read this VCF. RD is trained from `training_intervals` plus the depth
    # matrix, and its depth-only vs PESR split comes from separation constraints rather than from
    # the VCF, so a PESR-only training VCF is sufficient here.
    File training_vcf

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

    String? training_args
    String? genotype_args

    # Fail the workflow if SR frequency cutoff training produced a degenerate grid, which
    # silently resolves to zero cutoffs and disables SR background filtering. Set false to let
    # a diagnostic run complete and inspect genotyping_sr_cutoff_diagnostics anyway.
    Boolean fail_on_degenerate_sr_cutoffs = true

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_format
    RuntimeAttr? runtime_attr_train
    RuntimeAttr? runtime_attr_genotype
    RuntimeAttr? runtime_override_concat_shards
    RuntimeAttr? runtime_attr_regeno_coverage_medians
  }

  call FilterWhamDeletions {
    input:
      vcf = vcf,
      vcf_index = vcf + ".tbi",
      prefix = batch + ".wham_del_filtered",
      sv_base_mini_docker = sv_base_mini_docker
  }

  # Apply the same wham-deletion exclusion to the training sites. Training on variant classes that
  # are then filtered out of the genotyped output would calibrate the cutoffs against records that
  # never get genotyped.
  call FilterWhamDeletions as FilterWhamDeletionsTraining {
    input:
      vcf = training_vcf,
      vcf_index = training_vcf + ".tbi",
      prefix = batch + ".training.wham_del_filtered",
      sv_base_mini_docker = sv_base_mini_docker
  }

  call TrainSVGenotyping {
    input:
      vcf = FilterWhamDeletionsTraining.filtered_vcf,
      vcf_index = FilterWhamDeletionsTraining.filtered_vcf_index,
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
      training_args = training_args,
      gatk_docker = gatk_docker,
      runtime_attr_override = runtime_attr_train
  }

  # TrainSVGenotyping cannot fail on a degenerate SR cutoff grid: Cromwell delocalizes task
  # outputs only on success, so failing there would strand the diagnostics needed to explain it.
  # Enforcement therefore happens here, after the diagnostics have been copied out. Gates the
  # genotyping scatter so a doomed run does not pay for every contig shard.
  call ValidateSRCutoffs {
    input:
      sr_table = TrainSVGenotyping.sr_table,
      sr_cutoff_diagnostics = TrainSVGenotyping.sr_cutoff_diagnostics,
      fail_on_degenerate_sr_cutoffs = fail_on_degenerate_sr_cutoffs,
      sv_base_mini_docker = sv_base_mini_docker
  }

  scatter (contig in read_lines(contig_list)) {
    call GenotypeSVs {
      input:
        vcf = FilterWhamDeletions.filtered_vcf,
        vcf_index = FilterWhamDeletions.filtered_vcf_index,
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
        rd_depth_table = TrainSVGenotyping.rd_depth_table,
        rd_pesr_table = TrainSVGenotyping.rd_pesr_table,
        pe_table = TrainSVGenotyping.pe_table,
        # Via ValidateSRCutoffs so the scatter cannot start on a rejected cutoff grid
        sr_table = ValidateSRCutoffs.validated_sr_table,
        genotype_args = genotype_args,
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
    File genotyping_rd_depth_table = TrainSVGenotyping.rd_depth_table
    File genotyping_rd_pesr_table = TrainSVGenotyping.rd_pesr_table
    File genotyping_pe_table = TrainSVGenotyping.pe_table
    File genotyping_sr_table = TrainSVGenotyping.sr_table
    File genotyping_sr_cutoff_diagnostics = TrainSVGenotyping.sr_cutoff_diagnostics
    File regeno_coverage_medians = GenerateRegenoCoverageMedians.out
  }
}

task FilterWhamDeletions {
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
    set -euo pipefail

    bcftools view \
      -e 'INFO/SVTYPE="DEL" && COUNT(INFO/ALGORITHMS)=1 && INFO/ALGORITHMS[0]="wham"' \
      ~{vcf} \
      -Oz \
      -o ~{prefix}.vcf.gz

    tabix -p vcf ~{prefix}.vcf.gz
  >>>

  output {
    File filtered_vcf = "~{prefix}.vcf.gz"
    File filtered_vcf_index = "~{prefix}.vcf.gz.tbi"
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
    String? training_args
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

    PEQ=$(awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
               {if ($col["metric"]=="PEQ") print $col["cutoff"]}' \
          ~{rf_cutoffs} | head -n 1)
    SRQ=$(awk -F '\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
               {if ($col["metric"]=="SRQ") print $col["cutoff"]}' \
          ~{rf_cutoffs} | head -n 1)

    PESR_SEP=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
                 {
                   if (toupper($col["algtype"])=="PESR" &&
                     ($col["min_svsize"] + 0)==1000 &&
                     toupper($col["metric"])=="RD_MEDIAN_SEPARATION") {
                     print $col["cutoff"]
                   }
                 }' \
            ~{rf_cutoffs} | sort -nr | head -n 1)
    DEPTH_SEP=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
                {
                  if (toupper($col["algtype"])=="DEPTH" &&
                    toupper($col["metric"])=="RD_MEDIAN_SEPARATION") {
                    print $col["cutoff"]
                  }
                }' \
             ~{rf_cutoffs} | sort -nr | head -n 1)

    for required_var in PEQ SRQ PESR_SEP DEPTH_SEP; do
      if [[ -z "${!required_var}" ]]; then
      echo "Failed to extract required cutoff ${required_var} from ~{rf_cutoffs}" >&2
      exit 1
      fi
    done

    function getJavaMem() {
      cat /proc/meminfo | awk -v MEM_FIELD="$1" '{
        f[substr($1, 1, length($1)-1)] = $2
      } END {
        printf "%dM", f[MEM_FIELD] * ~{default="0.85" java_mem_fraction} / 1024
      }'
    }
    JVM_MAX_MEM=$(getJavaMem MemTotal)

    # TrainSVGenotyping writes SR cutoff diagnostics on a best-effort basis. Pre-create the
    # file so a diagnostics failure cannot fail the task on a missing output.
    touch ~{output_name}.sr_cutoff_diagnostics.txt

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
      --rd-depth-min-separation ${DEPTH_SEP} \
      --rd-pesr-min-separation ${PESR_SEP} \
      --output-dir ./ \
      --output-name ~{output_name} \
      ~{training_args}
  >>>

  output {
    File rd_depth_table = "~{output_name}.rd_depth_geno_params.tsv"
    File rd_pesr_table = "~{output_name}.rd_pesr_geno_params.tsv"
    File pe_table = "~{output_name}.pe_geno_params.tsv"
    File sr_table = "~{output_name}.sr_geno_params.tsv"
    File sr_cutoff_diagnostics = "~{output_name}.sr_cutoff_diagnostics.txt"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: gatk_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

# Enforce that SR frequency cutoff training found real signal.
#
# TrainSVGenotyping only reports a degenerate cutoff grid, it does not fail on one, because
# Cromwell delocalizes task outputs solely on success and failing there would leave the
# diagnostics report on the worker. This task runs after that report has been copied out, so it
# can fail loudly while the evidence remains retrievable. It passes sr_table through so the
# genotyping scatter depends on it and cannot start on rejected cutoffs.
task ValidateSRCutoffs {
  input {
    File sr_table
    File sr_cutoff_diagnostics
    Boolean fail_on_degenerate_sr_cutoffs
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    cp ~{sr_table} validated.sr_geno_params.tsv

    # Status keys are written by SplitReadEvidenceGenotyper.appendSelectionStatus
    REJECTED=$(awk -F'\t' '$1 ~ /_selection_status$/ && $2 != "OK" {print $1"="$2}' \
                 ~{sr_cutoff_diagnostics} || true)

    if [ -n "${REJECTED}" ]; then
      echo "SR frequency cutoff selection was rejected: ${REJECTED}" >&2
      awk -F'\t' '$1 ~ /_selection_rejection_reason$/ {print $2}' ~{sr_cutoff_diagnostics} >&2
      echo "Cutoffs fell back to 0.0, which makes the frequency predicate a tautology and" >&2
      echo "disables SR background filtering. Inspect the sr_cutoff_diagnostics output." >&2
      if ~{if fail_on_degenerate_sr_cutoffs then "true" else "false"}; then
        exit 1
      fi
      echo "fail_on_degenerate_sr_cutoffs is false; continuing with unfiltered SR genotypes." >&2
    else
      echo "SR frequency cutoff selection OK for both frequency bins."
    fi
  >>>

  output {
    File validated_sr_table = "validated.sr_geno_params.tsv"
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
    File rd_depth_table
    File rd_pesr_table
    File pe_table
    File sr_table

    String gatk_docker
    String? genotype_args
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

    # Localize RD and PE, but not SR due to interchromosomal BNDs that require access to all contigs' SR evidence
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

    gatk --java-options "-Xmx${JVM_MAX_MEM}" GenotypeSVs \
      -V ~{vcf} \
      -O ~{output_prefix}.vcf.gz \
      ~{"-L " + contig} \
      --median-coverage ~{median_coverage} \
      --rd-file local.rd.txt.gz \
      --discordant-pairs-file local.pe.txt.gz \
      --split-reads-file ~{sr_file} \
      --sequence-dictionary ~{reference_dict} \
      --ploidy-table ~{ploidy_table} \
      --pesr-exclusion-intervals ~{pesr_exclusion_intervals} \
      --depth-exclusion-intervals ~{depth_exclusion_intervals} \
      --rd-depth-table ~{rd_depth_table} \
      --rd-pesr-table ~{rd_pesr_table} \
      --pe-table ~{pe_table} \
      --sr-table ~{sr_table} \
      ~{genotype_args}
  >>>

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
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