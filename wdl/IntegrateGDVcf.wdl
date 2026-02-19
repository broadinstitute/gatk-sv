version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort

# Integrate Genomic Disorder CNV calls into the pipeline's final VCF, scattered
# across chromosomes for parallelism.
#
# This workflow extracts gd_cnv_calls.tsv.gz from each batch tarball and
# concatenates them into a combined GD calls table (PrepareGDCallsTask), then
# scatters the `gatk-sv-gd integrate` CLI across chromosomes. Each shard subsets
# the VCF (tabix), GD calls table (col 4 = chrom), and GD regions table (col 1 =
# CHROM) before running integration. Per-chromosome VCFs are merged with ConcatVcfs.
#
# The integration happens after MakeCohortVcf/FilterGenotypes/SampleQC filtering
# and before functional consequence and allele frequency annotations.
#
# Inputs:
#   vcf                 - The filtered VCF from the pipeline (cohort-level)
#   gd_output_tarballs  - One tarball per batch from CallGenomicDisorderCNVs
#   ploidy_tables       - One ploidy table per batch (wide format, one row per sample)
#   gd_table            - GD regions table (same file used for GD calling)
#   par_bed             - PAR regions BED
#   contig_list         - List of contigs to scatter over (one per line)

workflow IntegrateGDVcf {
  input {
    File vcf
    File vcf_index
    String prefix
    Array[File] gd_output_tarballs
    Array[File] ploidy_tables
    File gd_table
    File par_bed
    File contig_list
    String sv_pipeline_docker
    String sv_base_mini_docker
    String? integrate_args

    RuntimeAttr? runtime_attr_override_prepare
    RuntimeAttr? runtime_attr_override_integrate
    RuntimeAttr? runtime_attr_override_concat
  }

  call PrepareGDCallsTask {
    input:
      gd_output_tarballs = gd_output_tarballs,
      ploidy_tables = ploidy_tables,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_override_prepare
  }

  scatter (contig in read_lines(contig_list)) {
    call IntegrateGDVcfTask {
      input:
        vcf = vcf,
        vcf_index = vcf_index,
        prefix = prefix,
        contig = contig,
        combined_gd_calls = PrepareGDCallsTask.combined_gd_calls,
        combined_ploidy = PrepareGDCallsTask.combined_ploidy,
        gd_table = gd_table,
        par_bed = par_bed,
        sv_pipeline_docker = sv_pipeline_docker,
        integrate_args = integrate_args,
        runtime_attr_override = runtime_attr_override_integrate
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs = IntegrateGDVcfTask.integrated_vcf,
      vcfs_idx = IntegrateGDVcfTask.integrated_vcf_index,
      naive = true,
      outfile_prefix = prefix + ".integrate_gd",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_override_concat
  }

  output {
    File integrate_gd_vcf = ConcatVcfs.concat_vcf
    File integrate_gd_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}

task PrepareGDCallsTask {
  input {
    Array[File] gd_output_tarballs
    Array[File] ploidy_tables
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(gd_output_tarballs, "GiB")

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 8.0,
    disk_gb: ceil(200.0 + input_size),
    boot_disk_gb: 20,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # --- Extract & concatenate GD calls from all batch tarballs ---
    mkdir -p gd_calls
    for tarball in ~{sep=" " gd_output_tarballs}; do
      tar -xzf "$tarball" -C gd_calls/
    done

    # Concatenate all gd_cnv_calls.tsv.gz (preserve one header)
    FIRST=1
    cat /dev/null > combined_gd_calls.tsv
    for calls in $(find gd_calls -name "gd_cnv_calls.tsv.gz" | sort); do
      if [ $FIRST -eq 1 ]; then
        gunzip -c "$calls" >> combined_gd_calls.tsv
        FIRST=0
      else
        gunzip -c "$calls" | awk 'NR>1' >> combined_gd_calls.tsv
      fi
    done
    gzip combined_gd_calls.tsv

    # --- Combine per-batch ploidy tables (wide format, one row per sample) ---
    FIRST=1
    cat /dev/null > combined_ploidy.tsv
    for pt in ~{sep=" " ploidy_tables}; do
      if [ $FIRST -eq 1 ]; then
        cat "$pt" >> combined_ploidy.tsv
        FIRST=0
      else
        awk 'NR>1' "$pt" >> combined_ploidy.tsv
      fi
    done
  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }

  output {
    File combined_gd_calls = "combined_gd_calls.tsv.gz"
    File combined_ploidy = "combined_ploidy.tsv"
  }
}

task IntegrateGDVcfTask {
  input {
    File vcf
    File vcf_index
    String prefix
    String contig
    File combined_gd_calls
    File combined_ploidy
    File? gd_table
    File? par_bed
    String sv_pipeline_docker
    String? integrate_args

    RuntimeAttr? runtime_attr_override
  }

  Float vcf_size = size(vcf, "GiB")

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: ceil(50.0 + vcf_size * 2),
    boot_disk_gb: 20,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # Subset VCF to contig
    tabix -h ~{vcf} ~{contig} | bgzip -c > subset.vcf.gz
    tabix -p vcf subset.vcf.gz

    # Subset GD calls table to contig (col 4 = chrom)
    { zcat ~{combined_gd_calls} | awk 'NR==1'; \
      zcat ~{combined_gd_calls} | awk -v c="~{contig}" 'NR>1 && $4==c'; } \
      | gzip > combined_gd_calls.~{contig}.tsv.gz

    # Subset GD table to contig if defined (col 1 = CHROM, has header)
    GD_TABLE="~{select_first([gd_table, ""])}"
    if [ -n "${GD_TABLE}" ]; then
      { (zcat "${GD_TABLE}" 2>/dev/null || cat "${GD_TABLE}") | awk 'NR==1'; \
        (zcat "${GD_TABLE}" 2>/dev/null || cat "${GD_TABLE}") | awk -v c="~{contig}" 'NR>1 && $1==c'; } \
        | gzip > gd_table.~{contig}.tsv.gz
    fi

    # --- Run integration ---
    gatk-sv-gd integrate \
      --vcf subset.vcf.gz \
      --gd-calls combined_gd_calls.~{contig}.tsv.gz \
      ~{if defined(gd_table) then "--gd-table gd_table." + contig + ".tsv.gz" else ""} \
      ~{if defined(par_bed) then "--par-bed " + par_bed else ""} \
      --ploidy-table ~{combined_ploidy} \
      --out-vcf ~{prefix}.~{contig}.integrate_gd.vcf.gz \
      --temp-dir $(pwd) \
      ~{default="" integrate_args}

  >>>

  runtime {
    docker: sv_pipeline_docker
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }

  output {
    File integrated_vcf = "~{prefix}.~{contig}.integrate_gd.vcf.gz"
    File integrated_vcf_index = "~{prefix}.~{contig}.integrate_gd.vcf.gz.tbi"
  }
}
