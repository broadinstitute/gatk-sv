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
# The integration happens after CallGenomicDisorderCNVs 
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
    String sample_id
    Array[File] gd_output_tarballs
    Array[File] ploidy_tables  # TODO : use joint ploidy table from JoinRawCalls
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
        sample_id = sample_id,
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
    PLOIDY_HEADER=""
    for pt in ~{sep=" " ploidy_tables}; do
      HEADER=$(head -n 1 "$pt")
      if [ $FIRST -eq 1 ]; then
        PLOIDY_HEADER="$HEADER"
        cat "$pt" >> combined_ploidy.tsv
        FIRST=0
      else
        # Differing contig columns would append rows that do not line up with
        # the header the GD tools read, so fail here rather than emitting a
        # table whose trailing contigs silently read back as ploidy 2.
        if [ "$HEADER" != "$PLOIDY_HEADER" ]; then
          echo "ERROR: ploidy table $pt header does not match the first table" >&2
          exit 1
        fi
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
    String sample_id
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

    # --- Drop confident non-carriers, and carry EVIDENCE through ---------------
    # The integrator emits a row for every genomic-disorder locus it screened, and
    # for loci the case does not carry it emits a hom-ref row that also loses the
    # EVIDENCE INFO field present on every input row. A single-sample final VCF has
    # no use for a row saying the case is not a carrier, and Final_VCF_Metrics
    # (svtest) raises on any measured record missing EVIDENCE
    # (src/svtest/svtest/utils/VCFUtils.py:45).
    #
    # Drop hom-ref only. No-call rows (./.) are kept: they record a call whose
    # quality was questionable, not an absence of evidence, and discarding them
    # would erase real results. In this file hom-ref is exactly the 5 integrator
    # rows, while 301 no-call rows are not. GD calls are depth-derived, so the
    # evidence class stamped on any surviving GD row is RD -- "RD" and "RD,PE,SR"
    # are the only values this pipeline writes, and "RD" is within the allowed set
    # of Final_VCF_Metrics.
    #
    # The header is read via `bcftools view -h` into a file rather than through a
    # `gzip -cd | grep -q` pipeline: grep -q exits at the first match, the
    # upstream writer then takes SIGPIPE, and pipefail turns that into exit 141.
    # See 04fa5142; a `|| true` guard here would mask a genuine failure too.
    gd_out="~{prefix}.~{contig}.integrate_gd.vcf.gz"
    gd_filtered="~{prefix}.~{contig}.integrate_gd.filtered.vcf.gz"
    bcftools view -h "${gd_out}" > gd_header.txt
    if ! grep -q '^##INFO=<ID=EVIDENCE,' gd_header.txt; then
      echo "ERROR: ${gd_out} carries no EVIDENCE INFO header; refusing to stamp records" >&2
      exit 1
    fi
    sampleIndex=`bcftools view -h "${gd_out}" | grep '^#CHROM' | cut -f10- | tr "\t" "\n" | awk '$1 == "~{sample_id}" {found=1; print NR - 1} END { if (found != 1) { print "sample not found"; exit 1; }}'`
    bcftools view \
        -e "GT[${sampleIndex}]=\"ref\"" \
        -O v \
        "${gd_out}" \
    | awk \
        '$0 ~ /^#/ { print $0; next; }
        $8 ~ /EVIDENCE=/ { print $0; next; }
        { for(i=1; i<8; ++i) printf "%s\t", $i;
          printf "%s;EVIDENCE=RD", $8;
          for(i=9; i<=NF; ++i) printf "\t%s", $i;
          printf "\n"
        }' \
    | bgzip -c > "${gd_filtered}"
    tabix -p vcf "${gd_filtered}"

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
    File integrated_vcf = "~{prefix}.~{contig}.integrate_gd.filtered.vcf.gz"
    File integrated_vcf_index = "~{prefix}.~{contig}.integrate_gd.filtered.vcf.gz.tbi"
  }
}
