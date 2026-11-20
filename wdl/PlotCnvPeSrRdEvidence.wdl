version 1.0

import "Structs.wdl"

# Plots the PE, SR, and RD evidence of CNV carriers, one plot per variant and
# carrier sample. Adapted from VisualizeDeNovoSVs in
# https://github.com/talkowski-lab/gatk-sv-utils.
workflow PlotCnvPeSrRdEvidence {
  input {
    Array[File] vcfs
    # TSV of variants to plot instead of extracting them from vcfs. The first
    # line must be a header with columns named chr, start, end, svlen, vid,
    # svtype, and sample_id, in any order. There is one row per variant and
    # carrier sample. Only DEL and DUP variants are plotted.
    File? variants
    # variant IDs to plot from vcfs, one per line
    # the variants will still be restricted to DELs and DUPs
    File? variant_ids
    Int min_size
    Int max_carriers_per_variant = 10
    Int? max_svlen
    Int variants_per_shard = 25
    String plot_prefix

    # Batch (sample set) IDs and their sample lists, one sample ID per line.
    # All arrays of batch data must match the order of batch_name_list.
    # Evidence files are streamed from cloud storage, so their indexes must
    # be next to them with a ".tbi" suffix.
    Array[String] batch_name_list
    Array[File] batch_sample_lists
    Array[String] pe_files
    Array[String] sr_files
    Array[String] bincov_files
    Array[String] medians_files

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_extract_cnv_carriers
    RuntimeAttr? runtime_attr_make_manifests
    RuntimeAttr? runtime_attr_shard_variants
    RuntimeAttr? runtime_attr_make_plots
    RuntimeAttr? runtime_attr_merge_plot_tars
  }

  if (!defined(variants)) {
    scatter (vcf in vcfs) {
      call ExtractCnvCarriers {
        input:
          vcf = vcf,
          variant_ids = variant_ids,
          min_size = min_size,
          max_carriers_per_variant = max_carriers_per_variant,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_extract_cnv_carriers
      }
    }
  }

  call MakeManifests {
    input:
      batch_name_list = batch_name_list,
      batch_sample_lists = batch_sample_lists,
      pe_files = pe_files,
      sr_files = sr_files,
      bincov_files = bincov_files,
      medians_files = medians_files,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_make_manifests
  }

  call ShardVariants {
    input:
      variants = if defined(variants) then [select_first([variants])] else select_first([ExtractCnvCarriers.carriers]),
      variants_per_shard = variants_per_shard,
      max_svlen = max_svlen,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_shard_variants
  }

  scatter (i in range(length(ShardVariants.shards))) {
    call MakePlots {
      input:
        variants = ShardVariants.shards[i],
        max_svlen_file = ShardVariants.max_svlens[i],
        evidence_manifest = MakeManifests.evidence_manifest,
        sample_table = MakeManifests.sample_table,
        tar_prefix = "shard_~{i}",
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_make_plots
    }
  }

  call MergePlotTars {
    input:
      plot_tars = MakePlots.plots,
      plot_prefix = plot_prefix,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_plot_tars
  }

  output {
    File cnv_evidence_plots = MergePlotTars.plots
  }
}

task ExtractCnvCarriers {
  input {
    File vcf
    File? variant_ids
    Int min_size
    Int max_carriers_per_variant
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(vcf, "GB") * 3 + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 4.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    vcf='~{vcf}'
    min_size='~{min_size}'
    max_carriers='~{max_carriers_per_variant}'
    variant_ids='~{if defined(variant_ids) then variant_ids else ""}'

    format='%CHROM\t%POS\t%INFO/END\t%INFO/SVLEN\t%ID\t%INFO/SVTYPE\t[%SAMPLE,]\n'
    svtype_filter='(INFO/SVTYPE == "DEL" || INFO/SVTYPE == "DUP")'
    gt_filter='GT = "alt"'
    if [[ -n "${variant_ids:-}" ]]; then
      cp "${variant_ids}" variant_ids.list
      filter="${svtype_filter} & ID=@variant_ids.list & ${gt_filter}"
    else
      filter="${svtype_filter} & FILTER ~ \"PASS\" & INFO/SVLEN >= ${min_size} & ${gt_filter}"
    fi

    printf 'chr\tstart\tend\tsvlen\tvid\tsvtype\tsample_id\n' > carriers.tsv
    bcftools query --include "${filter}" --format "${format}" "${vcf}" \
      | awk -F'\t' -v OFS='\t' -v max_carriers="${max_carriers}" '{
          sub(/,$/, "", $7)
          n = split($7, samples, ",")
          if (max_carriers > 0 && n > max_carriers) {
            n = max_carriers
          }
          svlen = $4 < 0 ? -$4 : $4
          for (i = 1; i <= n; ++i) {
            print $1, $2, $3, svlen, $5, $6, samples[i]
          }
        }' >> carriers.tsv
  >>>

  output {
    File carriers = "carriers.tsv"
  }
}

task MakeManifests {
  input {
    Array[String] batch_name_list
    Array[File] batch_sample_lists
    Array[String] pe_files
    Array[String] sr_files
    Array[String] bincov_files
    Array[String] medians_files
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(batch_sample_lists, "GB") * 2 + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    batch_names='~{write_lines(batch_name_list)}'
    batch_sample_lists='~{write_lines(batch_sample_lists)}'
    pe_files='~{write_lines(pe_files)}'
    sr_files='~{write_lines(sr_files)}'
    bincov_files='~{write_lines(bincov_files)}'
    medians_files='~{write_lines(medians_files)}'

    n_batches=$(wc -l < "${batch_names}")
    if [[ "${n_batches}" -eq 0 ]]; then
      echo "ERROR: at least one batch is required" >&2
      exit 1
    fi
    for f in "${batch_sample_lists}" "${pe_files}" "${sr_files}" "${bincov_files}" "${medians_files}"; do
      if [[ $(wc -l < "${f}") -ne "${n_batches}" ]]; then
        echo "ERROR: all arrays of batch data must have the same length as batch_name_list (${n_batches})" >&2
        exit 1
      fi
    done

    # headerless TSV of sample ID and batch ID
    : > sample_table.tsv
    while IFS=$'\t' read -r batch samples_file; do
      awk -F'\t' -v batch="${batch}" 'NF > 0 { print $1 "\t" batch }' "${samples_file}" >> sample_table.tsv
    done < <(paste "${batch_names}" "${batch_sample_lists}")

    # headerless TSV of batch ID, PE, SR, RD, and median coverage paths
    paste "${batch_names}" "${pe_files}" "${sr_files}" "${bincov_files}" "${medians_files}" \
      > evidence_manifest.tsv
  >>>

  output {
    File sample_table = "sample_table.tsv"
    File evidence_manifest = "evidence_manifest.tsv"
  }
}

task ShardVariants {
  input {
    Array[File] variants
    Int variants_per_shard
    Int? max_svlen
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(variants, "GB") * 3 + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    variants='~{write_lines(variants)}'
    variants_per_shard='~{variants_per_shard}'
    max_svlen='~{select_first([max_svlen, -1])}'

    mkdir shards max_svlens
    # Rows of the same variant are kept in the same shard so the evidence of
    # each variant is only queried once. Each shard gets a file with the
    # largest SV length in the shard, used to size the plotting task.
    xargs -a "${variants}" gawk -F'\t' -v n="${variants_per_shard}" -v max_svlen="${max_svlen}" '
      BEGIN {
        OFS = "\t"
        nreq = split("chr start end svlen vid svtype sample_id", req, " ")
        header = "chr\tstart\tend\tsvlen\tvid\tsvtype\tsample_id"
      }
      FNR == 1 {
        delete col
        for (i = 1; i <= NF; ++i) {
          col[$i] = i
        }
        for (i = 1; i <= nreq; ++i) {
          if (!(req[i] in col)) {
            printf("ERROR: %s does not have a \"%s\" column\n", FILENAME, req[i]) > "/dev/stderr"
            exit 1
          }
        }
        next
      }
      {
        vid = $col["vid"]
        svlen = $col["svlen"] + 0
      }
      max_svlen > 0 && svlen > max_svlen { next }
      vid != last_vid && nvids % n == 0 {
        if (outpath) {
          close(outpath)
          print shard_max_svlen > svlen_path
          close(svlen_path)
          shard_max_svlen = 0
        }
        outpath = sprintf("shards/%06d.tsv", k)
        svlen_path = sprintf("max_svlens/%06d", k)
        ++k
        print header > outpath
      }
      {
        if (vid != last_vid) {
          ++nvids
          last_vid = vid
        }
        line = $col[req[1]]
        for (i = 2; i <= nreq; ++i) {
          line = line OFS $col[req[i]]
        }
        print line > outpath
        shard_max_svlen = shard_max_svlen >= svlen ? shard_max_svlen : svlen
      }
      END {
        if (outpath) {
          print shard_max_svlen > svlen_path
        }
      }'

    if [[ -z "$(ls shards)" ]]; then
      echo "ERROR: no variants to plot" >&2
      exit 1
    fi
  >>>

  output {
    # The shards and max SV length files share the same names so the relative
    # ordering of the globs is the same.
    Array[File] shards = glob("shards/*.tsv")
    Array[File] max_svlens = glob("max_svlens/*")
  }
}

task MakePlots {
  input {
    File variants
    File max_svlen_file
    File evidence_manifest
    File sample_table
    String tar_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Int max_svlen = read_int(max_svlen_file)
  Int default_cpus = if max_svlen >= 10000000 then if max_svlen >= 50000000 then 8 else 4 else 2
  Float default_mem_gb = if max_svlen >= 10000000 then if max_svlen >= 50000000 then 64.0 else 32.0 else 16.0

  RuntimeAttr default_attr = object {
    cpu_cores: default_cpus,
    mem_gb: default_mem_gb,
    disk_gb: 64,
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    # gorilla streams the evidence with the system tabix, which must be built
    # with GCS support and recent fixes for reading remote files
    export PATH="/opt/htslib/bin:${PATH}"

    variants='~{variants}'
    evidence_manifest='~{evidence_manifest}'
    sample_table='~{sample_table}'
    tar_prefix='~{tar_prefix}'

    Rscript /opt/sv-pipeline/scripts/visualize_cnv.R \
      "${variants}" "${evidence_manifest}" "${sample_table}" \
      "${tar_prefix}" exclusions.tsv

    mv exclusions.tsv "${tar_prefix}/exclusions.tsv"
    tar -cf "${tar_prefix}.tar" "${tar_prefix}"
  >>>

  output {
    File plots = "~{tar_prefix}.tar"
  }
}

task MergePlotTars {
  input {
    Array[File] plot_tars
    String plot_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float disk_size = size(plot_tars, "GB") * 3 + 16

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2.0,
    disk_gb: ceil(disk_size),
    boot_disk_gb: 8,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

  command <<<
    set -o errexit
    set -o nounset
    set -o pipefail

    plot_tars='~{write_lines(plot_tars)}'
    plot_prefix='~{plot_prefix}'

    mkdir store
    xargs -a "${plot_tars}" -I '{}' tar -xf '{}' -C store

    # sort the plots by SV type; plot names end in _<svtype>_<sample>.png
    mkdir -p "${plot_prefix}/DEL" "${plot_prefix}/DUP"
    find store -type f -name '*_DEL_*.png' -exec mv -t "${plot_prefix}/DEL" '{}' +
    find store -type f -name '*_DUP_*.png' -exec mv -t "${plot_prefix}/DUP" '{}' +

    find store -type f -name 'exclusions.tsv' -exec cat '{}' + \
      | awk 'NR == 1 || $0 !~ /^sample_id\t/' > "${plot_prefix}/exclusions.tsv"

    tar -czf "${plot_prefix}.tar.gz" "${plot_prefix}"
  >>>

  output {
    File plots = "~{plot_prefix}.tar.gz"
  }
}
