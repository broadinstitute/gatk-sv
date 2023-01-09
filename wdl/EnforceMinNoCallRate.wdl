#### Workflow to compute and filter no-call GT rates from a GATK-SV VCF


version 1.0


import "TasksMakeCohortVcf.wdl" as Tasks


workflow EnforceMinNoCallRate {
  input {
    File vcf
    File vcf_idx
    Int records_per_shard
    Boolean always_shard_vcf = false
    Boolean? reannotate_ncrs_in_vcf
    Array[String]? sample_subset_prefixes
    Array[File]? sample_subset_lists
    File? min_ncr_table  # Note: currently this input does nothing
    Float? global_max_ncr
    String global_ncr_filter_field = "NCR"
    String? chrx_ncr_filter_field
    String? chry_ncr_filter_field

    String sv_base_mini_docker
    String sv_pipeline_docker
  }

  # Step 1: determine if NCR is already defined in the VCF header
  call CheckHeader {
    input:
      vcf=vcf,
      vcf_idx=vcf_idx,
      subset_prefixes=sample_subset_prefixes,
      sv_base_mini_docker=sv_base_mini_docker
  }
  Boolean vcf_is_annotated = if (defined(reannotate_ncrs_in_vcf)) then !reannotate_ncrs_in_vcf else CheckHeader.result

  # Step 2: shard VCF if necessary
  if ( always_shard_vcf || defined(min_ncr_table) || defined(global_max_ncr) || !vcf_is_annotated ) {
    call Tasks.SplitVcf {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        prefix=basename(vcf, ".vcf.gz") + "_sharded",
        min_vars_per_shard=records_per_shard,
        sv_base_mini_docker=sv_base_mini_docker
    }
  }
  Array[File] vcf_shards = select_first([SplitVcf.vcf_shards, [vcf]])

  # # Step 2a: if necessary, annotate all VCF shards with NCR per record
  if ( !vcf_is_annotated ) {
    scatter ( shard in vcf_shards ) {
      call AnnotateNCRs {
        input:
          vcf=shard,
          subset_prefixes=sample_subset_prefixes,
          subset_lists=sample_subset_lists,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }
    call Tasks.ConcatVcfs as ConcatAnnotatedVcfs {
      input:
        vcfs=AnnotateNCRs.annotated_vcf,
        outfile_prefix=basename(vcf, ".vcf.gz") + ".NCR_annotated",
        sv_base_mini_docker=sv_base_mini_docker
    }
  }
  # Step 2b: otherwise, simply extract a table of NCRs for all records
  if ( vcf_is_annotated ) {
    Array[File] vcfs_for_collection = select_first([vcf_shards, [vcf]])
    scatter ( shard in vcfs_for_collection ) {
      call CollectNCRTable {
        input:
          vcf=shard,
          subset_prefixes=sample_subset_prefixes,
          sv_base_mini_docker=sv_base_mini_docker
      }
    }
  }
  call ConcatNCRTables {
    input:
      shards=select_first([AnnotateNCRs.ncr_table, CollectNCRTable.ncr_table]),
      subset_prefixes=sample_subset_prefixes,
      prefix=basename(vcf, ".vcf.gz"),
      sv_base_mini_docker=sv_base_mini_docker
  }

  # Step 3: if min_ncr_table is provided or global_max_ncr is defined, apply min NCR thresholds
  if ( defined(min_ncr_table) || defined(global_max_ncr) ) {
    scatter ( shard in select_first([AnnotateNCRs.annotated_vcf]) ) {
      call ApplyNCRFilter {
        input:
          vcf=shard,
          min_ncr_table=min_ncr_table,
          global_max_ncr=global_max_ncr,
          chrx_ncr_filter_field = chrx_ncr_filter_field,
          chry_ncr_filter_field = chry_ncr_filter_field,
          global_ncr_filter_field=global_ncr_filter_field,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }
    call Tasks.ConcatVcfs as ConcatFilteredVcfs {
      input:
        vcfs=ApplyNCRFilter.filtered_vcf,
        outfile_prefix=basename(vcf, ".vcf.gz") + ".NCR_filtered",
        sv_base_mini_docker=sv_base_mini_docker
    }
  }

  # # Outputs:
  # #   1. VCF with NCR annotated
  # #   2. Table of NCRs for all records prior to filtering
  # #   3. VCF with FILTER labels enforced (only if min_ncr_table was provided)
  output {
    File ncr_annotated_vcf = select_first([ConcatAnnotatedVcfs.concat_vcf, vcf])
    File ncr_annotated_vcf_idx = select_first([ConcatAnnotatedVcfs.concat_vcf_idx, vcf_idx])
    File ncr_table = ConcatNCRTables.merged_table
    File? ncr_filtered_vcf = ConcatFilteredVcfs.concat_vcf
    File? ncr_filtered_vcf_idx = ConcatFilteredVcfs.concat_vcf_idx
  }
}


# Checks VCF header to see if all NCR fields are already defined. Returns Boolean.
task CheckHeader {
  input {
    File vcf
    File vcf_idx
    Array[String]? subset_prefixes
    String sv_base_mini_docker
  }

  Float input_size = size(vcf, "GiB")
  Float compression_factor = 2.0
  Float base_disk_gb = 10.0
  Int total_disk_gb = ceil(base_disk_gb + (compression_factor * input_size))
  runtime {
    memory: "2.0 GiB"
    disks: "local-disk ~{total_disk_gb} HDD"
    cpu: 1
    preemptible: 1
    maxRetries: 1
    docker: sv_base_mini_docker
    bootDiskSizeGb: 10
  }

  command <<<
    set -eu -o pipefail

    cat ~{write_lines(select_first([subset_prefixes, []]))} > subsets.tsv
    cat subsets.tsv
    cut -f1 subsets.tsv | awk '{ print $1"_NCR" }' | cat <( echo "NCR" ) - \
    | awk '{ print "ID="$1 }' > mandatory_infos.txt
    cat mandatory_infos.txt
    n_infos=$( cat mandatory_infos.txt | wc -l )

    if [ $( tabix -H ~{vcf} | fgrep -wf mandatory_infos.txt | wc -l ) -eq $n_infos ]; then
      echo "true" > result.txt
    else
      echo "false" > result.txt
    fi
  >>>

  output {
    Boolean result = read_boolean("result.txt")
  }
}


# Annotates NCRs for all samples and specified subsets for all records in input VCF
task AnnotateNCRs {
  input {
    File vcf
    Array[String]? subset_prefixes
    Array[File]? subset_lists
    String sv_pipeline_docker
  }
  String prefix = basename(vcf, ".vcf.gz")

  Float input_size = size(vcf, "GiB")
  Float compression_factor = 2.0
  Float base_disk_gb = 10.0
  Int total_disk_gb = ceil(base_disk_gb + (2 * compression_factor * input_size))
  runtime {
    memory: "2.0 GiB"
    disks: "local-disk ~{total_disk_gb} HDD"
    cpu: 1
    preemptible: 1
    maxRetries: 1
    docker: sv_pipeline_docker
    bootDiskSizeGb: 10
  }

  command <<<
    set -eu -o pipefail

    # Prepare sample subsets
    if [ ~{defined(subset_prefixes)} == "true" ] && [ ~{defined(subset_lists)} == "true" ]; then
      paste \
        ~{write_lines(select_first([subset_prefixes, []]))} \
        ~{write_lines(select_first([subset_lists, []]))} \
      > subsets.tsv
      annotate_options="--sample-subsets subsets.tsv"
    else
      annotate_options=""
    fi

    # Annotate
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/annotate_nocall_rates.py \
      $annotate_options \
      --stats-tsv ~{prefix}.NCR_stats.tsv \
      ~{vcf} \
      ~{prefix}.NCR_annotated.vcf.gz
    gzip -f ~{prefix}.NCR_stats.tsv
  >>>

  output {
    File annotated_vcf = "~{prefix}.NCR_annotated.vcf.gz"
    File ncr_table = "~{prefix}.NCR_stats.tsv.gz"
  }
}


# Queries NCR (and subset NCRs) for all records in VCF and writes to table
task CollectNCRTable {
  input {
    File vcf
    Array[String]? subset_prefixes
    String sv_base_mini_docker
  }
  String prefix = basename(vcf, ".vcf.gz")

  Float input_size = size(vcf, "GiB")
  Float compression_factor = 2.0
  Float base_disk_gb = 10.0
  Int total_disk_gb = ceil(base_disk_gb + (compression_factor * input_size))
  runtime {
    memory: "2.0 GiB"
    disks: "local-disk ~{total_disk_gb} HDD"
    cpu: 1
    preemptible: 1
    maxRetries: 1
    docker: sv_base_mini_docker
    bootDiskSizeGb: 10
  }

  command <<<
    set -eu -o pipefail

    infos=$( cat ~{write_lines(select_first([subset_prefixes, []]))} \
             | awk -v OFS="_" '{ print $1, "NCR" }' \
             | cat <( echo "NCR" ) - \
             | awk -v OFS="/" -v ORS="\\\t" '{ print "%INFO", $1 }' \
             | sed 's/\\t$/\\n/g' )

    bcftools query -f "%ID\t$infos" ~{vcf} | gzip -c > ~{prefix}.NCRs.tsv.gz
  >>>

  output {
    File ncr_table = "~{prefix}.NCRs.tsv.gz"
  }
}


# Concatenate and reheader NCR tables
task ConcatNCRTables {
  input {
    Array[File] shards
    Array[String]? subset_prefixes
    String prefix
    String sv_base_mini_docker
  }

  Float input_size = size(shards, "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 10.0
  Int total_disk_gb = ceil(base_disk_gb + (2 * compression_factor * input_size))
  runtime {
    memory: "3.5 GiB"
    disks: "local-disk ~{total_disk_gb} HDD"
    cpu: 1
    preemptible: 1
    maxRetries: 1
    docker: sv_base_mini_docker
    bootDiskSizeGb: 10
  }

  command <<<
    set -eu -o pipefail

    # Build header
    cat ~{write_lines(select_first([subset_prefixes, []]))} \
    | awk '{ print $1"_NCR" } ' \
    | cat <( echo -e "#VID\nNCR" ) - \
    | paste -s \
    > ~{prefix}.NCR_stats.tsv

    # Append values & compress
    zcat ~{sep=" " shards} | grep -ve '^#' >> ~{prefix}.NCR_stats.tsv
    gzip ~{prefix}.NCR_stats.tsv
  >>>

  output {
    File merged_table = "~{prefix}.NCR_stats.tsv.gz"
  }
}


# Apply filter status to records in VCF based on NCR
task ApplyNCRFilter {
  input {
    File vcf
    File? vcf_idx
    File? min_ncr_table # Note: currently this input does nothing / feature not yet implemented
    Float? global_max_ncr
    String global_ncr_filter_field = "NCR"
    String? chrx_ncr_filter_field
    String? chry_ncr_filter_field
    String sv_pipeline_docker
  }
  String prefix = basename(vcf, ".vcf.gz")

  Float input_size = size(vcf, "GiB")
  Float compression_factor = 2.0
  Float base_disk_gb = 10.0
  Int total_disk_gb = ceil(base_disk_gb + (2 * compression_factor * input_size))
  runtime {
    memory: "2.0 GiB"
    disks: "local-disk ~{total_disk_gb} HDD"
    cpu: 1
    preemptible: 1
    maxRetries: 1
    docker: sv_pipeline_docker
    bootDiskSizeGb: 10
  }

  command <<<
    set -eu -o pipefail

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -f ~{vcf}
    fi

    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/nocall_rate_filter.py \
      --verbose \
      --global-filter-on ~{global_ncr_filter_field} \
      ~{"--chrx_filter_on " + chrx_ncr_filter_field} \
      ~{"--chry_filter_on " + chry_ncr_filter_field} \
      ~{"--global-max-ncr " + global_max_ncr} \
      ~{vcf} \
      ~{prefix}.NCR_filtered.vcf.gz
  >>>

  output {
    File filtered_vcf = "~{prefix}.NCR_filtered.vcf.gz"
  }
}

