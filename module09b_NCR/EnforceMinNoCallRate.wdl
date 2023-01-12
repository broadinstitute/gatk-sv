#### Copyright (C) 2022 Ryan Collins & the Talkowski Laboratory
####
#### Workflow to compute and filter no-call GT rates from a GATK-SV VCF


version 1.0

import "Structs.wdl"
import "Utils.wdl" as Utils
import "ApplyNoCallRateCutoffs.wdl" as apply_no_call_rate_cutoffs

workflow EnforceMinNoCallRate {
  input {
    File vcf
    File vcf_idx
    Int records_per_shard
    Boolean always_shard_vcf = false
    Boolean? reannotate_ncrs_in_vcf
    Boolean exclude_CTX = true

    Array[String]? sample_subset_prefixes
    Array[File]? sample_subset_lists

    Array[String]? info_col_to_remove

    Array[String]? svtype_list
    Array[String]? ncr_filter_field
    Array[Float]? NCR_cff_list
    String? global_ncr_filter_field
    Float? global_max_ncr

    String sv_base_mini_docker
    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_index_vcf
    RuntimeAttr? runtime_attr_shard_vcf
    RuntimeAttr? runtime_attr_split_vcf_by_type
    RuntimeAttr? runtime_attr_apply_ncr_filter
    RuntimeAttr? runtime_attr_concat_filtered_vcfs
    RuntimeAttr? runtime_attr_exclude_type_from_vcf
    RuntimeAttr? runtime_attr_clean_vcf_info_column


  }
  
  # Step 1: determine if NCR is already defined in the VCF header
  if(defined(info_col_to_remove)){
    call CleanVcfInfoColumn{
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        info_col_to_remove = select_first([info_col_to_remove]),
        sv_pipeline_base_docker = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_clean_vcf_info_column
    }
  }

  File input_vcf = select_first([CleanVcfInfoColumn.cleaned_vcf, vcf])
  File input_vcf_idx = select_first([CleanVcfInfoColumn.cleaned_vcf_idx, vcf_idx])

  call CheckHeader {
    input:
      vcf=input_vcf,
      vcf_idx=input_vcf_idx,
      subset_prefixes=sample_subset_prefixes,
      sv_base_mini_docker=sv_base_mini_docker
  }

  Boolean vcf_is_annotated = select_first([!reannotate_ncrs_in_vcf, CheckHeader.result])

  # Step 2: shard VCF if necessary
  if ( always_shard_vcf ) {
    call Utils.ShardVcf as ShardVcf {
      input:
        vcf=input_vcf,
        vcf_idx=input_vcf_idx,
        prefix=basename(input_vcf, ".vcf.gz") + "_sharded",
        min_vars_per_shard=records_per_shard,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override = runtime_attr_shard_vcf
    }
  }
  Array[File] vcf_shards = select_first([ShardVcf.vcf_shards, [input_vcf]])

  # Step 3: annotate all VCF shards with NCR per record
  if ( !vcf_is_annotated ) {
    scatter ( shard in vcf_shards ) {
      call AnnotateNCRs {
        input:
          vcf=shard,
          subset_prefixes=sample_subset_prefixes,
          subset_lists=sample_subset_lists,
          sv_pipeline_base_docker=sv_pipeline_base_docker
      }
    }

    call Utils.ConcatVcfs as ConcatAnnotatedVcfs {
      input:
        vcfs=AnnotateNCRs.annotated_vcf,
        outfile_prefix=basename(vcf, ".vcf.gz") + ".NCR_annotated",
        sv_base_mini_docker=sv_base_mini_docker
    }
  }

  # Step4: if ncr_cff_list is defined, apply min NCR thresholds

  Array[File] vcf_shards_Step2 = select_first([AnnotateNCRs.annotated_vcf, vcf_shards])

  if ( defined(NCR_cff_list) ){
    scatter (shard_vcf in vcf_shards_Step2){
      call apply_no_call_rate_cutoffs.ApplyNoCallRateCutoffs{
        input:
          vcf = shard_vcf,
          svtype_list  = select_first([svtype_list]),
          ncr_filter_field = select_first([ncr_filter_field]),
          NCR_cff_list = select_first([NCR_cff_list]),
          global_max_ncr = global_max_ncr,
          global_ncr_filter_field = global_ncr_filter_field,
          exclude_CTX = exclude_CTX,

          sv_pipeline_base_docker = sv_pipeline_base_docker,
          sv_base_mini_docker = sv_base_mini_docker,

          runtime_attr_index_vcf = runtime_attr_index_vcf,
          runtime_attr_split_vcf_by_type = runtime_attr_split_vcf_by_type,
          runtime_attr_apply_ncr_filter = runtime_attr_apply_ncr_filter,
          runtime_attr_concat_filtered_vcfs = runtime_attr_concat_filtered_vcfs, 
          runtime_attr_exclude_type_from_vcf = runtime_attr_exclude_type_from_vcf
      }
    }

    call Utils.ConcatVcfs as ConcatFilteredVcfs {
      input:
        vcfs = ApplyNoCallRateCutoffs.ncr_filtered_vcf,
        outfile_prefix = basename(vcf, ".vcf.gz") + ".NCR_filtered",
        sv_base_mini_docker = sv_base_mini_docker
    }
  }


  # Step 5: otherwise, simply extract a table of NCRs for all records
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
    String sv_pipeline_base_docker
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
    docker: sv_pipeline_base_docker
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
    Float? global_max_ncr
    String global_ncr_filter_field
    String sv_pipeline_base_docker
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
    docker: sv_pipeline_base_docker
    bootDiskSizeGb: 10
  }

  command <<<
    set -eu -o pipefail

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -f ~{vcf}
    fi

    script_options=""
    if [ ~{defined(global_max_ncr)} == "true" ]; then
      script_options="$script_options --global-max-ncr ~{global_max_ncr}"
    fi

    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/nocall_rate_filter.py \
      --verbose \
      --global-filter-on ~{global_ncr_filter_field} \
      $script_options \
      ~{vcf} \
      ~{prefix}.NCR_filtered.vcf.gz
  >>>

  output {
    File filtered_vcf = "~{prefix}.NCR_filtered.vcf.gz"
  }
}


task CleanVcfInfoColumn{
    input{
        File vcf
        File vcf_idx
        Array[String] info_col_to_remove
        String sv_pipeline_base_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: ceil(10.0 +  size(vcf, "GB")*2),
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    String prefix = basename(vcf,"vcf.gz")

    command<<<
        set -euo pipefail
        python <<CODE
        
        import os
        import pysam
        fa = open("~{write_lines(info_col_to_remove)}")
        info_col_to_remove = []
        for line in fa:
          pin=line.strip().split()
          info_col_to_remove+=pin
        fa.close()

        fin=pysam.VariantFile("~{vcf}")
        header = fin.header
        for key in header.info.keys():
          if key in info_col_to_remove:
            header.info.remove_header(key)

        fo=pysam.VariantFile("~{prefix}.cleaned.vcf.gz",'w', header = header)
        for record in fin:
          for key in info_col_to_remove:
            if key in record.info.keys():
              del record.info[key]
          fo.write(record)

        fin.close()
        fo.close()
        CODE

        tabix -p vcf "~{prefix}.cleaned.vcf.gz"

    >>>

    output{
        File cleaned_vcf = "~{prefix}.cleaned.vcf.gz"
        File cleaned_vcf_idx = "~{prefix}.cleaned.vcf.gz.tbi"
    }

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_base_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}



