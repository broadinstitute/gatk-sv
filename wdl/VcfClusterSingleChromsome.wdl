version 1.0
# based on snapshot 11
# https://portal.firecloud.org/#methods/Talkowski-SV/04b_vcfcluster_single_chrom/11/wdl

# Copyright (c) 2018 Talkowski Lab

# Contact Ryan Collins <rlcollins@g.harvard.edu>

# Distributed under terms of the MIT License

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks
import "ClusterSingleChromosome.wdl" as VcfClusterTasks

# Workflow to run parallelized vcf clustering for a single chromosome
workflow VcfClusterSingleChrom {
  input {
    Array[File] vcfs
    String prefix
    Int dist
    Float frac
    Float sample_overlap
    File? blacklist
    Array[String] batches
    Int sv_size
    Array[String] sv_types
    String contig
    Int max_shards_per_chrom_svtype
    Int min_variants_per_shard_per_chrom_svtype
    Boolean subset_sr_lists
    File bothside_pass
    File background_fail
    File empty_file

    String sv_pipeline_docker
    String sv_base_mini_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_join_vcfs

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail

    # overrides for VcfClusterTasks
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_concat_sv_types
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_shards
  }
  
  #Remote tabix each vcf & join into a single vcf
  call JoinContigFromRemoteVcfs as JoinVcfs {
    input:
      vcfs=vcfs,
      batches=batches,
      contig=contig,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_join_vcfs
  }

  #Run vcfcluster per chromosome
  call VcfClusterTasks.ClusterSingleChrom as ClusterSingleChrom {
    input:
      vcf=JoinVcfs.joined_vcf,
      contig=contig,
      prefix=prefix,
      max_shards=max_shards_per_chrom_svtype,
      min_per_shard=min_variants_per_shard_per_chrom_svtype,
      dist=dist,
      frac=frac,
      sample_overlap=sample_overlap,
      blacklist=blacklist,
      sv_size=sv_size,
      sv_types=sv_types,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_override_subset_sv_type=runtime_override_subset_sv_type,
      runtime_override_concat_sv_types=runtime_override_concat_sv_types,
      runtime_override_shard_vcf_precluster=runtime_override_shard_vcf_precluster,
      runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
      runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
      runtime_override_concat_shards=runtime_override_concat_shards
  }

  String filtered_bothside_pass_name = prefix + "." + contig + ".pass.VIDs.list"
  String filtered_background_fail_name = prefix + "." + contig + ".fail.VIDs.list"
  if(subset_sr_lists) {
    #Subset bothside_pass & background_fail to chromosome of interest
    call MiniTasks.SubsetVariantList as SubsetBothsidePass {
      input:
        vid_list=bothside_pass,
        vcf=JoinVcfs.joined_vcf,
        outfile_name=prefix + "." + contig + ".pass.VIDs.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_bothside_pass
    }
    call MiniTasks.SubsetVariantList as SubsetBackgroundFail {
      input:
        vid_list=background_fail,
        vcf=JoinVcfs.joined_vcf,
        outfile_name=prefix + "." + contig + ".fail.VIDs.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_background_fail
    }
  }

  output {
    File clustered_vcf = ClusterSingleChrom.clustered_vcf
    File clustered_vcf_idx = ClusterSingleChrom.clustered_vcf_idx
    File filtered_bothside_pass = select_first([SubsetBothsidePass.filtered_vid_list, empty_file])
    File filtered_background_fail = select_first([SubsetBackgroundFail.filtered_vid_list, empty_file])
  }
}


# Task to remote tabix a single chromosome for all VCFs, then merge row-wise
task JoinContigFromRemoteVcfs {
  input {
    Array[File] vcfs
    Array[String] batches
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    vcfs: {
      localization_optional: true
    }
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Int num_vcfs = length(vcfs)
  Float max_vcf_size_gb = 0.5
  Float input_size = max_vcf_size_gb * num_vcfs
  #Float input_size = size([vcf_list, batches_list], "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + 2.0 * compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
  
  command <<<
    set -eu -o pipefail
    
    #Remote tabix all vcfs to chromosome of interest
    1>&2 echo "REMOTE TABIXING VCFs"
    
    # needed for tabix to operate on remote files
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    paste ~{write_lines(batches)} ~{write_lines(vcfs)} | while read BATCH VCF_PATH; do
      1>&2 echo "BATCH=$BATCH"
      1>&2 echo "VCF_PATH=$VCF_PATH"
      if gsutil ls "$VCF_PATH*" | grep -q '\.tbi$' || false; then
        INDEX_PRESENT=1
      else
        INDEX_PRESENT=0
      fi
      if [ $INDEX_PRESENT == 1 ]; then
        1>&2 echo "Found index: $VCF_PATH.tbi"
        TABIX_VCF="$VCF_PATH"
        1>&2 echo "USING $TABIX_VCF"
      else
        1>&2 echo -e "WARNING: no index file for $VCF_PATH\n\tlocalizing and indexing"
        TABIX_VCF=$(basename "$VCF_PATH")
        1>&2 echo "Using $TABIX_VCF"
        gsutil -m cp "$VCF_PATH" .
        tabix -p vcf "$TABIX_VCF"
      fi
      BATCH_VCF="$BATCH.~{contig}.vcf"
      BATCH_VCF=${BATCH_VCF//[[:space:]]/_}
      tabix -h "$TABIX_VCF" "~{contig}:0-300000000"|sed "s/AN=[0-9]*;//g"|sed "s/AC=[0-9]*;//g" > "$BATCH_VCF"
      bgzip -f "$BATCH_VCF"
      if [ $INDEX_PRESENT == 0 ]; then
        rm $TABIX_VCF
      fi

      # echo and pipe batch vcf name to subsetted_vcfs.list, and echo to stderr for debugging purposes
      echo "$BATCH_VCF.gz"
      1>&2 echo "Made $BATCH_VCF.gz"
    done > subsetted_vcfs.list

    1>&2 echo "SANITY CHECK"

    #Sanity check to make sure all subsetted VCFs have same number of records
    # crazy ' || printf ""' statement to avoid pipefail if grep encounters no matching lines
    while read VCF; do
      zcat "$VCF" | (grep -Ev "^#" || printf "") | wc -l
    done < subsetted_vcfs.list \
      > records_per_vcf.txt

    if [ $( sort records_per_vcf.txt | uniq | wc -l ) -gt 1 ]; then
      1>&2 echo "ERROR: INCONSISTENT NUMBER OF RECORDS PER VCF DETECTED"
      cat records_per_vcf.txt
      exit 1
    fi

    1>&2 echo "CALL join_vcfs_paste_implementation.sh"

    #Join vcfs
    /opt/sv-pipeline/04_variant_resolution/scripts/join_vcfs_paste_implementation.sh \
      subsetted_vcfs.list "~{prefix}.joined"

    # more debugging output
    echo "FINISHED join_vcfs_parallel_implementation.sh; RESULTS:"
    find . -name "~{prefix}.joined.vcf*"

    /opt/sv-pipeline/04_variant_resolution/scripts/make_concordant_multiallelic_alts.py \
      $( find . -name "~{prefix}.joined.vcf.gz" ) \
      subsetted_vcfs.list \
      ~{prefix}.unclustered.vcf

    cat ~{prefix}.unclustered.vcf \
      | sed -e 's/:RD,PE,SR/:7/g' \
      | sed -e 's/:PE,SR/:6/g' \
      | sed -e 's/:RD,SR/:5/g' \
      | sed -e 's/:RD,PE/:3/g' \
      | sed -e 's/:PE\t/:2\t/g' -e 's/:SR\t/:4\t/g' -e 's/:RD\t/:1\t/g' \
      | sed -e 's/ID=EV,Number=.,Type=String/ID=EV,Number=1,Type=Integer/g' \
      | bgzip -c > ~{prefix}.unclustered.vcf.gz

    tabix -f -p vcf "~{prefix}.unclustered.vcf.gz"
  >>>

  output {
    File joined_vcf = "~{prefix}.unclustered.vcf.gz"
    File joined_vcf_idx = "~{prefix}.unclustered.vcf.gz.tbi"
  }
}
