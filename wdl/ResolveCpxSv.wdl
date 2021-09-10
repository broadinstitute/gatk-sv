version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "TasksMakeCohortVcf.wdl" as MiniTasks

#Resolve complex SV for a single chromosome
workflow ResolveComplexSv {
  input {
    File vcf
    String prefix
    String contig
    Int max_shard_size
    File cytobands
    File mei_bed
    Array[File] disc_files
    Array[File] rf_cutoff_files
    File pe_exclude_list
    Boolean inv_only
    File ref_dict

    Int precluster_distance
    Float precluster_overlap_frac

    String sv_pipeline_docker
    String sv_base_mini_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_get_se_cutoff
    RuntimeAttr? runtime_override_shard_vcf_cpx
    RuntimeAttr? runtime_override_shard_vids
    RuntimeAttr? runtime_override_resolve_prep
    RuntimeAttr? runtime_override_resolve_cpx_per_shard
    RuntimeAttr? runtime_override_restore_unresolved_cnv_per_shard
    RuntimeAttr? runtime_override_merge_resolve_inner

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_concat_resolved_per_shard
  }

  File vcf_idx = vcf + ".tbi"
  File pe_exclude_list_idx = pe_exclude_list + ".tbi"
  File cytobands_idx = cytobands + ".tbi"
  scatter (i in range(length(disc_files))) {
    File disc_files_idx = disc_files[i]
  }

  # Get SR count cutoff from RF metrics to use in single-ender rescan procedure

  #Shard vcf for complex resolution
  #Note: as of Nov 2, 2018, return lists of variant IDs for each shard. This should
  # dramatically improve sharding speed
  call ShardVcfCpx {
    input:
      vcf=vcf,
      dist=precluster_distance,
      frac=precluster_overlap_frac,
      records_per_shard=max_shard_size,
      prefix="~{prefix}.shard_cpx",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_shard_vcf_cpx
  }

  call MiniTasks.ShardVids {
    input:
      clustered_vcf=ShardVcfCpx.out,
      prefix=prefix,
      records_per_shard=max_shard_size,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_shard_vids
  }

  if (length(ShardVids.out) > 0) {

    call GetSeCutoff {
      input:
        rf_cutoffs=rf_cutoff_files,
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_override_get_se_cutoff
    }

    #Scatter over shards and resolve variants per shard
    scatter ( i in range(length(ShardVids.out)) ) {

      #Prep files for svtk resolve using bucket streaming
      call ResolvePrep {
        input:
          vcf=vcf,
          VIDs_list=ShardVids.out[i],
          chrom=contig,
          disc_files=disc_files,
          disc_files_index=disc_files_idx,
          ref_dict=ref_dict,
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_override_resolve_prep
      }

      #Run svtk resolve
      call SvtkResolve {
        input:
          noref_vcf=ResolvePrep.noref_vcf,
          prefix="~{prefix}.svtk_resolve.shard_~{i}",
          chrom=contig,
          cytobands=cytobands,
          cytobands_idx=cytobands_idx,
          mei_bed=mei_bed,
          pe_exclude_list=pe_exclude_list,
          pe_exclude_list_idx=pe_exclude_list_idx,
          se_pe_cutoff=GetSeCutoff.median_PE_cutoff,
          merged_discfile=ResolvePrep.merged_discfile,
          merged_discfile_idx=ResolvePrep.merged_discfile_idx,
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_override_resolve_cpx_per_shard
      }

      call MergeResolve {
        input:
          full_vcf=ResolvePrep.subsetted_vcf,
          resolved_vcf=SvtkResolve.rs_vcf,
          prefix="~{prefix}.merge_resolve.shard_~{i}",
          noref_vids=ResolvePrep.noref_vids,
          sv_base_mini_docker=sv_base_mini_docker,
          runtime_attr_override=runtime_override_merge_resolve_inner
      }

      #Add unresolved variants back into resolved VCF
      call RestoreUnresolvedCnv as RestoreUnresolvedCnvPerShard {
        input:
          resolved_vcf=MergeResolve.out,
          unresolved_vcf=SvtkResolve.un_vcf,
          prefix="~{prefix}.restore_unresolved.shard_~{i}",
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_override_restore_unresolved_cnv_per_shard
      }
    }

    #Merge across shards
    call MiniTasks.ConcatVcfs as ConcatResolvedPerShard {
      input:
        vcfs=RestoreUnresolvedCnvPerShard.res,
        vcfs_idx=RestoreUnresolvedCnvPerShard.res_idx,
        allow_overlaps=true,
        outfile_prefix=prefix + ".resolved",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_concat_resolved_per_shard
    }
  }

  output {
    File resolved_vcf_merged = select_first([ConcatResolvedPerShard.concat_vcf, vcf])
    File resolved_vcf_merged_idx = select_first([ConcatResolvedPerShard.concat_vcf_idx, vcf_idx])
  }
}


#Get SE cutoff: first quartile of PE cutoff from SR random forest across all batches
#Defaults to 4 if first quartile < 4
task GetSeCutoff {
  input {
    Array[File] rf_cutoffs
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(rf_cutoffs, "GiB")
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + input_size,
    disk_gb: ceil(base_disk_gb + input_size * 2.0),
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

    while read FILE; do
      /opt/sv-pipeline/04_variant_resolution/scripts/convert_poisson_p.py \
        $( awk -F '\t' '{if ( $5=="PE_log_pval") print $2 }' $FILE | head -n1 )
    done < ~{write_lines(rf_cutoffs)} \
      | Rscript -e "cat(max(c(4,floor(quantile(as.numeric(scan('stdin',quiet=T)),probs=0.25)))),sep='\n')" \
      > median_cutoff.txt
  >>>

  output {
    Int median_PE_cutoff = read_lines("median_cutoff.txt")[0]
  }
}


task ShardVcfCpx {
  input {
    File vcf
    String prefix
    Int dist
    Float frac
    Int records_per_shard
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  Float base_disk_gb = 10.0
  Float input_disk_scale = 10.0
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
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
    set -euo pipefail
    bcftools view -G ~{vcf} -Oz -o sites_only.vcf.gz
    svtk vcfcluster <(echo "sites_only.vcf.gz") ~{prefix}.vcf \
      -d ~{dist} \
      -f ~{frac} \
      -p candidate_complex_clusters \
      --svtypes DEL,DUP,INS,INV,BND \
      --ignore-svtypes \
      -o 0 \
      --preserve-header \
      --preserve-ids \
      --skip-merge
    bgzip ~{prefix}.vcf
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}

#Prep files for svtk resolve
task ResolvePrep {
  input {
    File vcf
    File VIDs_list
    File ref_dict
    String chrom
    Array[File] disc_files
    Array[File]? disc_files_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    disc_files: {
      localization_optional: true
    }
    disc_files_index: {
      localization_optional: true
    }
  }

  # sections of disc_files are straemed, but the every operation in this task is record-by-record except
  # bedtools merge, which should only need to keep a few records in memory at a time.
  # assuming memory overhead is fixed
  # assuming disk overhead is input size (accounting for compression) + sum(size of disc_files)
  #  (this is an over-estimate because we only take chunks overlapping VIDs from vcf, but the disk files are not *THAT*
  #   big and disk is cheap)
  Float compressed_input_size = size(vcf, "GiB")
  Float uncompressed_input_size = size([VIDs_list], "GiB")
  Float compression_factor = 30.0
  Float base_disk_gb = 10.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + uncompressed_input_size + compression_factor * compressed_input_size),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  Float mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  runtime {
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    
    # First, subset VCF to variants of interest
    # -uncompress vcf
    zcat "~{vcf}" > uncompressed.vcf
    # -Extract vcf header:
    #   search for first line not starting with '#', stop immediately,
    #   take everything up to that point, then remove last line
    ONLY_HEADER=false
    grep -B9999999999 -m1 -Ev "^#" uncompressed.vcf | sed '$ d' > header.vcf \
      || ONLY_HEADER=true

    if $ONLY_HEADER; then
      # filter is trivial, just copy the vcf
      mv "~{vcf}" input.vcf.gz
    else
      rm -f "~{vcf}"

      N_HEADER=$(wc -l < header.vcf)
      # filter records, concatenate and zip
      tail -n+$((N_HEADER+1)) uncompressed.vcf \
        | { fgrep -wf ~{VIDs_list} || true; } \
        | cat header.vcf - \
        | bgzip -c \
        > input.vcf.gz
      rm -f uncompressed.vcf
    fi

    #Second, extract all-ref variants from VCF. These break svtk resolve with
    # remote tabixing enabled
    svtk vcf2bed input.vcf.gz input.bed
    { grep -Ev "^#" input.bed || true ; } \
      | awk -v FS="\t" '{ if ($6!="") print $4 }' \
      > noref.VIDs.list

    {
      cat header.vcf;
      zcat input.vcf.gz | fgrep -wf noref.VIDs.list || true;
    } \
      | vcf-sort \
      | bgzip -c \
      > noref.vcf.gz
    rm -f header.vcf

    #Third, use GATK to pull down the discfile chunks within ±2kb of all
    # INVERSION breakpoints, and bgzip / tabix
    echo "Forming regions.bed"
    { grep -Ev "^#" input.bed || true; } \
      | (fgrep INV || printf "") \
      | awk -v OFS="\t" -v buffer=2000 \
        '{ print $1, $2-buffer, $2+buffer"\n"$1, $3-buffer, $3+buffer }' \
      | awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3 }' \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | bedtools merge -i - \
      > regions.bed

    if [ -s regions.bed ]; then
      echo "Localizing all discfile shards..."
      DISC_FILE_NUM=0
      while read GS_PATH_TO_DISC_FILE; do
        ((++DISC_FILE_NUM))
        SLICE="disc"$DISC_FILE_NUM"shard"

        if [ -s regions.bed ]; then
          java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
                --skip-header \
                --sequence-dictionary ~{ref_dict} \
                --evidence-file $GS_PATH_TO_DISC_FILE \
                -L regions.bed \
                -O ${SLICE}.PE.txt
        else
          touch ${SLICE}.PE.txt
        fi

        cat ${SLICE}.PE.txt \
          | awk '{ if ($1==$4 && $3==$6) print }' \
          | bgzip -c \
          > ${SLICE}.PE.txt.gz
        rm ${SLICE}.PE.txt
      done < ~{write_lines(disc_files)}
    
      #Fourth, merge PE files and add one artificial pair corresponding to the chromosome of interest
      #This makes it so that svtk doesn't break downstream
      echo "Merging PE files"
      {
        zcat disc*shard.PE.txt.gz;
        echo -e "~{chrom}\t1\t+\t~{chrom}\t2\t+\tDUMMY_SAMPLE_IGNORE";
        echo -e "chr~{chrom}\t1\t+\t~{chrom}\t2\t+\tDUMMY_SAMPLE_IGNORE";
      } \
        | sort -Vk1,1 -k2,2n -k5,5n -k7,7 \
        | bgzip -c \
        > discfile.PE.txt.gz

      rm disc*shard.PE.txt.gz
    else
      echo "No regions to retrieve, making dummy-sample discfile"
      {
        echo -e "~{chrom}\t1\t+\t~{chrom}\t2\t+\tDUMMY_SAMPLE_IGNORE";
        echo -e "chr~{chrom}\t1\t+\t~{chrom}\t2\t+\tDUMMY_SAMPLE_IGNORE";
      } \
        | sort -Vk1,1 -k2,2n -k5,5n -k7,7 \
        | bgzip -c \
        > discfile.PE.txt.gz
    fi

    tabix -s 1 -b 2 -e 2 -f discfile.PE.txt.gz
  >>>

  output {
    File subsetted_vcf = "input.vcf.gz"
    File noref_vcf = "noref.vcf.gz"
    File noref_vids = "noref.VIDs.list"
    File merged_discfile = "discfile.PE.txt.gz"
    File merged_discfile_idx = "discfile.PE.txt.gz.tbi"
  }
}

#Resolve complex SV
task SvtkResolve {
  input {
    File noref_vcf
    String prefix
    String chrom
    File cytobands
    File cytobands_idx
    File mei_bed
    File pe_exclude_list
    File pe_exclude_list_idx
    Int se_pe_cutoff
    File merged_discfile
    File merged_discfile_idx
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  
  String resolved_vcf = "~{prefix}.resolved.unmerged.vcf"
  String unresolved_vcf = "~{prefix}.unresolved.vcf"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(
    [noref_vcf, cytobands, mei_bed, pe_exclude_list, pe_exclude_list_idx,merged_discfile], "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 3 + input_size * 10,
    disk_gb: ceil(10 + input_size * 12),
    cpu_cores: 1,
    preemptible_tries: 1,
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

    #Run svtk resolve on variants after all-ref exclusion
    svtk resolve \
      ~{noref_vcf} \
      ~{resolved_vcf} \
      -p AllBatches_CPX_~{chrom} \
      -u ~{unresolved_vcf} \
      --discfile ~{merged_discfile} \
      --mei-bed ~{mei_bed} \
      --cytobands ~{cytobands} \
      --min-rescan-pe-support ~{se_pe_cutoff} \
      -x ~{pe_exclude_list}

    bgzip ~{resolved_vcf}
    bgzip ~{unresolved_vcf}
  >>>

  output {
    File rs_vcf = resolved_vcf + ".gz"
    File un_vcf = unresolved_vcf + ".gz"
  }  
}

task MergeResolve {
  input {
    File full_vcf
    File resolved_vcf
    String prefix
    File noref_vids
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String out_vcf = "~{prefix}.resolved.vcf.gz"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size([full_vcf, resolved_vcf, noref_vids], "GiB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10 + input_size * 15),
                                  cpu_cores: 1,
                                  preemptible_tries: 1,
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
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    #Add all-ref variants back into resolved VCF
    #Note: requires modifying the INFO field with sed & awk given pysam C bug
    zcat ~{full_vcf} \
      | grep -Ev "^#" \
      | awk 'ARGIND==1{inFileA[$1]; next} {if (!($3 in inFileA)) print }' ~{noref_vids} - OFS='\t' \
      | sed -e 's/;MEMBERS=[^\t]*\t/\t/g' \
      | awk -v OFS="\t" '{ $8=$8";MEMBERS="$3; print }' \
      | cat <(zcat ~{resolved_vcf}) - \
      | vcf-sort \
      | bgzip \
      > ~{out_vcf}
  >>>

  output {
    File out = out_vcf
  }
}

#Restore unresolved CNVs to resolved VCF
task RestoreUnresolvedCnv {
  input {
    File resolved_vcf
    File unresolved_vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String resolved_plus_cnv = "~{prefix}.vcf.gz"

  # straightforward filtering via grep and python script (line-by-line processing)
  #   -> means essentially no memory or disk overhead
  Float input_size = size([resolved_vcf, unresolved_vcf], "GiB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10 + input_size * 20),
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

    # get unresolved records
    zcat ~{unresolved_vcf} \
      | (grep -v "^#" || printf "") \
      > unresolved_records.vcf
    rm "~{unresolved_vcf}"

    # avoid possible obliteration of input file during later processing by writing
    # to temporary file (and postCPX_cleanup.py writing final result to output name)
    zcat ~{resolved_vcf} > ~{resolved_plus_cnv}.tmp
    rm ~{resolved_vcf}

    #Add unresolved CNVs to resolved VCF and wipe unresolved status
    cat unresolved_records.vcf \
      | (fgrep -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" || printf "") \
      | sed -r -e 's/;EVENT=[^;]*;/;/' -e 's/;UNRESOLVED[^;]*;/;/g' \
      | sed -r -e 's/;UNRESOLVED_TYPE[^;]*;/;/g' -e 's/;UNRESOLVED_TYPE[^\t]*\t/\t/g' \
      >> ~{resolved_plus_cnv}.tmp

    #Add other unresolved variants & retain unresolved status (except for inversion single enders)
    cat unresolved_records.vcf \
      | (fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" \
                  -e "INVERSION_SINGLE_ENDER" || printf "") \
      >> ~{resolved_plus_cnv}.tmp

    #Add inversion single enders as SVTYPE=BND
    cat unresolved_records.vcf \
      | (fgrep -v -e "<DEL>" -e "<DUP>" -e "SVTYPE=DEL" -e "SVTYPE=DUP" -e "SVTYPE=CNV" -e "SVTYPE=MCNV" || printf "") \
      | (fgrep -e "INVERSION_SINGLE_ENDER" || printf "") \
      | sed -e 's/SVTYPE=INV/SVTYPE=BND/g' \
      | sed -e 's/END=\([0-9]*\)/END=\1;END2=\1/' \
      >> ~{resolved_plus_cnv}.tmp
    rm unresolved_records.vcf

    #Sort, clean, and compress
    cat ~{resolved_plus_cnv}.tmp \
      | vcf-sort -c \
      | /opt/sv-pipeline/04_variant_resolution/scripts/postCPX_cleanup.py \
        /dev/stdin /dev/stdout \
      | bgzip -c \
      > ~{resolved_plus_cnv}
    tabix ~{resolved_plus_cnv}
  >>>

  output {
    File res = resolved_plus_cnv
    File res_idx = resolved_plus_cnv + ".tbi"
  }
}
