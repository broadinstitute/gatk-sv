version 1.0

import "MinGQRocOpt.wdl" as roc_opt_sub
import "CalcAF.wdl" as calcAF
import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "ReviseSVtypeINStoMEI.wdl" as ReviseSVtype

workflow Module07MinGQ {
  input {
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_updates_docker
    File vcf
    File vcf_idx
    String prefix
    File contiglist
    File trios_famfile
    String optimize_minSizes
    String optimize_maxSizes
    String optimize_minFreqs
    String optimize_maxFreqs
    String optimize_includeSVTYPEs
    String optimize_includeFILTERs
    String optimize_excludeFILTERs
    String optimize_includeEV
    String optimize_excludeEV
    Int optimize_maxSVperTrio
    Float roc_max_fdr_PCRMINUS
    Float roc_max_fdr_PCRPLUS
    Int roc_min_gq
    Int roc_max_gq
    Int roc_step_gq
    Int roc_shards
    Int min_sv_per_proband_per_condition
    Int max_shards_per_chrom_step1
    Int min_records_per_shard_step1
    Float max_noCallRate
    Int global_minGQ
    String ref_build
    File? sanders_2015_tarball
    File? collins_2017_tarball
    File? werling_2018_tarball
    File? pcrplus_samples_list
    File? pcrminus_filter_lookup_table
    Boolean MingqTraining=! defined(pcrminus_filter_lookup_table)

    # overrides for local tasks
    RuntimeAttr? runtime_attr_CombineVcfs
    RuntimeAttr? runtime_attr_GatherTrioData
    RuntimeAttr? runtime_attr_ReviseSVtypeMEI
    RuntimeAttr? runtime_override_split_vcf_to_clean

    RuntimeAttr? runtime_attr_compute_shard_af
    RuntimeAttr? runtime_attr_scatter_vcf_calcaf
    RuntimeAttr? runtime_attr_combine_sharded_vcfs_calcaf
  }

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Get svtype of MEI
  call ReviseSVtype.ReviseSVtypeINStoMEI as ReviseSVtypeMEI {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      sv_base_mini_docker = sv_base_mini_docker,
      sv_pipeline_docker = sv_pipeline_docker,
      prefix = prefix,
      contiglist = contiglist,
      max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
      min_records_per_shard_step1 = min_records_per_shard_step1,
      runtime_attr_ReviseSVtypeMEI = runtime_attr_ReviseSVtypeMEI,
      runtime_override_split_vcf_to_clean=runtime_override_split_vcf_to_clean
  }

  # Get list of PCRMINUS samples
  call GetSampleLists {
    input:
      vcf = ReviseSVtypeMEI.updated_vcf,
      vcf_idx = ReviseSVtypeMEI.updated_vcf_idx,
      pcrplus_samples_list = pcrplus_samples_list,
      prefix = prefix,
      sv_base_mini_docker = sv_base_mini_docker
  }

  # Shard VCF per-chromosome and add AF annotation
  scatter ( contig in contigs ) {
    #Split VCF into PCR+ and PCR-
    call calcAF.CalcAF as getAFs {
      input:
        vcf=ReviseSVtypeMEI.updated_vcf,
        vcf_idx=ReviseSVtypeMEI.updated_vcf_idx,
        sv_per_shard=1000,
        prefix="~{prefix}.~{contig[0]}",
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_updates_docker=sv_pipeline_updates_docker,
        runtime_attr_scatter_vcf = runtime_attr_scatter_vcf_calcaf,
        runtime_attr_compute_shard_af = runtime_attr_compute_shard_af,
        runtime_attr_combine_sharded_vcfs = runtime_attr_combine_sharded_vcfs_calcaf
    }
    if (defined(pcrplus_samples_list)) {
      call SplitPcrVcf {
        input:
          vcf=getAFs.vcf_wAFs,
          prefix="~{prefix}.~{contig[0]}",
          pcrplus_samples_list=select_first([pcrplus_samples_list]),
          sv_base_mini_docker=sv_base_mini_docker
      }
    }
    File pcr_minus_vcf = select_first([SplitPcrVcf.PCRMINUS_vcf, getAFs.vcf_wAFs])

    # Dev note Feb 18 2021: the output from cat_AF_table_PCRMINUS is a required
    # input to Module07XfBatchEffect.wdl, so the subsequent three tasks always 
    # need to be generated (even if passing a precomputed minGQ cutoff table)

    # Annotate PCR-specific AFs
    call calcAF.CalcAF as getAFs_byPCR {
      input:
        vcf=ReviseSVtypeMEI.updated_vcf,
        vcf_idx=ReviseSVtypeMEI.updated_vcf_idx,
        sv_per_shard=1000,
        prefix="~{prefix}.~{contig[0]}",
        sample_pop_assignments=GetSampleLists.sample_PCR_labels,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_updates_docker=sv_pipeline_updates_docker,
        runtime_attr_scatter_vcf = runtime_attr_scatter_vcf_calcaf,
        runtime_attr_compute_shard_af = runtime_attr_compute_shard_af,
        runtime_attr_combine_sharded_vcfs = runtime_attr_combine_sharded_vcfs_calcaf
    }
    # Gather table of AC/AN/AF for PCRPLUS and PCRMINUS samples
    call GetAfTables {
      input:
        vcf=getAFs_byPCR.vcf_wAFs,
        pcrplus_samples_list=pcrplus_samples_list,
        vcf_idx=getAFs_byPCR.vcf_wAFs_idx,
        prefix="~{prefix}.~{contig[0]}",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }
  call CombineRocOptResults as cat_AF_table_PCRMINUS {
    input:
      shards=GetAfTables.PCRMINUS_AF_table,
      outfile="~{prefix}.PCRMINUS.AF_preMinGQ.txt",
      sv_base_mini_docker=sv_base_mini_docker,
  }


  if (MingqTraining) {
    ###PCRMINUS
    call SplitFamfile as SplitFamfile_PCRMINUS {
      input:
        vcf=pcr_minus_vcf[0],
        vcf_idx=pcr_minus_vcf[0] + ".tbi",
        famfile=trios_famfile,
        fams_per_shard=1,
        prefix="~{prefix}.PCRMINUS",
        sv_base_mini_docker=sv_base_mini_docker
    }
    scatter ( fam in SplitFamfile_PCRMINUS.famfile_shards ) {
      call CollectTrioSVdat as CollectTrioSVdat_PCRMINUS {
        input:
          vcf_shards=pcr_minus_vcf,
          famfile=fam,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }
    call GatherTrioData as GatherTrioData_PCRMINUS {
      input:
        files=CollectTrioSVdat_PCRMINUS.trio_SVdata,
        prefix="~{prefix}.PCRMINUS",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_GatherTrioData
    }


    # Get table of all conditions to evaluate
    call EnumerateConditions {
      input:
        prefix=prefix,
        condition_shards=roc_shards,
        optimize_minSizes=optimize_minSizes,
        optimize_maxSizes=optimize_maxSizes,
        optimize_minFreqs=optimize_minFreqs,
        optimize_maxFreqs=optimize_maxFreqs,
        optimize_includeSVTYPEs=optimize_includeSVTYPEs,
        optimize_includeFILTERs=optimize_includeFILTERs,
        optimize_excludeFILTERs=optimize_excludeFILTERs,
        optimize_includeEV=optimize_includeEV,
        optimize_excludeEV=optimize_excludeEV,
        sv_pipeline_docker=sv_pipeline_docker
    }


    # Scatter over each shard of conditions and send the trio data for ROC optimization
    scatter ( shard in EnumerateConditions.minGQ_conditions_table_noHeader_shards ) {
      ### PCRMINUS
      call roc_opt_sub.MinGQRocOpt as roc_opt_PCRMINUS {
        input:
          trio_tarball=GatherTrioData_PCRMINUS.tarball,
          prefix="~{prefix}.PCRMINUS",
          trios_list=SplitFamfile_PCRMINUS.cleaned_trios_famfile,
          conditions_table=shard,
          maxSVperTrio=optimize_maxSVperTrio,
          roc_max_fdr=roc_max_fdr_PCRMINUS,
          roc_min_gq=roc_min_gq,
          roc_max_gq=roc_max_gq,
           roc_step_gq=roc_step_gq,
          min_sv_per_proband_per_condition=min_sv_per_proband_per_condition,
          sv_base_mini_docker=sv_base_mini_docker,
          sv_pipeline_docker=sv_pipeline_docker
      }
    }


    # Merge ROC results to build minGQ filtering lookup tree
    ###PCRMINUS
    call CombineRocOptResults as combine_roc_optimal_PCRMINUS {
      input:
        shards=roc_opt_PCRMINUS.roc_optimal_merged,
        outfile="~{prefix}.PCRMINUS.minGQ_condition_opts.txt",
        sv_base_mini_docker=sv_base_mini_docker
    }
    call CombineRocOptResults as combine_roc_stats_PCRMINUS {
      input:
        shards=roc_opt_PCRMINUS.distrib_stats_merged,
        outfile="~{prefix}.minGQ_condition_distrib_stats.txt",
        sv_base_mini_docker=sv_base_mini_docker
    }


    # Create final minGQ filtering tree
    ###PCRMINUS
    call BuildFilterTree as build_tree_PCRMINUS {
      input:
        conditions_table=EnumerateConditions.minGQ_conditions_table,
        condition_optimizations=combine_roc_optimal_PCRMINUS.merged_file,
        condition_distrib_stats=combine_roc_stats_PCRMINUS.merged_file,
        prefix="~{prefix}.PCRMINUS",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  # Apply filter per chromosome
  ###PCRMINUS
  scatter ( vcf_shard in pcr_minus_vcf ) {
    call ApplyMinGQFilter as apply_filter_PCRMINUS {
      input:
        vcf=vcf_shard,
        minGQ_lookup_table=select_first([pcrminus_filter_lookup_table,build_tree_PCRMINUS.filter_lookup_table]),
        prefix="~{prefix}.PCRMINUS",
        PCR_status="PCRMINUS",
        maxNCR=max_noCallRate,
        global_minGQ=global_minGQ,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }


  call MiniTasks.ConcatVcfs as CombineVcfs {
    input:
      vcfs=apply_filter_PCRMINUS.filtered_vcf,
      naive=true,
      outfile_prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_CombineVcfs
  }


  # Final output
  output {
    File filtered_vcf = CombineVcfs.concat_vcf
    File filtered_vcf_idx = CombineVcfs.concat_vcf_idx
    File? AF_table_preMinGQ_PCRMINUS = cat_AF_table_PCRMINUS.merged_file
    File? filter_lookup_table = build_tree_PCRMINUS.filter_lookup_table
  }
}

# Get lists of PCRPLUS and PCRMINUS samples present in input VCF
task GetSampleLists {
  input{
    String sv_base_mini_docker
    File vcf
    File vcf_idx
    File? pcrplus_samples_list
    String prefix
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10 + size(vcf, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    bcftools query -l ~{vcf} > all_samples.list
    if ~{defined(pcrplus_samples_list)}; then
      awk -v OFS="\t" 'ARGIND==1{inFileA[$1]; next} {if($1 in inFileA){print $1,"PCRPLUS"}else{print $1,"PCRMINUS"}}' ~{pcrplus_samples_list} all_samples.list \
        > ~{prefix}.PCR_status_assignments.txt
    else
      awk -v OFS="\t" '{ print $1, "PCRMINUS" }' all_samples.list \
        > ~{prefix}.PCR_status_assignments.txt
    fi
  >>>

  output {
    File sample_PCR_labels = "~{prefix}.PCR_status_assignments.txt"
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


# Split a VCF into two parts, corresponding to PCR+ and PCR-
task SplitPcrVcf {
  input{
    File vcf
    String prefix
    File pcrplus_samples_list
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    bcftools query -l ~{vcf} > all_samples.list
    awk 'ARGIND==1{inFileA[$1]; next} !($1 in inFileA)' ~{pcrplus_samples_list} all_samples.list \
      > pcrminus_samples.list
    bcftools reheader -s pcrminus_samples.list -Oz -o ~{prefix}.PCRMINUS.vcf.gz
    tabix ~{prefix}.PCRMINUS.vcf.gz
  >>>

  output {
    File PCRMINUS_vcf = "~{prefix}.PCRMINUS.vcf.gz"
    File PCRMINUS_vcf_idx = "~{prefix}.PCRMINUS.vcf.gz.tbi"
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


# Get a simple table with ID/AC/AN/AF per variant, prior to minGQ
task GetAfTables {
  input{
    File vcf
    File vcf_idx
    String prefix
    String sv_pipeline_docker
    File? pcrplus_samples_list
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: ceil(10 + size(vcf, "GB") * 3),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    #Run vcf2bed
    svtk vcf2bed --info ALL --no-samples ~{vcf} "~{prefix}.vcf2bed.bed"
    #Cut to necessary columns
    idxs=$( sed -n '1p' "~{prefix}.vcf2bed.bed" \
            | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ print $1, NR }' \
            | grep -e 'name\|SVLEN\|SVTYPE\|_AC\|_AN\|_CN_NONREF_COUNT\|_CN_NUMBER' \
            | fgrep -v "OTH" \
            | cut -f2 \
            | paste -s -d\, || true )
    cut -f"$idxs" "~{prefix}.vcf2bed.bed" \
      | sed 's/^name/\#VID/g' \
      | gzip -c \
      > "~{prefix}.frequencies.preclean.txt.gz"
    if [ ! -z "~{pcrplus_samples_list}" ]; then
      echo -e "dummy\tPCRMINUS\ndummy2\tPCRPLUS" > dummy.tsv
    else
      echo -e "dummy\tPCRMINUS" > dummy.tsv
    fi
    #Clean frequencies
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/clean_frequencies_table.R \
      "~{prefix}.frequencies.preclean.txt.gz" \
      dummy.tsv \
      "~{prefix}.frequencies.txt"
    for PCR in $( cut -f2 dummy.tsv | sort | uniq ); do
      AC_idx=$( zcat "~{prefix}.frequencies.txt.gz" | sed -n '1p' | sed 's/\t/\n/g' | awk -v PCR="$PCR" '{ if ($1==PCR"_AC") print NR }' )
      AN_idx=$( zcat "~{prefix}.frequencies.txt.gz" | sed -n '1p' | sed 's/\t/\n/g' | awk -v PCR="$PCR" '{ if ($1==PCR"_AN") print NR }' )
      zcat "~{prefix}.frequencies.txt.gz" \
        | sed '1d' \
        | awk -v FS="\t" -v OFS="\t" -v AC="$AC_idx" -v AN="$AN_idx" \
          '{ print $1, $(AC), $(AN) }' \
      > ~{prefix}."$PCR".AF_preMinGQ.txt
    done
    if [ ! -z ~{prefix}.PCRPLUS.AF_preMinGQ.txt ]; then
      touch ~{prefix}.PCRPLUS.AF_preMinGQ.txt
    fi
  >>>

  output {
    File PCRPLUS_AF_table = "~{prefix}.PCRPLUS.AF_preMinGQ.txt"
    File PCRMINUS_AF_table = "~{prefix}.PCRMINUS.AF_preMinGQ.txt"
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

# Shard a trio famfile to keep only trios that are all represented in the vcf header
task SplitFamfile {
  input{
    File vcf
    File vcf_idx
    File famfile
    String prefix
    String sv_base_mini_docker
    Int fams_per_shard
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    #Get list of sample IDs & column numbers from VCF header
    tabix -H ~{vcf} | fgrep -v "##" | sed 's/\t/\n/g' \
      | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt
    #Iterate over families & subset VCF
    while read famID pro fa mo prosex pheno; do
      pro_idx=$( awk -v ID=$pro '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      fa_idx=$( awk -v ID=$fa '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      mo_idx=$( awk -v ID=$mo '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
      if ! [ -z $pro_idx ] && ! [ -z $fa_idx ] && ! [ -z $mo_idx ]; then
        fgrep -w "$famID" ~{famfile} || true
      fi
    done < ~{famfile} \
      | awk -v FS="\t" -v OFS="\t" '{ if ($2!="0" && $3!="0" && $4!="0") print $0 }' \
    > "~{prefix}.cleaned_trios.fam"  # proband-only
    split -l ~{fams_per_shard} --numeric-suffixes=00001 -a 5 ~{prefix}.cleaned_trios.fam famfile_shard_
  >>>

  output {
    File cleaned_trios_famfile = "~{prefix}.cleaned_trios.fam"
    Array[File] famfile_shards = glob("famfile_shard_*")
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

# Collect a single table of all relevant variants for a single family
task CollectTrioSVdat {
  input{
    Array[File] vcf_shards
    File famfile
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


  command <<<
    for wrapper in 1; do
      #Write header
      echo -e "#famID\tVID\tSVLEN\tAF\tSVTYPE\tFILTER\tpro_EV\tpro_AC\tfa_AC\tmo_AC\tpro_GQ\tfa_GQ\tmo_GQ"
      #Iterate over list of VCF shards
      while read vcf; do
        #Get list of sample IDs & column numbers from VCF header
        zfgrep "#" $vcf | fgrep -v "##" | head -n1000 |sed 's/\t/\n/g' \
          | awk -v OFS="\t" '{ print $1, NR }' > vcf_header_columns.txt
        #Iterate over families & subset VCF
        while read famID pro fa mo prosex pheno; do
          pro_idx=$( awk -v ID=$pro '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
          fa_idx=$( awk -v ID=$fa '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
          mo_idx=$( awk -v ID=$mo '{ if ($1==ID) print $2 }' vcf_header_columns.txt )
          if ! [ -z $pro_idx ] && ! [ -z $fa_idx ] && ! [ -z $mo_idx ]; then
            #Subset vcf to only multiallelic sites in teh family
            zcat "$vcf" | cut -f1-9,"$pro_idx","$fa_idx","$mo_idx" \
            | grep -e '\#\|[0-1]\/1\|MULTIALLELIC' \
            | bgzip -c > $famID.vcf.gz
            #Get list of CNVs in proband that are ≥5kb have ≥50% coverage in either parent
            svtk vcf2bed -i SVTYPE --no-header $famID.vcf.gz stdout \
            | awk -v OFS="\t" '{ if ($NF ~ /DEL|DUP|CNV/) print $1, $2, $3, $4, $NF, $6 }' \
            > $famID.CNVs.bed
            fgrep -w $pro $famID.CNVs.bed \
            | awk -v OFS="\t" '{ if ($3-$2>=5000 && $5!="CNV") print $1, $2, $3, $4, $5 }' \
            > $pro.CNVs.gt5kb.bed
            fgrep -w $fa $famID.CNVs.bed > $fa.CNVs.bed
            fgrep -w $mo $famID.CNVs.bed > $mo.CNVs.bed
            #Deletions
            awk -v OFS="\t" '{ if ($NF=="DEL") print $0, "1" }' $pro.CNVs.gt5kb.bed \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DEL|CNV/) print $0 }' $fa.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $1, $2, $3, $4, $5, $6, $NF }' \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DEL|CNV/) print $0 }' $mo.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $4, $6, $7, $NF }' \
            > $famID.RD_genotype_update.txt
            #Duplications
            awk -v OFS="\t" '{ if ($NF=="DUP") print $0, "1" }' $pro.CNVs.gt5kb.bed \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DUP|CNV/) print $0 }' $fa.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $1, $2, $3, $4, $5, $6, $NF }' \
            | bedtools coverage -a - \
              -b <( awk '{ if ($5 ~ /DUP|CNV/) print $0 }' $mo.CNVs.bed ) \
            | awk -v OFS="\t" '{ if ($NF>=0.5) $NF=1; else $NF=0; print $4, $6, $7, $NF }' \
            >> $famID.RD_genotype_update.txt
            #Get variant stats
            /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/gather_trio_genos.py \
              --ac-adj $famID.RD_genotype_update.txt \
              --no-header \
              $famID.vcf.gz stdout "$pro" "$fa" "$mo" \
            | awk -v famID="$famID" -v OFS="\t" '{ print famID, $0 }'
          fi
        done < ~{famfile}
      done < ~{write_lines(vcf_shards)}
    done | bgzip -c > "trio_variant_info.txt.gz"
  >>>

  output {
    File trio_SVdata = "trio_variant_info.txt.gz"
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


# Gather all trio SV data into a single tarball (helps with Cromwell file localization)
task GatherTrioData {
  input{
    Array[File] files
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    tar -czvf ~{prefix}.tar.gz -T ~{write_lines(files)}
  >>>

  output {
    File tarball = "~{prefix}.tar.gz"
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


# Enumerate all minGQ conditions to test
task EnumerateConditions {
  input{
    String prefix
    Int condition_shards
    String optimize_minSizes
    String optimize_maxSizes
    String optimize_minFreqs
    String optimize_maxFreqs
    String optimize_includeSVTYPEs
    String optimize_includeFILTERs
    String optimize_excludeFILTERs
    String optimize_includeEV
    String optimize_excludeEV
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/create_minGQ_tranches_table.R \
      --min.sizes "~{optimize_minSizes}" \
      --max.sizes "~{optimize_maxSizes}" \
      --min.freqs "~{optimize_minFreqs}" \
      --max.freqs "~{optimize_maxFreqs}" \
      --svtype.include "~{optimize_includeSVTYPEs}" \
      --filter.include "~{optimize_includeFILTERs}" \
      --filter.exclude "~{optimize_excludeFILTERs}" \
      --ev.include "~{optimize_includeEV}" \
      --ev.exclude "~{optimize_excludeEV}" \
      "~{prefix}.minGQ_conditions.txt"
    fgrep -v "#" "~{prefix}.minGQ_conditions.txt" \
      > "~{prefix}.minGQ_conditions.noHeader.txt"
    /opt/sv-pipeline/04_variant_resolution/scripts/evenSplitter.R \
      -S ~{condition_shards} \
      "~{prefix}.minGQ_conditions.noHeader.txt" \
      "~{prefix}.minGQ_conditions.noHeader.shard"
  >>>

  output {
    File minGQ_conditions_table = "~{prefix}.minGQ_conditions.txt"
    File minGQ_conditions_table_noHeader = "~{prefix}.minGQ_conditions.noHeader.txt"
    Array[File] minGQ_conditions_table_noHeader_shards = glob("~{prefix}.minGQ_conditions.noHeader.shard*")
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


# Merge ROC optimal cutoffs or stats
task CombineRocOptResults {
  input{
    Array[File] shards
    String outfile
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    cat ~{write_lines(shards)} | xargs -I {} fgrep -v "#" {}  | sort -Vk1,1 > ~{outfile}
  >>>

  output {
    File merged_file = "~{outfile}"
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


# Build final minGQ filtering tree
task BuildFilterTree {
  input{
    File conditions_table
    File condition_optimizations
    File condition_distrib_stats
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/create_minGQ_lookup_table.R \
      "~{conditions_table}" \
      "~{condition_distrib_stats}" \
      "~{condition_optimizations}" \
      "~{prefix}.minGQ.ordered_tree_hierarchy.txt" \
      "~{prefix}.minGQ.filter_lookup_table.txt"
  >>>


  output {
    File ordered_tree_hierarchy = "~{prefix}.minGQ.ordered_tree_hierarchy.txt"
    File filter_lookup_table = "~{prefix}.minGQ.filter_lookup_table.txt"
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


# Apply minGQ filter to VCF
task ApplyMinGQFilter {
  input{
    File vcf
    File minGQ_lookup_table
    String prefix
    String PCR_status
    Float maxNCR
    Int global_minGQ
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override    
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/apply_minGQ_filter.py \
      --minGQ "~{global_minGQ}" \
      --maxNCR "~{maxNCR}" \
      --simplify-INS-SVTYPEs \
      --cleanAFinfo \
      --prefix "~{PCR_status}" \
      "~{vcf}" \
      "~{minGQ_lookup_table}" \
      stdout \
    | fgrep -v "##INFO=<ID=AN," \
    | fgrep -v "##INFO=<ID=AC," \
    | fgrep -v "##INFO=<ID=AF," \
    | fgrep -v "##INFO=<ID=N_BI_GENOS," \
    | fgrep -v "##INFO=<ID=N_HOMREF," \
    | fgrep -v "##INFO=<ID=N_HET," \
    | fgrep -v "##INFO=<ID=N_HOMALT," \
    | fgrep -v "##INFO=<ID=FREQ_HOMREF," \
    | fgrep -v "##INFO=<ID=FREQ_HET," \
    | fgrep -v "##INFO=<ID=FREQ_HOMALT," \
    | fgrep -v "##INFO=<ID=CN_NUMBER," \
    | fgrep -v "##INFO=<ID=CN_COUNT," \
    | fgrep -v "##INFO=<ID=CN_FREQ," \
    | fgrep -v "##INFO=<ID=CN_NONREF_COUNT," \
    | fgrep -v "##INFO=<ID=CN_NONREF_FREQ," \
    | bgzip -c \
    > "~{prefix}.minGQ_filtered.vcf.gz"
  >>>

  output {
    File filtered_vcf = "~{prefix}.minGQ_filtered.vcf.gz"
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


# Merge PCRPLUS and PCRMINUS VCFs for a single chromosome
task MergePcrVCFs {
  input{
    File? PCRPLUS_vcf
    File PCRMINUS_vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 4,
    disk_gb: 20,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    if [ ! -z "~{PCRPLUS_vcf}" ];then
      #Sanitize FILTER columns
      zcat "~{PCRPLUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRPLUS_filters.txt
      zcat "~{PCRMINUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRMINUS_filters.txt
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/merge_filter_columns.py \
        PCRPLUS_filters.txt \
        PCRMINUS_filters.txt \
        merged_filters.txt
      #Write new VCF header
      zgrep -e '^##' ~{PCRPLUS_vcf} > "~{prefix}.minGQ_filtered.vcf"
      zgrep -e '^##' ~{PCRMINUS_vcf} | fgrep "NOCALL_RATE" >> "~{prefix}.minGQ_filtered.vcf"
      #Column-wise merger
      paste \
        <( zgrep -ve '^##' "~{PCRPLUS_vcf}" | cut -f1-6 ) \
        <( cat <( echo -e "FILTER" ) merged_filters.txt ) \
        <( zgrep -ve '^##' "~{PCRPLUS_vcf}" | cut -f8- ) \
        <( zgrep -ve '^##' "~{PCRMINUS_vcf}"  | cut -f10- ) \
      >> "~{prefix}.minGQ_filtered.vcf"
      /opt/sv-pipeline/scripts/drop_empty_records.py \
        "~{prefix}.minGQ_filtered.vcf" \
        "~{prefix}.minGQ_filtered.no_blanks.vcf"
      #Bgzip & tabix
      bgzip -f "~{prefix}.minGQ_filtered.no_blanks.vcf"
    else
      #Sanitize FILTER columns
      zcat "~{PCRMINUS_vcf}" | cut -f7 | grep -ve '^#' | sed '1d' > PCRMINUS_filters.txt
      #Write new VCF header
      zcat "~{PCRMINUS_vcf}" | sed -n '1,1000p' | grep -e '^##' > "~{prefix}.minGQ_filtered.vcf"
      #Column-wise merger
      paste \
        <( zcat "~{PCRMINUS_vcf}" | grep -ve '^##' | cut -f1-6 ) \
        <( cat <( echo -e "FILTER" ) PCRMINUS_filters.txt ) \
        <( zcat "~{PCRMINUS_vcf}" | grep -ve '^##' | cut -f8- ) \
      >> "~{prefix}.minGQ_filtered.vcf"
      /opt/sv-pipeline/scripts/drop_empty_records.py \
        "~{prefix}.minGQ_filtered.vcf" \
        "~{prefix}.minGQ_filtered.no_blanks.vcf"
      #Bgzip & tabix
      bgzip -f "~{prefix}.minGQ_filtered.no_blanks.vcf"
    fi
  >>>

  output {
    File merged_vcf = "~{prefix}.minGQ_filtered.no_blanks.vcf.gz"
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
