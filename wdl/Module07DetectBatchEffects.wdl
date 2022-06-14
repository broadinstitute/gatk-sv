##########################
## EXPERIMENTAL WORKFLOW
##########################

## Base script:   https://portal.firecloud.org/#workspaces/gnomad-v3-sv/1KGP_2504/method-configs/Talkowski-SV/check_batch_effects_CMC_copy

version 1.0

import "prune_add_af.wdl" as calcAF
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow DetectBatchEffects {
  input{
    File vcf
    File vcf_idx
    File sample_batch_assignments
    File batches_list
    File pcrminus_batches_list
    File sample_pop_assignments
    File excludesamples_list #empty file if need be
    File famfile
    File contiglist
    File? par_bed
    Int variants_per_shard
    String prefix
    File af_pcrminus_premingq

    # Optional inputs if PCR+ samples are in callset    
    File? pcrplus_batches_list
    File? af_pcrplus_premingq
    File? pcrplus_samples_list

    String sv_pipeline_docker
    String sv_pipeline_updates_docker
    String sv_pipeline_base_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_merge_labeled_vcfs
  }
    Array[String] batches = read_lines(batches_list)
    Array[Array[String]] contigs = read_tsv(contiglist)
    
  # Shard VCF per batch, compute pops-specific AFs, and convert to table of VID & AF stats
  scatter ( batch in batches ) {
    # Get list of samples to include & exclude per batch
    call GetBatchSamplesList {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        batch=batch,
        sample_batch_assignments=sample_batch_assignments,
        probands_list=excludesamples_list,
        sv_pipeline_docker=sv_pipeline_docker
    }
    # Prune VCF to samples 
    call calcAF.prune_and_add_vafs as getAFs {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        prefix=batch,
        sample_pop_assignments=sample_pop_assignments,
        prune_list=GetBatchSamplesList.exclude_samples_list,
        famfile=famfile,
        sv_per_shard=25000,
        contiglist=contiglist,
        drop_empty_records="FALSE",
        par_bed=par_bed,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_updates_docker=sv_pipeline_updates_docker
    }
    # Get minimal table of AF data per batch, split by ancestry
    call GetFreqTable {
      input:
        vcf=getAFs.output_vcf,
        sample_pop_assignments=sample_pop_assignments,
        prefix=batch,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  # Merge frequency results per batch into a single table of all variants with AF data across batches
  call MergeFreqTables {
    input:
      tables=GetFreqTable.freq_data,
      batches_list=batches_list,
      rows_per_shard=variants_per_shard,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }
  call MergeFreqTables as MergeFreqTables_allPops {
    input:
      tables=GetFreqTable.freq_data_allPops,
      batches_list=batches_list,
      rows_per_shard=variants_per_shard,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Test 1: compare frequencies before and after minGQ, and generate list 
  # of variants that are significantly different between the steps
  call CompareFreqsPrePostMinGQ {
    input:
      af_pcrminus_premingq=af_pcrminus_premingq,
      af_pcrplus_premingq=af_pcrplus_premingq,
      AF_postMinGQ_table=MergeFreqTables_allPops.merged_table,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Compute AF stats per pair of batches & determine variants with batch effects
  scatter ( af_table_shard in MergeFreqTables.merged_table_shards ) {

    # Test 2: cross-batch comparison of frequencies for PCR- batches
    call CrossBatchPropTest as MinusPropTest {
      input:
        freq_table=af_table_shard,
        batches_list=pcrminus_batches_list,
        prefix="~{prefix}.PCRMINUS_cross_batch_test",
        sv_pipeline_base_docker=sv_pipeline_base_docker
    }

    # Tests 3 & 4 only run if pcrplus_batches_list is provided
    if ( defined(pcrplus_batches_list) ) {

        # Test 3: cross-batch comparison of frequencies for PCR- batches
        call CrossBatchPropTest as PlusPropTest {
          input:
            freq_table=af_table_shard,
            batches_list=select_first([pcrplus_batches_list]),
            prefix="~{prefix}.PCRPLUS_cross_batch_test",
            sv_pipeline_base_docker=sv_pipeline_base_docker
        }

        # Test 4: Chi-square test of aggregated AFs between PCR+ and PCR- batches
        call ComparePlusMinusFreqs {
          input:
            freq_table=af_table_shard,
            minus_batches_list=pcrminus_batches_list,
            plus_batches_list=select_first([pcrminus_batches_list]),
            prefix="~{prefix}.aggregated_PCRPLUS_vs_PCRMINUS_test",
            sv_pipeline_base_docker=sv_pipeline_base_docker
        }
    }
  }
  # TODO: Combine results from tests 1-4 into single tsv with P-values per test
  # Need to add task here


  # Test 5: Perform one-vs-all comparison of AFs per batch to find batch-specific sites
  call PrepOneVsAllPairsList as PrepMinusPairs {
    input:
      batches_list=pcrminus_batches_list,
      prefix="~{prefix}.PCRMINUS_one_vs_all",
      sv_base_mini_docker=sv_base_mini_docker
  }
  if ( defined(pcrplus_batches_list) ) {
    call PrepOneVsAllPairsList as PrepPlusPairs {
      input:
        batches_list=select_first([pcrplus_batches_list]),
        prefix="~{prefix}.PCRPLUS_one_vs_all",
        sv_base_mini_docker=sv_base_mini_docker
    }
  }
  scatter ( af_table_shard in MergeFreqTables.merged_table_shards ) {
    call CompareBatchPairs as OneVsAllMinus {
      input:
        freq_table=af_table_shard,
        batch_pairs_list=PrepMinusPairs.pairs_list,
        prefix="~{prefix}.PCRMINUS_one_vs_all",
        sv_pipeline_base_docker=sv_pipeline_base_docker
    }
    if ( defined(pcrplus_batches_list) ) {
      call CompareBatchPairs as OneVsAllPlus {
        input:
          freq_table=af_table_shard,
          batch_pairs_list=select_first([PrepPlusPairs.pairs_list]),
          prefix="~{prefix}.PCRPLUS_one_vs_all",
          sv_pipeline_base_docker=sv_pipeline_base_docker
      }
    }
  }
  # TODO: Collect results from one-vs-all pairwise batch effect detection
  # Need to add task here

  # TODO: need to update this task based on the revised outputs of the above tasks
  # # Distill final table of variants to be reclassified
  # call MakeReclassificationTable {
  #   input:
  #     freq_table=MergeFreqTables.merged_table,
  #     pairwise_fails=merge_pairwise_checks.fails_per_variant_all,
  #     onevsall_fails=merge_one_vs_all_checks.fails_per_variant_all,
  #     prefix=prefix,
  #     pairwise_cutoff = pairwise_cutoff,
  #     onevsall_cutoff = onevsall_cutoff,
  #     sv_pipeline_docker=sv_pipeline_docker
  # }

  # TODO: need to update this task based on the revised outputs of the above tasks
  # Apply batch effect labels
  # scatter ( contig in contigs ) {
  #   call ApplyBatchEffectLabels as apply_labels_perContig {
  #     input:
  #       vcf=vcf,
  #       vcf_idx=vcf_idx,
  #       contig=contig[0],
  #       reclassification_table=MakeReclassificationTable.reclassification_table,
  #       mingq_prePost_pcrminus_fails=CompareFreqsPrePostMinGQ.pcrminus_fails,
  #       mingq_prePost_pcrplus_fails=CompareFreqsPrePostMinGQ.pcrplus_fails,
  #       pcrplus_sample_ids=pcrplus_samples_list,
  #       prefix="~{prefix}.~{contig[0]}",
  #       sv_pipeline_docker=sv_pipeline_docker
  #   }
  # }
  # call MiniTasks.ConcatVcfs as merge_labeled_vcfs {
  #   input:
  #     vcfs=apply_labels_perContig.labeled_vcf,
  #     vcfs_idx=apply_labels_perContig.labeled_vcf_idx,
  #     outfile_prefix="~{prefix}.batch_effects_labeled_merged",
  #     sv_base_mini_docker=sv_pipeline_docker,
  #     runtime_attr_override=runtime_attr_merge_labeled_vcfs
  # }

  output {
    # File labeled_vcf = merge_labeled_vcfs.concat_vcf
    # File labeled_vcf_idx = merge_labeled_vcfs.concat_vcf_idx
  }
}


# Get list of samples to include & exclude per batch
# Always exclude probands from all batches
task GetBatchSamplesList {
  input{
    File vcf
    File vcf_idx
    String batch
    File sample_batch_assignments
    File probands_list
    String sv_pipeline_docker
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
    set -euo pipefail
    # Get list of all samples present in VCF header
    tabix -H ~{vcf} | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' | sort -Vk1,1 \
    > all_samples.list
    # Get list of samples in batch
    fgrep -w ~{batch} ~{sample_batch_assignments} | cut -f1 \
    | fgrep -wf - all_samples.list \
    | fgrep -wvf ~{probands_list} \
    > "~{batch}.samples.list" || true
    # Get list of samples not in batch
    fgrep -wv ~{batch} ~{sample_batch_assignments} | cut -f1 \
    | cat - ~{probands_list} | sort -Vk1,1 | uniq \
    | fgrep -wf - all_samples.list \
    > "~{batch}.exclude_samples.list" || true
  >>>

  output {
    File include_samples_list = "~{batch}.samples.list"
    File exclude_samples_list = "~{batch}.exclude_samples.list"
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


# Run vcf2bed and subset to just include VID, SVTYPE, SVLEN, _AC, and _AN
task GetFreqTable {
  input{
    File vcf
    File sample_pop_assignments
    String prefix
    String sv_pipeline_docker
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
    set -euo pipefail
    #Run vcf2bed
    svtk vcf2bed \
      --info ALL \
      --no-samples \
      ~{vcf} "~{prefix}.vcf2bed.bed"
    ### Create table of freqs by ancestry
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
    #Clean frequencies (note that this script automatically gzips the output file taken as the last argument)
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/clean_frequencies_table.R \
      "~{prefix}.frequencies.preclean.txt.gz" \
      "~{sample_pop_assignments}" \
      "~{prefix}.frequencies.txt"
    ### Create table of freqs, irrespective of ancestry
    #Cut to necessary columns
    idxs=$( sed -n '1p' "~{prefix}.vcf2bed.bed" \
            | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ if ($1=="name" || $1=="SVLEN" || $1=="SVTYPE" || $1=="AC" || $1=="AN" || $1=="CN_NUMBER" || $1=="CN_NONREF_COUNT") print NR }' \
            | paste -s -d\, || true )
    cut -f"$idxs" "~{prefix}.vcf2bed.bed" > minfreq.subset.bed
    svtype_idx=$( sed -n '1p' minfreq.subset.bed \
                  | sed 's/\t/\n/g' \
                  | awk -v OFS="\t" '{ if ($1=="SVTYPE") print NR }' || true )
    awk -v OFS="\t" -v sidx="$svtype_idx" '{ if ($sidx=="CNV" || $sidx=="MCNV") print $1, $2, $3, $6, $7; else print $1, $2, $3, $4, $5 }' minfreq.subset.bed \
    | sed 's/^name/\#VID/g' \
    | gzip -c \
    > "~{prefix}.frequencies.allPops.txt.gz"
  >>>

  output {
    File freq_data = "~{prefix}.frequencies.txt.gz"
    File freq_data_allPops = "~{prefix}.frequencies.allPops.txt.gz"
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


# Combine frequency data across batches
task MergeFreqTables {
  input{
    Array[File] tables
    File batches_list
    Int rows_per_shard
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 100 + (5 * size(tables, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail

    #Get list of batch IDs and batch table paths
    while read batch; do
      echo "$batch"
      find ./ -name "$batch.frequencies*txt.gz" | sed -n '1p'
    done < ~{batches_list} | paste - - \
    > input.list

    #Make sure all input files have the same number of lines
    while read batch file; do
      zcat "$file" | wc -l
    done < input.list > nlines.list
    nlines=$( sort nlines.list | uniq | wc -l )
    if [ "$nlines" -gt 1 ]; then
      echo "AT LEAST ONE INPUT FILE HAS A DIFFERENT NUMBER OF LINES"
      exit 0
    fi

    #Prep files for paste joining
    echo "PREPPING FILES FOR MERGING"
    while read batch file; do
        #Header
        zcat "$file" | sed -n '1p' | cut -f1-3
        #Body
        zcat "$file" | sed '1d' \
        | sort -Vk1,1 \
        | cut -f1-3
    done < <( sed -n '1p' input.list ) \
    > header.txt
    while read batch file; do
      for wrapper in 1; do
        #Header
        zcat "$file" | sed -n '1p' \
        | cut -f4- | sed 's/\t/\n/g' \
        | awk -v batch="$batch" '{ print $1"."batch }' \
        | paste -s
        #Body
        zcat "$file" | sed '1d' \
        | sort -Vk1,1 \
        | cut -f4-
      done > "$batch.prepped.txt" 
    done < input.list

    #Join files with simple paste
    paste \
      header.txt \
      $( awk -v ORS=" " '{ print $1".prepped.txt" }' input.list ) \
    | gzip -c \
    > "~{prefix}.merged_AF_table.txt.gz"

    # Shard column-wise joined file by # of records
    zcat "~{prefix}.merged_AF_table.txt.gz" | head -n1 > shard_header.tsv || true
    zcat "~{prefix}.merged_AF_table.txt.gz" | sed '1d' > shard_body.tsv
    /opt/sv-pipeline/04_variant_resolution/scripts/evenSplitter.R \
      -L ${rows_per_shard} \
      shard_body.tsv \
      "~{prefix}.sharded_"
    n_shards=$( find ./ -name "~{prefix}.sharded_*" | wc -l )
    for i in $( seq 1 $n_shards ); do
      cat shard_header.tsv "~{prefix}.sharded_$i" | gzip -c \
      > "~{prefix}.merged_AF_table.shard_$i.txt.gz"
    done

    # Get total number of records (required for downstream steps)
    cut -f1 shard_body.tsv | wc -l > n_records.txt
  >>>

  output {
    File merged_table = "~{prefix}.merged_AF_table.txt.gz"
    Array[File] merged_table_shards = glob("~{prefix}.merged_AF_table.shard_*.txt.gz")
    Int n_records = read_int("n_records.txt")
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


# Compare frequencies before and after minGQ
task CompareFreqsPrePostMinGQ {
  input{
  File af_pcrminus_premingq
  File? af_pcrplus_premingq
  File AF_postMinGQ_table
  String prefix
  String sv_pipeline_docker
  RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 8,
    disk_gb: 30,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  
  command <<<
    set -euo pipefail

    # Different behavior whether or not af_pcrplus_premingq is provided
    if [ "~{defined(af_pcrplus_premingq)}" == "true" ]; then
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/compare_freqs_pre_post_minGQ.R \
        ~{af_pcrminus_premingq} \
        ~{af_pcrplus_premingq} \
        ~{AF_postMinGQ_table} \
        ./ \
        "~{prefix}."
    else
      /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/compare_freqs_pre_post_minGQ.PCRMinus_only.R \
        ~{af_pcrminus_premingq} \
        ~{AF_postMinGQ_table} \
        ./ \
        "~{prefix}."
    fi
  >>>

  output {
    File pcrminus_fails = "~{prefix}.PCRMINUS_minGQ_AF_prePost_fails.VIDs.list"
    File? pcrplus_fails = "~{prefix}.PCRPLUS_minGQ_AF_prePost_fails.VIDs.list"
    File minGQ_prePost_comparison_data = "~{prefix}.minGQ_AF_prePost_comparison.data.txt.gz"
    File minGQ_prePost_comparison_plot = "~{prefix}.minGQ_AF_prePost_comparison.plot.png"
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


# Test of equality of AFs across batches
task CrossBatchPropTest {
  input {
    File freq_table
    File batches_list
    String prefix
    String sv_pipeline_base_docker
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

    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/cross_batch_AF_equality_test.R \
      ~{freq_table} \
      ~{batches_list} \
      "~{prefix}.results.tsv"
    gzip "~{prefix}.results.tsv"
  >>>

  output {
    File results = "~{prefix}.results.tsv.gz"
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


# Test of equality of AFs summed for PCR+ batches and for PCR- batches
task ComparePlusMinusFreqs {
  input {
    File freq_table
    File minus_batches_list
    File plus_batches_list
    String prefix
    String sv_pipeline_base_docker
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

    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/compare_AFs_between_PCR_batches.R \
      ~{freq_table} \
      ~{minus_batches_list} \
      ~{plus_batches_list} \
      "~{prefix}.results.tsv"
    gzip "~{prefix}.results.tsv"
  >>>

  output {
    File results = "~{prefix}.results.tsv.gz"
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


# Prepare batch pairs list for one-vs-all comparisons
task PrepOneVsAllPairsList {
  input {
    File batches_list
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 2.5,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    awk -v OFS="\t" '{ print $1, "ALL_OTHERS" }' ~{batches_list} > ~{prefix}.batch_pairs.tsv
  >>>

  output {
    File pairs_list = "~{prefix}.batch_pairs.tsv"
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


# Compare AF stats per variant between a list of batch pairs
task CompareBatchPairs {
  input{
    File freq_table
    File batch_pairs_list
    Int total_n_records
    String prefix
    String sv_pipeline_base_docker
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

    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/pairwise_batch_effect_test.R \
      ~{freq_table} \
      ~{batch_pairs_list} \
      "~{prefix}.results.txt"
    gzip "~{prefix}.results.txt"
  >>>

  output {
    File results = "~{prefix}.results.txt.gz"
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



# Merge lists of batch effect checks and count total number of times each variant failed
task MergeVariantFailureLists {
  input{
    Array[File] fail_variant_lists
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
    #Write list of paths to all batch effect variant lists
    #Get master list of PCR+ to PCR+ failures   #removed from the PCR- only projects
    #Get master list of PCR- to PCR- failures   #removed from the PCR- only projects
    #Get master list of PCR+ to PCR- failures   #removed from the PCR- only projects
    #Get master list of all possible failures
    cat ~{write_lines(fail_variant_lists)} \
    | xargs -I {} cat {} \
    | sort -Vk1,1 | uniq -c \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > "~{prefix}.all.failures.txt" || true
  >>>

  output {
    File fails_per_variant_all = "~{prefix}.all.failures.txt"
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


# Consolidate all batch effect check results into a single table with reclassification per variant
task MakeReclassificationTable {
  input{
    File freq_table
    File pairwise_fails
    File onevsall_fails
    String prefix
    Int? pairwise_cutoff
    Int? onevsall_cutoff
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 8,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/make_batch_effect_reclassification_table.PCRMinus_only.R \
      ~{freq_table} \
      ~{pairwise_fails} \
      ~{onevsall_fails} \
      "~{prefix}.batch_effect_reclassification_table.txt" \
      ~{pairwise_cutoff} \
      ~{onevsall_cutoff}
  >>>

  output {
    File reclassification_table = "~{prefix}.batch_effect_reclassification_table.txt"
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


# Apply batch effect labels to VCF
task ApplyBatchEffectLabels {
  input{
    File vcf
    File vcf_idx
    String contig
    File reclassification_table
    File mingq_prePost_pcrminus_fails
    File? mingq_prePost_pcrplus_fails
    File? pcrplus_sample_ids
    String prefix
    String sv_pipeline_docker
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
    set -euo pipefail

    if [ "~{defined(pcrplus_sample_ids)}" == "false" ]; then
      touch pcrplus.samples.list
      plus_samples="pcrplus.samples.list"
    else
      plus_samples="~{pcrplus_sample_ids}"
    fi

    if [ "~{defined(mingq_prePost_pcrplus_fails)}" == "false" ]; then
      touch pcrplus.fails.list
      plus_fails="pcrplus.fails.list"
    else
      plus_fails="~{mingq_prePost_pcrplus_fails}"
    fi

    tabix -h ~{vcf} ~{contig} \
    | /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/label_batch_effects.py \
        --unstable-af-pcrminus "~{mingq_prePost_pcrminus_fails}" \
        --unstable-af-pcrplus "$plus_fails" \
        stdin \
        ~{reclassification_table} \
        "$plus_samples" \
        stdout \
        | bgzip -c \
        > "~{prefix}.batch_effects_labeled.vcf.gz"
    tabix -p vcf -f "~{prefix}.batch_effects_labeled.vcf.gz"
  >>>

  output {
    File labeled_vcf = "~{prefix}.batch_effects_labeled.vcf.gz"
    File labeled_vcf_idx = "~{prefix}.batch_effects_labeled.vcf.gz.tbi"
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
