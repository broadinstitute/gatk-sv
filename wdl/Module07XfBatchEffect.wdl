##########################
## EXPERIMENTAL WORKFLOW
##########################

version 1.0

import "CalcAF.wdl" as calcAF
import "batch_effect_helper.wdl" as helper
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as util

workflow XfBatchEffect {
  input {
    File vcf
    File vcf_idx
    File sample_batch_assignments
    File batches_list
    Array[File] batch_include_lists
    File sample_pop_assignments
    File famfile
    File? par_bed
    Int? onevsall_cutoff=2
    String prefix
    File? af_pcrmins_premingq

    Int records_per_shard_compare
    Int records_per_shard_calc_af
    Int records_per_shard_label
    String sv_pipeline_docker
    String sv_base_mini_docker

    RuntimeAttr? runtime_attr_merge_labeled_vcfs
  }
  Array[String] batches = read_lines(batches_list)

  # Shard VCF per batch, compute pops-specific AFs, and convert to table of VID & AF stats
  scatter ( i in range(length(batches)) ) {

    call util.SubsetVcfBySamplesList {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        outfile_name = "~{prefix}.~{batches[i]}.subset.vcf.gz",
        list_of_samples = batch_include_lists[i],
        remove_private_sites = false,
        sv_base_mini_docker = sv_base_mini_docker
    }
    # Prune VCF to samples
    call calcAF.CalcAF {
      input:
        vcf=SubsetVcfBySamplesList.vcf_subset,
        vcf_idx=SubsetVcfBySamplesList.vcf_subset_index,
        prefix="~{prefix}.~{batches[i]}",
        sample_pop_assignments=sample_pop_assignments,
        famfile=famfile,
        drop_empty_records="FALSE",
        sv_per_shard=records_per_shard_calc_af,
        par_bed=par_bed,
        sv_pipeline_docker=sv_pipeline_docker,
        sv_pipeline_updates_docker=sv_pipeline_docker
    }
    # Get minimal table of AF data per batch, split by ancestry
    call GetFreqTable {
      input:
        vcf=CalcAF.vcf_wAFs,
        sample_pop_assignments=sample_pop_assignments,
        prefix=batches[i],
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  # Merge frequency results per batch into a single table of all variants with AF data across batches
  call MergeFreqTables {
    input:
      tables=GetFreqTable.freq_data,
      batches_list=batches_list,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }
  call MergeFreqTables as MergeFreqTables_allPops {
    input:
      tables=GetFreqTable.freq_data_allPops,
      batches_list=batches_list,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Compare frequencies before and after minGQ, and generate list of variants
  # that are significantly different between the steps
  if (defined(af_pcrmins_premingq)) {
    call CompareFreqsPrePostMinGQPcrminus {
      input:
        af_pcrmins_premingq=select_first([af_pcrmins_premingq]),
        AF_postMinGQ_table=MergeFreqTables_allPops.merged_table,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  # Perform one-vs-all comparison of AFs per batch to find batch-specific sites
  scatter ( batch in batches ) {
    call helper.check_batch_effects as one_vs_all_comparison {
      input:
        freq_table=MergeFreqTables.merged_table,
        batch1=batch,
        batch2="ALL_OTHERS",
        prefix=prefix,
        variants_per_shard=records_per_shard_compare,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }
  # Collect results from pairwise batch effect detection
  call MergeVariantFailureLists as merge_one_vs_all_checks {
    input:
      fail_variant_lists=one_vs_all_comparison.batch_effect_variants,
      prefix="~{prefix}.one_vs_all_comparisons",
      sv_pipeline_docker=sv_pipeline_docker
  }

  # Distill final table of variants to be reclassified
  call MakeReclassificationTable {
    input:
      freq_table=MergeFreqTables.merged_table,
      onevsall_fails=merge_one_vs_all_checks.fails_per_variant_all,
      prefix=prefix,
      onevsall_cutoff = onevsall_cutoff,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call MiniTasks.ScatterVcf {
    input:
      vcf=vcf,
      prefix=prefix,
      records_per_shard=records_per_shard_label,
      sv_pipeline_docker=sv_pipeline_docker
  }

  scatter (i in range(length(ScatterVcf.shards))) {
    call ApplyBatchEffectLabels {
      input:
        vcf=ScatterVcf.shards[i],
        reclassification_table=MakeReclassificationTable.reclassification_table,
        mingq_prePost_pcrminus_fails=CompareFreqsPrePostMinGQPcrminus.pcrminus_fails,
        prefix="~{prefix}.shard_~{i}",
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  call MiniTasks.ConcatVcfs as merge_labeled_vcfs {
    input:
      vcfs=ApplyBatchEffectLabels.labeled_vcf,
      naive=true,
      outfile_prefix="~{prefix}.batch_effects_labeled_merged",
      sv_base_mini_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_merge_labeled_vcfs
  }

  output {
    File labeled_vcf = merge_labeled_vcfs.concat_vcf
    File labeled_vcf_idx = merge_labeled_vcfs.concat_vcf_idx
    File reclassification_table = MakeReclassificationTable.reclassification_table
  }
}


# Get list of samples to include & exclude per batch
# Always exclude probands from all batches
task GetBatchSamplesList {
  input {
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
  input {
    File vcf
    File sample_pop_assignments
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 6,
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
      | grep -e 'name\|SVLEN\|SVTYPE\|_AC\|_AN\|AC_\|AN_\|_CN_NONREF_COUNT\|_CN_NUMBER' \
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
  input {
    Array[File] tables
    File batches_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 16,
                               disk_gb: 100,
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
  >>>

  output {
    File merged_table = "~{prefix}.merged_AF_table.txt.gz"
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


# Compare
task CompareFreqsPrePostMinGQPcrminus {
  input {
    File af_pcrmins_premingq
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
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/compare_freqs_pre_post_minGQ.PCRMinus_only.R \
      ~{af_pcrmins_premingq} \
      ~{AF_postMinGQ_table} \
      ./ \
      "~{prefix}."
  >>>

  output {
    File pcrminus_fails = "~{prefix}.PCRMINUS_minGQ_AF_prePost_fails.VIDs.list"
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


# Calculate & plot cross-batch correlation coefficient matrixes
task MakeCorrelationMatrices {
  input {
    File freq_table
    String pop
    File batches_list
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 8,
                               disk_gb: 50,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/correlate_batches_singlePop.R \
      ~{batches_list} \
      ~{freq_table} \
      "~{pop}" \
      "~{prefix}.~{pop}"
  >>>
  output {
    Array[File] corr_matrixes = glob("~{prefix}.~{pop}.*.R2_matrix.txt")
    Array[File] heat_maps = glob("~{prefix}.~{pop}.*heatmap*.pdf")
    Array[File] dot_plots = glob("~{prefix}.~{pop}.*perBatch_R2_sina_plot.pdf")
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

# Merge lists of batch effect checks and count total number of times each variant failed
task MergeVariantFailureLists {
  input {
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
  input {
    File freq_table
    File onevsall_fails
    String prefix
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
      ~{onevsall_fails} \
      "~{prefix}.batch_effect_reclassification_table.txt" \
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
  input {
    File vcf
    File reclassification_table
    File? mingq_prePost_pcrminus_fails
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
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/label_batch_effects.PCRMinus_only.py \
      ~{vcf} \
      ~{reclassification_table} \
      ~{"--unstable-af-pcrminus " + mingq_prePost_pcrminus_fails} \
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