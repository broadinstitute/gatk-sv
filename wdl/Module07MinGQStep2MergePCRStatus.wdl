# Author: Ryan Collins <rlcollins@g.harvard.edu>

# This is an analysis WDL that wraps three steps in the Talkowski SV pipeline:
#  1) minGQ optimization
#  2) minGQ filter application
#  3) post-filter VCF QC

# This is the second build of this workflow, which enumerates many more fine-grained
#  minGQ filtering conditions, and may not be optimized for small cohorts with fewer
#  variants

version 1.0

import "MasterVcfQc.wdl" as QC


workflow MinGQStep2MergePcrStatus {
  input{
    String prefix
    File vcf_pcrplus
    File vcf_pcrplus_idx
    File vcf_pcrminus
    File vcf_pcrminus_idx
    Array[File]? thousand_genomes_benchmark_calls
    Array[File]? hgsv_benchmark_calls
    Array[File]? asc_benchmark_calls
    File? sanders_2015_tarball
    File? collins_2017_tarball
    File? werling_2018_tarball
    Int? random_seed

    File trios_famfile
    File contiglist

    String sv_pipeline_qc_docker
    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_override_collect_vids_per_sample
  }

  Array[String] contigs = transpose(read_tsv(contiglist))[0]


  scatter (i in range(length(contigs))) {

    #split PCR Minus vcf into each contig
    call Split_Vcf_by_contig as split_pcr_plus{
      input:
        vcf = vcf_pcrplus,
        vcf_idx = vcf_pcrplus_idx,
        contig = contigs[i],
        sv_pipeline_docker = sv_pipeline_docker
    }

    call Split_Vcf_by_contig as split_pcr_minux{
      input:
        vcf = vcf_pcrminus,
        vcf_idx = vcf_pcrminus_idx,
        contig = contigs[i],
        sv_pipeline_docker = sv_pipeline_docker
    }


    # Merge filtered VCFs by PCR status & across chromosomes
    call merge_PCR_VCFs {
      input:
        PCRPLUS_vcf=split_pcr_plus.vcf_out,
        PCRMINUS_vcf=split_pcr_minux.vcf_out,
        prefix=prefix
    }
  }

  call combine_vcfs {
    input:
      vcfs=merge_PCR_VCFs.merged_vcf,
      prefix=prefix
  }


  # Run QC on filtered VCF
  call QC.MasterVcfQc as filtered_VCF_QC {
    input:
      vcf=combine_vcfs.vcf,
      vcf_idx=combine_vcfs.vcf_idx,
      ped_file=trios_famfile,
      prefix="${prefix}",
      sv_per_shard=10000,
      samples_per_shard=100,
      thousand_genomes_benchmark_calls=thousand_genomes_benchmark_calls,
      hgsv_benchmark_calls=hgsv_benchmark_calls,
      asc_benchmark_calls=asc_benchmark_calls,
      sanders_2015_tarball=sanders_2015_tarball,
      collins_2017_tarball=collins_2017_tarball,
      werling_2018_tarball=werling_2018_tarball,
      contigs=contigs,
      random_seed=random_seed,
      sv_pipeline_qc_docker=sv_pipeline_qc_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_override_collect_vids_per_sample=runtime_override_collect_vids_per_sample
  }

  # Final output
  output {
    File filtered_VCF = combine_vcfs.vcf
    File filtered_VCF_idx = combine_vcfs.vcf_idx
    File filtered_VCF_QC_output = filtered_VCF_QC.sv_vcf_qc_output
  }
}


# Get lists of PCRPLUS and PCRMINUS samples present in input VCF
task get_sample_lists {
  input{
    File vcf
    File vcf_idx
    File PCRPLUS_samples_list
    String prefix
  }

  command <<<
    set -euo pipefail
    tabix -H ~{vcf} | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' > all_samples.list
    fgrep -wf ~{PCRPLUS_samples_list} all_samples.list > "~{prefix}.PCRPLUS.samples.list" || true
    fgrep -wvf ~{PCRPLUS_samples_list} all_samples.list > "~{prefix}.PCRMINUS.samples.list" || true
    cat \
      <( awk -v OFS="\t" '{ print $1, "PCRPLUS" }' "~{prefix}.PCRPLUS.samples.list" || true ) \
      <( awk -v OFS="\t" '{ print $1, "PCRMINUS" }' "~{prefix}.PCRMINUS.samples.list" || true ) \
    > "~{prefix}.PCR_status_assignments.txt"
  >>>

  output {
    File updated_PCRPLUS_samples_list = "~{prefix}.PCRPLUS.samples.list"
    File updated_PCRMINUS_samples_list = "~{prefix}.PCRMINUS.samples.list"
    File sample_PCR_labels = "~{prefix}.PCR_status_assignments.txt"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 50 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


# Split a VCF into two parts, corresponding to PCR+ and PCR-
task split_PCR_vcf {
  input{
    File vcf
    String prefix
    File PCRPLUS_samples_list
  }

  command <<<
    set -euo pipefail
    #Get index of PCR+ samples
    PCRPLUS_idxs=$( zcat ~{vcf} | sed -n '1,500p' | fgrep "#" | fgrep -v "##" \
                    | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print NR, $1 }' \
                    | fgrep -wf ~{PCRPLUS_samples_list} | cut -f1 | paste -s -d, )
    #Get PCR+ VCF
    zcat ~{vcf} \
    | cut -f1-9,"$PCRPLUS_idxs" \
    | bgzip -c \
    > "~{prefix}.PCRPLUS.vcf.gz"
    tabix -f "~{prefix}.PCRPLUS.vcf.gz"
    #Get PCR- VCF
    zcat ~{vcf} \
    | cut --complement -f"$PCRPLUS_idxs" \
    | bgzip -c \
    > "~{prefix}.PCRMINUS.vcf.gz"
    tabix -f "~{prefix}.PCRMINUS.vcf.gz"
  >>>

  output {
    File PCRPLUS_vcf = "~{prefix}.PCRPLUS.vcf.gz"
    File PCRPLUS_vcf_idx = "~{prefix}.PCRPLUS.vcf.gz.tbi"
    File PCRMINUS_vcf = "~{prefix}.PCRMINUS.vcf.gz"
    File PCRMINUS_vcf_idx = "~{prefix}.PCRMINUS.vcf.gz.tbi"
  }

  runtime {
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 50 HDD"
    preemptible: 1
    maxRetries: 1
  }
}


task Split_Vcf_by_contig {
  input{
    File vcf
    File vcf_idx
    String contig
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
      set -euo pipefail
      #Tabix chromosome of interest
      tabix -h ~{vcf} ~{contig} | bgzip -c > ~{contig}.vcf.gz
      tabix -p vcf ~{contig}.vcf.gz
  >>>

  output {
    File vcf_out = "~{contig}.vcf.gz"
    File vcf_out_index = "~{contig}.vcf.gz.tbi"
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
task merge_PCR_VCFs {
  input{
    File PCRPLUS_vcf
    File PCRMINUS_vcf
    String prefix
  }

  command <<<
    #Sanitize FILTER columns
    zcat ~{PCRPLUS_vcf} | cut -f7 | grep -ve '^#' | sed '1d' > PCRPLUS_filters.txt
    zcat ~{PCRMINUS_vcf} | cut -f7 | grep -ve '^#' | sed '1d' > PCRMINUS_filters.txt
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/merge_filter_columns.py \
      PCRPLUS_filters.txt \
      PCRMINUS_filters.txt \
      merged_filters.txt
    #Write new VCF header
    zcat ~{PCRPLUS_vcf} | sed -n '1,1000p' | grep -e '^##' > "~{prefix}.minGQ_filtered.vcf"
    zcat ~{PCRMINUS_vcf} | sed -n '1,1000p' | grep -e '^##' | fgrep "NOCALL_RATE" >> "~{prefix}.minGQ_filtered.vcf"
    #Column-wise merger
    paste \
      <( zcat ~{PCRPLUS_vcf} | grep -ve '^##' | cut -f1-6 ) \
      <( cat <( echo -e "FILTER" ) merged_filters.txt ) \
      <( zcat ~{PCRPLUS_vcf} | grep -ve '^##' | cut -f8- ) \
      <( zcat ~{PCRMINUS_vcf} | grep -ve '^##' | cut -f10- ) \
    >> "~{prefix}.minGQ_filtered.vcf"
    /opt/sv-pipeline/scripts/drop_empty_records.py \
      "~{prefix}.minGQ_filtered.vcf" \
      "~{prefix}.minGQ_filtered.no_blanks.vcf"

    #extract and add MCNVs:
    awk '{if ($7=="MULTIALLELIC") print}'  "~{prefix}.minGQ_filtered.vcf" >> "~{prefix}.minGQ_filtered.no_blanks.vcf"

    #Bgzip & tabix
    vcf-sort  "~{prefix}.minGQ_filtered.no_blanks.vcf" | bgzip > "~{prefix}.minGQ_filtered.no_blanks.vcf.gz"
  >>>

  output {
    File merged_vcf = "~{prefix}.minGQ_filtered.no_blanks.vcf.gz"
  }

  runtime {
    preemptible: 1
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 250 SSD"
    memory: "4 GB"
  }
}

# Merge per-chromosome VCF shards
task combine_vcfs {
  input{
    Array[File] vcfs
    String prefix

  }
  
  command <<<
    vcf-concat ~{sep=" " vcfs} | vcf-sort | bgzip -c > ~{prefix}.minGQ_filtered.vcf.gz;
    tabix -p vcf ~{prefix}.minGQ_filtered.vcf.gz
  >>>

  runtime {
    preemptible: 0
    maxRetries: 1
    docker: "talkowski/sv-pipeline@sha256:193d18c26100fdd603c569346722513f5796685e990ec3abcaeb4be887062a1a"
    disks: "local-disk 250 SSD"
    memory: "4 GB"
  }

  output {
    File vcf="~{prefix}.minGQ_filtered.vcf.gz"
    File vcf_idx="~{prefix}.minGQ_filtered.vcf.gz.tbi"
  }
}





