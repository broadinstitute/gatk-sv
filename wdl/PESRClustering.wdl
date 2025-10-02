version 1.0

import "Structs.wdl"
import "TasksClusterBatch.wdl" as tasks_cluster
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as util

workflow ClusterPESR {
  input {
    File vcf_tar

    File ploidy_table
    String batch
    String caller

    # For testing
    File? contig_subset_list
    File? svtk_to_gatk_script
    File? gatk_to_svtk_script

    Int min_size
    File exclude_intervals
    File contig_list

    Float pesr_interval_overlap
    Int pesr_breakend_window
    String? clustering_algorithm

    File reference_fasta
    File reference_fasta_fai
    File reference_dict

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    Float? java_mem_fraction

    RuntimeAttr? runtime_attr_prepare_pesr_vcfs
    RuntimeAttr? runtime_attr_svcluster
    RuntimeAttr? runtime_attr_exclude_intervals_pesr
    RuntimeAttr? runtime_override_concat_vcfs_pesr
    RuntimeAttr? runtime_attr_gatk_to_svtk_vcf
  }

  call PreparePESRVcfs {
    input:
      vcf_tar=vcf_tar,
      reference_fasta_fai=reference_fasta_fai,
      exclude_intervals=exclude_intervals,
      exclude_intervals_index=exclude_intervals + ".tbi",
      ploidy_table=ploidy_table,
      min_size=min_size,
      output_prefix="~{batch}.cluster_batch.~{caller}.prep_vcfs",
      script=svtk_to_gatk_script,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_prepare_pesr_vcfs
  }

  Array[String] contigs = transpose(read_tsv(select_first([contig_subset_list, contig_list])))[0]
  scatter (contig in contigs) {
    call tasks_cluster.SVCluster {
      input:
        vcfs_tar=PreparePESRVcfs.out,
        ploidy_table=ploidy_table,
        output_prefix="~{batch}.cluster_batch.~{caller}.~{contig}.clustered",
        contig=contig,
        fast_mode=true,
        algorithm=clustering_algorithm,
        pesr_sample_overlap=0,
        pesr_interval_overlap=pesr_interval_overlap,
        pesr_breakend_window=pesr_breakend_window,
        reference_fasta=reference_fasta,
        reference_fasta_fai=reference_fasta_fai,
        reference_dict=reference_dict,
        java_mem_fraction=java_mem_fraction,
        variant_prefix="~{batch}_~{caller}_~{contig}_",
        gatk_docker=gatk_docker,
        runtime_attr_override=runtime_attr_svcluster
    }
    call tasks_cluster.ExcludeIntervalsByEndpoints {
      input:
        vcf=SVCluster.out,
        reference_fasta_fai=reference_fasta_fai,
        output_prefix="~{batch}.cluster_batch.~{caller}.~{contig}.exclude_intervals",
        intervals=exclude_intervals,
        intervals_index=exclude_intervals + ".tbi",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_exclude_intervals_pesr
    }
    call tasks_cluster.GatkToSvtkVcf {
      input:
        vcf=ExcludeIntervalsByEndpoints.out,
        output_prefix="~{batch}.cluster_batch.~{caller}.~{contig}.svtk_formatted",
        script=gatk_to_svtk_script,
        source=caller,
        contig_list=contig_list,
        remove_formats="CN",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_gatk_to_svtk_vcf
    }
  }

  call tasks_cohort.ConcatVcfs {
    input:
      vcfs=GatkToSvtkVcf.out,
      vcfs_idx=GatkToSvtkVcf.out_index,
      naive=true,
      outfile_prefix="~{batch}.cluster_batch.~{caller}",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_vcfs_pesr
  }

  output {
    File clustered_vcf = ConcatVcfs.concat_vcf
    File clustered_vcf_index = ConcatVcfs.concat_vcf_idx
  }
}


task PreparePESRVcfs {
  input {
    File vcf_tar
    File reference_fasta_fai
    File exclude_intervals
    File exclude_intervals_index
    File ploidy_table
    Int min_size
    File? script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: 200,
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.tar.gz"
  }
  command <<<
    set -euo pipefail
    cut -f1,2 ~{reference_fasta_fai} > genome.file
    mkdir in/ out/
    tar xzf ~{vcf_tar} -C in/
    i=0
    for VCF in in/*.vcf.gz; do
      NAME=$(basename $VCF .vcf.gz)
      SAMPLE_NUM=`printf %05d $i`
      # Convert format
      python ~{default="/opt/sv-pipeline/scripts/format_svtk_vcf_for_gatk.py" script} \
        --vcf $VCF \
        --out tmp.vcf.gz \
        --ploidy-table ~{ploidy_table}
      # Interval, contig, and size filtering
      bcftools query -f '%CHROM\t%POS\t%POS\t%ID\t%SVTYPE\n%CHROM\t%END\t%END\t%ID\t%SVTYPE\n%CHR2\t%END2\t%END2\t%ID\t%SVTYPE\n' tmp.vcf.gz \
        | awk '$1!="." && $2!="."' \
        | sort -k1,1V -k2,2n -k3,3n \
        > ends.bed
      bedtools intersect -sorted -u -wa -g genome.file -wa -a ends.bed -b ~{exclude_intervals} | cut -f4 | sort | uniq \
        > excluded_vids.list
      bcftools view -i 'ID!=@excluded_vids.list && (INFO/SVLEN="." || INFO/SVLEN=-1 || INFO/SVLEN>=~{min_size})' tmp.vcf.gz \
        -Oz -o out/$SAMPLE_NUM.$NAME.vcf.gz
      tabix out/$SAMPLE_NUM.$NAME.vcf.gz
      i=$((i+1))
    done
    tar czf ~{output_prefix}.tar.gz -C out/ .
  >>>
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