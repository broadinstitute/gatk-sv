##########################
## EXPERIMENTAL WORKFLOW
##########################

version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks


workflow FilterCleanupQualRecalibration {
  input{
    File vcf
    File vcf_idx
    File? pcrplus_samples_list
    File famfile
    Float min_callrate_global
    Float min_callrate_smallDels
    File contiglist
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_ConcatVcfs
  }
  Array[Array[String]] contigs = read_tsv(contiglist)

  scatter ( contig in contigs ) {
    call Cleanup {
      input:
        vcf=vcf,
        vcf_idx=vcf_idx,
        contig=contig[0],
        pcrplus_samples_list=pcrplus_samples_list,
        famfile=famfile,
        min_callrate_global=min_callrate_global,
        min_callrate_smallDels=min_callrate_smallDels,
        prefix=prefix,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }

  call MiniTasks.ConcatVcfs as ConcatVcfs {
    input:
      vcfs=Cleanup.out_vcf,
      naive=true,
      outfile_prefix="~{prefix}.cleaned_filters_qual_recalibrated",
      sv_base_mini_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_ConcatVcfs
  }

  output {
    File cleaned_vcf = ConcatVcfs.concat_vcf
    File cleaned_vcf_idx =ConcatVcfs.concat_vcf_idx
  }
}


# Applies filters & cleanup to VCF for a single chromosome
task Cleanup {
  input{
    File vcf
    File vcf_idx
    String contig
    File? pcrplus_samples_list
    File famfile
    Float min_callrate_global
    Float min_callrate_smallDels
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    
    set -euo pipefail
    #Subset to chromosome of interest
    tabix -h ~{vcf} ~{contig} | bgzip -c > input.vcf.gz
    #Get list of PCR- samples
    tabix -H ~{vcf} | fgrep -v "##" | cut -f10- | sed 's/\t/\n/g' \
    > all.samples.list
    if [ ! -z "~{pcrplus_samples_list}" ];then
      fgrep -wvf ~{pcrplus_samples_list} all.samples.list \
      > pcrminus.samples.list
    else
      cp all.samples.list pcrminus.samples.list
    fi
    #Restrict famfiles
    awk -F "\t" 'NR==FNR{c[$1]++;next};c[$2] > 0' all.samples.list ~{famfile}  > revised.pre
    awk 'NR==FNR{o[FNR]=$1; next} {t[$2]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' all.samples.list revised.pre > revised.fam
    fgrep -wf pcrminus.samples.list revised.fam > revised.pcrminus.fam
    #Compute fraction of missing genotypes per variant
    zcat input.vcf.gz \
    | awk '{ if ($7 !~ /MULTIALLELIC/) print $0 }' \
    | bgzip -c \
    > input.noMCNV.vcf.gz
    plink2 \
      --missing variant-only \
      --max-alleles 2 \
      --keep-fam revised.pcrminus.fam \
      --fam revised.fam \
      --vcf input.noMCNV.vcf.gz
    fgrep -v "#" plink2.vmiss \
    | awk -v OFS="\t" '{ print $2, 1-$NF }' \
    > callrates.txt
    #Clean up VCF
    /opt/sv-pipeline/scripts/downstream_analysis_and_filtering/filter_cleanup_and_QUAL_recalibration.PCRMinus_only.py \
      --callrate-table callrates.txt \
      --min-callrate-global ~{min_callrate_global} \
      --min-callrate-smallDels ~{min_callrate_smallDels} \
      input.vcf.gz \
      stdout \
    | bgzip -c \
    > "~{prefix}.~{contig}.cleaned_filters_qual_recalibrated.vcf.gz"
  >>>

  output {
    File out_vcf = "~{prefix}.~{contig}.cleaned_filters_qual_recalibrated.vcf.gz"
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

