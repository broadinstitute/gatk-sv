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

  call RemoveMCNVs{
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      sv_pipeline_docker=sv_pipeline_docker
  }

  scatter ( contig in contigs ) {
    call Cleanup {
      input:
        vcf=RemoveMCNVs.no_mcnv_vcf,
        vcf_idx=RemoveMCNVs.no_mcnv_idx,
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
      outfile_prefix="~{prefix}.cleaned_filters_qual_recalibrated",
      naive=true,
      sv_base_mini_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_ConcatVcfs
  }

  call MiniTasks.ConcatVcfs as MergeMCNV {
    input:
      vcfs= [ConcatVcfs.concat_vcf, RemoveMCNVs.mcnv_vcf],
      vcfs_idx = [ConcatVcfs.concat_vcf_idx,RemoveMCNVs.mcnv_idx],
      allow_overlaps = true,
      outfile_prefix = "~{prefix}.cleaned_filters_qual_recali",
      sv_base_mini_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_ConcatVcfs
      }

    output {
      File cleaned_vcf = MergeMCNV.concat_vcf
      File cleaned_vcf_idx = MergeMCNV.concat_vcf_idx
    }
}

#remove mCNV from the vcf, which will be added back to the output:

task RemoveMCNVs{
  input{
    File vcf
    File vcf_idx
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
    set -euo pipefail
    zcat ~{vcf} | awk '{if ($7!="MULTIALLELIC") print}' | bgzip > no_MCNV.vcf.gz
    tabix no_MCNV.vcf.gz
    zcat ~{vcf} | grep '#' > MCNV.vcf
    zcat ~{vcf} | awk '{if ($7=="MULTIALLELIC") print}' >> MCNV.vcf
    bgzip MCNV.vcf
    tabix MCNV.vcf.gz
  >>>

  output{
    File mcnv_vcf = "MCNV.vcf.gz"
    File mcnv_idx = "MCNV.vcf.gz.tbi"
    File no_mcnv_vcf = "no_MCNV.vcf.gz"
    File no_mcnv_idx = "no_MCNV.vcf.gz.tbi"
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

task MergeMCNV{
  input{
    File vcf
    File vcf_idx
    File mcnv
    File mcnv_idx
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
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    #zcat ~{mcnv} |uniq | bgzip > mcnv.vcf.gz
    #tabix mcnv.vcf.gz
    vcf-concat ~{vcf} ~{mcnv} | vcf-sort | bgzip > ~{prefix}.cleaned_filters_qual_recali.vcf.gz
    tabix ~{prefix}.cleaned_filters_qual_recali.vcf.gz
  >>>

  output{
    File with_mcnv_vcf = "~{prefix}.cleaned_filters_qual_recali.vcf.gz"
    File with_mcnv_idx = "~{prefix}.cleaned_filters_qual_recali.vcf.gz.tbi"
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
    #while read ptn; do fgrep -w $ptn ~{famfile}; done < all.samples.list > revised.fam
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
    # tabix -p vcf -f "~{prefix}.cleaned_filters_qual_recalibrated.vcf.gz"
  >>>

  output {
    File out_vcf = "~{prefix}.~{contig}.cleaned_filters_qual_recalibrated.vcf.gz"
    # File out_vcf_idx = "~{prefix}.cleaned_filters_qual_recalibrated.vcf.gz.tbi"
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

