version 1.0

import "Genotype_2.wdl" as g2
import "CombineReassess.wdl" as creassess
import "Genotype_3.wdl" as g3

workflow Module04b{
  input{
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_base_docker
    String sv_pipeline_rdtest_docker
    Array[File] depth_vcfs
    File cohort_depth_vcf
    Array[File] batch_depth_vcfs
    Array[File] coveragefiles
    Array[File] coveragefile_idxs
    Array[File] medianfiles
    Array[File] famfiles
    Array[File] RD_depth_sepcutoffs
    Int n_per_split
    String sv_base_mini_docker
    Int n_RdTest_bins
    Array[String] batches
    Array[File] samples_lists
    File regeno_sample_ids_lookup
    File regeno_sample_counts_lookup 
    File regeno_raw_combined_depth 
    Array[File] regeno_coverage_medians
    Float regeno_max_allele_freq = 0.01 
    Int regeno_allele_count_threshold = 3 
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_merge_list
    RuntimeAttr? runtime_attr_get_count_cohort_samplelist
    RuntimeAttr? runtime_attr_get_regeno
    RuntimeAttr? runtime_attr_get_median_subset
    RuntimeAttr? runtime_attr_median_intersect
  }
  call GetAndCountCohortSampleList {
    input:
      batch_sample_lists = samples_lists,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_get_count_cohort_samplelist
  }
  scatter(i in range(length(batches))) {
    call GetRegenotype {
      input:
        depth_genotyped_vcf = depth_vcfs[i],
        Batch = batches[i],
        regeno_sample_counts_lookup = regeno_sample_counts_lookup,
        regeno_raw_combined_depth = regeno_raw_combined_depth,
        n_samples_cohort = GetAndCountCohortSampleList.n_samples_cohort,
        regeno_max_allele_freq = regeno_max_allele_freq, 
        regeno_allele_count_threshold = regeno_allele_count_threshold, 
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_get_regeno
    }
    call GetMedianSubset {
      input: 
        medians = regeno_coverage_medians,
        batch = batches[i],
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_get_median_subset
    }
    call MedianIntersect {
      input: 
        median = GetMedianSubset.regeno_median,
        regeno_list = GetRegenotype.regeno_bed,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_median_intersect
    }
  }
  call MergeList{
    input:
      regeno_beds = MedianIntersect.regeno_bed,
      prefix="master_regeno",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_merge_list
  }
  scatter (i in range(length(batches))) {
    call g2.Regenotype as Genotype_2 {
      input:
        depth_vcf=depth_vcfs[i],
        regeno_bed= MergeList.master_regeno,
        cohort_depth_vcf=cohort_depth_vcf,
        batch_depth_vcf=batch_depth_vcfs[i],
        coveragefile=coveragefiles[i],
        coveragefile_idx=coveragefile_idxs[i],
        medianfile=medianfiles[i],
        famfile=famfiles[i],
        RD_depth_sepcutoff=RD_depth_sepcutoffs[i],
        n_per_split=n_per_split,
        n_RdTest_bins=n_RdTest_bins,
        batch=batches[i],
        samples_list=samples_lists[i],
        sv_pipeline_docker=sv_pipeline_docker,
        sv_base_mini_docker=sv_base_mini_docker,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker
    }
  }

  call creassess.CombineReassess as CombineReassess {
    input:
      samplelist=GetAndCountCohortSampleList.cohort_samplelist,
      regeno_file=MergeList.master_regeno,
      regeno_sample_ids_lookup=regeno_sample_ids_lookup,
      vcfs=Genotype_2.genotyped_vcf,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_pipeline_base_docker=sv_pipeline_base_docker,
      runtime_attr_vcf2bed = runtime_attr_vcf2bed
  }
    
  scatter (i in range(length(Genotype_2.genotyped_vcf))) {
    call ConcatRegenotypedVcfs{
      input:
        depth_vcf=depth_vcfs[i],
        batch=batches[i],
        regeno_vcf=Genotype_2.genotyped_vcf[i],
        regeno_variants=CombineReassess.regeno_variants,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }
  output{
    Array[File] regenotyped_depth_vcfs = ConcatRegenotypedVcfs.genotyped_vcf
  }
}

task GetAndCountCohortSampleList {
  input {
    Array[File] batch_sample_lists
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
    set -euxo pipefail
    # concatenate batch sample lists (one sample ID per line) into cohort sample list
    while read batch_list; do
      cat $batch_list >> cohort.sample.list
    done < ~{write_lines(batch_sample_lists)}
    # count lines = number of samples
    wc -l < cohort.sample.list | tr -d ' ' > cohort_sample_count.txt
  >>>
  output {
    File cohort_samplelist = "cohort.sample.list"
    File cohort_samplecount_file = "cohort_sample_count.txt"
    Int n_samples_cohort = read_int("cohort_sample_count.txt")
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

task GetRegenotype{
  input {
    File depth_genotyped_vcf
    File regeno_sample_counts_lookup
    File regeno_raw_combined_depth
    Int n_samples_cohort
    Float regeno_max_allele_freq # default = 0.01 set in Module04b.wdl
    Int regeno_allele_count_threshold # default = 3 set in Module04b.wdl
    String Batch
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
    set -euxo pipefail
    n_samp=~{n_samples_cohort}
    max_af=~{regeno_max_allele_freq}
    min_count=~{regeno_allele_count_threshold}
    svtk vcf2bed ~{depth_genotyped_vcf} ~{Batch}.bed
    sample=$(fgrep -v "#" ~{Batch}.bed|awk '{if($6!="" )print $6}' |head -n1|cut -d"," -f1)||true
    fgrep "~{Batch}"_ ~{Batch}.bed|sort -k4,4> ~{Batch}.origin.depth.bed
    sort -k4,4 ~{regeno_raw_combined_depth} > sort.bed
    join ~{Batch}.origin.depth.bed sort.bed -1 4 -2 4 -o 1.1,1.2,1.3,1.4,1.5,1.6,2.6 > f3.txt
    while read line;do
        string=$(echo "$line" | cut -d" " -f6 |sed 's/,/|/g')
        string2=$(echo "$line" | cut -d" " -f7 |sed 's/,/|/g')
        echo "$line" |sed -r "s/$string|,//g" |awk -v OFS="\t" '{if ($3-$2>10000 && $6!="" && $5!="CN0")print $1,$2,$3,$4,$5}' >> missing_RF.bed
        echo "$line" |sed -r "s/$string2|,//g" |awk -v OFS="\t" '{if ($3-$2>10000 && $6!="" && $5!="CN0")print $1,$2,$3,$4,$5}' >> missing_RF.bed
    done<f3.txt
    sort -u missing_RF.bed > test.txt;mv test.txt missing_RF.bed
    # VF<0.01 (or provided value) filter for missing sample variants, or AC<=3 (or provided value) if cohort is small
    while read chr start end variant type; do
     varid=$variant: #varid should be single variants
     num=$(fgrep -m1 $varid ~{regeno_sample_counts_lookup}|cut -f8)
     if python -c "exit(0 if (((float($num) / float($n_samp)) < float($max_af)) or (int($num) <= int($min_count))) else 1)"; then
        printf "$chr\t$start\t$end\t$variant\t$sample\t$type\n">>~{Batch}.to_regeno.bed # modify this to suit intput for Rdtest
     fi
    done < missing_RF.bed 
  >>>
  output {
    File regeno_bed="~{Batch}.to_regeno.bed"
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

task GetMedianSubset{
  input {
    String batch
    Array[File] medians
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
  command {
    cat ~{sep=' ' medians} |fgrep -v "start"> medians.txt
    R -e 'd<-read.table("medians.txt",header=T,check.names=F);result<-d[((apply(as.matrix(d[,5:length(d)]),FUN=median,1)<0.97 & apply(as.matrix(d[,5:length(d)]),FUN=median,1)>0.85)|(apply(as.matrix(d[,5:length(d)]),FUN=median,1)>1.03 & apply(as.matrix(d[,5:length(d)]),FUN=median,1)<1.65)),1:4];write.table(result,file="~{batch}_to_regeno.bed",sep="\t",quote=F,row.names=F,col.names=F)'
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
  output {
    File regeno_median = "~{batch}_to_regeno.bed"
    File median="medians.txt"
  }
}

task MedianIntersect{
  input {
    File median
    File regeno_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 2
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    cut -f 4 ~{median} > medmissingID.txt
    { fgrep -w -f medmissingID.txt ~{regeno_list} || true; }|sort -k1,1V -k2,2n -k3,3n -u > median_missing_sorted.bed
  >>>
  output {
    File regeno_bed="median_missing_sorted.bed"
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

task MergeList {
  input {
    String prefix
    Array[File] regeno_beds
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

  command {
    cat ~{sep=' ' regeno_beds} |sort -k1,1V -k2,2n -k3,3n > ~{prefix}.bed
  }
  output {
    File master_regeno="~{prefix}.bed"
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

task ConcatRegenotypedVcfs {
  input {
    String batch
    File regeno_variants
    File depth_vcf
    File regeno_vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 16,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euo pipefail
    zcat ~{regeno_vcf} |fgrep "#" > head.txt
    zcat ~{regeno_vcf} |fgrep -f ~{regeno_variants} >body.txt
    cat head.txt body.txt|bgzip -c > regeno.vcf.gz
    zcat ~{depth_vcf} |fgrep -f ~{regeno_variants} -v |bgzip -c > no_variant.vcf.gz
    vcf-concat regeno.vcf.gz no_variant.vcf.gz \
      | vcf-sort -c \
      | bgzip -c > ~{batch}.depth.regeno_final.vcf.gz
    tabix ~{batch}.depth.regeno_final.vcf.gz
  >>>
  output {
    File genotyped_vcf = "~{batch}.depth.regeno_final.vcf.gz"
    File genotyped_vcf_idx = "~{batch}.depth.regeno_final.vcf.gz.tbi"
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
