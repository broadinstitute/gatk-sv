##########################################################################################

## Base script:   https://portal.firecloud.org/#methods/Talkowski-SV/04_v2_genotype_depth_part2/14/wdl

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

version 1.0

import "Structs.wdl"
import "Tasks04.wdl" as tasks04

workflow GenotypeDepthPart2 {
  input {
    File bin_exclude
    File cohort_vcf
    File RD_pesr_sepcutoff
    File RD_depth_sepcutoff
    Int n_per_split
    Int n_RdTest_bins
    String batch
    File? regeno_sample_counts_lookup # required if doing regenotyping
    File? regeno_raw_combined_depth # required if doing regenotyping
    Int? n_samples_cohort # required if doing regenotyping
    Float regeno_max_allele_freq # default = 0.01 set in Module04.wdl
    Int regeno_allele_count_threshold # default = 3 set in Module04.wdl
    File ref_dict
    File medianfile
    File famfile
    Array[String] samples

    File coveragefile

    String sv_pipeline_docker
    String sv_base_mini_docker
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_rdtest_genotype
    RuntimeAttr? runtime_attr_make_subset_vcf
    RuntimeAttr? runtime_attr_integrate_depth_gq
    RuntimeAttr? runtime_attr_add_genotypes
    RuntimeAttr? runtime_attr_concat_vcfs
  }

  File bin_exclude_idx = bin_exclude + ".tbi"

  call tasks04.SplitVariants as SplitVariants {
  input:
      vcf = cohort_vcf,
      n_per_split = n_per_split,
      generate_bca = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_split_variants
  }

  scatter (gt5kb_bed in SplitVariants.gt5kb_beds) {

    call tasks04.MakeSubsetVcf as MakeSubsetVcfOver5kb {
    input:
    vcf = cohort_vcf,
    bed = gt5kb_bed,
    sv_base_mini_docker = sv_base_mini_docker,
    runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasks04.RDTestGenotype as RDTestGenotypeOver5kb {
    input:
    bin_exclude=bin_exclude,
    bin_exclude_idx=bin_exclude_idx,
    bed = gt5kb_bed,
    coveragefile = coveragefile,
    medianfile = medianfile,
    famfile = famfile,
    samples = samples,
    gt_cutoffs = RD_depth_sepcutoff,
    n_bins = n_RdTest_bins,
    prefix = basename(gt5kb_bed),
    generate_melted_genotypes = true,
    ref_dict = ref_dict,
    sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
    runtime_attr_override = runtime_attr_rdtest_genotype
    }

    call tasks04.IntegrateDepthGq as IntegrateDepthGqOver5kb {
    input:
    vcf = MakeSubsetVcfOver5kb.subset_vcf,
    RD_melted_genotypes = RDTestGenotypeOver5kb.melted_genotypes,
    RD_vargq = RDTestGenotypeOver5kb.varGQ,
    sv_pipeline_docker = sv_pipeline_docker,
    runtime_attr_override = runtime_attr_integrate_depth_gq
    }

    call tasks04.AddGenotypes as AddGenotypesOver5kb {
    input:
    vcf = MakeSubsetVcfOver5kb.subset_vcf,
    genotypes = IntegrateDepthGqOver5kb.genotypes,
    varGQ = IntegrateDepthGqOver5kb.varGQ,
    prefix = basename(gt5kb_bed, ".bed"),
    sv_pipeline_docker = sv_pipeline_docker,
    runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  scatter (lt5kb_bed in SplitVariants.lt5kb_beds) {

    call tasks04.RDTestGenotype as RDTestGenotypeUnder5kb {
    input:
    bin_exclude=bin_exclude,
    bin_exclude_idx=bin_exclude_idx,
    bed = lt5kb_bed,
    coveragefile = coveragefile,
    medianfile = medianfile,
    famfile = famfile,
    samples = samples,
    gt_cutoffs = RD_pesr_sepcutoff,
    n_bins = n_RdTest_bins,
    prefix = basename(lt5kb_bed, ".bed"),
    generate_melted_genotypes = true,
    ref_dict = ref_dict,
    sv_pipeline_rdtest_docker = sv_pipeline_rdtest_docker,
    runtime_attr_override = runtime_attr_rdtest_genotype
    }

    call tasks04.MakeSubsetVcf as MakeSubsetVcfUnder5kb {
    input:
    vcf = cohort_vcf,
    bed = lt5kb_bed,
    sv_base_mini_docker = sv_base_mini_docker,
    runtime_attr_override = runtime_attr_make_subset_vcf
    }

    call tasks04.IntegrateDepthGq as IntegrateDepthGqUnder5kb {
    input:
    vcf = MakeSubsetVcfUnder5kb.subset_vcf,
    RD_melted_genotypes = RDTestGenotypeUnder5kb.melted_genotypes,
    RD_vargq = RDTestGenotypeUnder5kb.varGQ,
    sv_pipeline_docker = sv_pipeline_docker,
    runtime_attr_override = runtime_attr_integrate_depth_gq
    }
    call tasks04.AddGenotypes as AddGenotypesUnder5kb {
    input:
    vcf = MakeSubsetVcfUnder5kb.subset_vcf,
    genotypes = IntegrateDepthGqUnder5kb.genotypes,
    varGQ = IntegrateDepthGqUnder5kb.varGQ,
    prefix = basename(lt5kb_bed, ".bed"),
    sv_pipeline_docker = sv_pipeline_docker,
    runtime_attr_override = runtime_attr_add_genotypes
    }
  }

  call tasks04.ConcatGenotypedVcfs as ConcatGenotypedVcfs {
  input:
      lt5kb_vcfs = AddGenotypesUnder5kb.genotyped_vcf,
      gt5kb_vcfs = AddGenotypesOver5kb.genotyped_vcf,
      bca_vcfs = [],
      batch = batch,
      evidence_type = "depth",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_concat_vcfs
  }
  Boolean regenotype=defined(regeno_sample_counts_lookup) && defined(regeno_raw_combined_depth) && defined(n_samples_cohort)
  if (regenotype){
    call GetRegenotype{input:
      depth_genotyped_vcf=ConcatGenotypedVcfs.genotyped_vcf,
      Batch=batch,
      regeno_sample_counts_lookup = select_first([regeno_sample_counts_lookup]),
      regeno_raw_combined_depth = select_first([regeno_raw_combined_depth]),
      n_samples_cohort = select_first([n_samples_cohort]),
      regeno_max_allele_freq = regeno_max_allele_freq, 
      regeno_allele_count_threshold = regeno_allele_count_threshold, 
      sv_pipeline_docker=sv_pipeline_docker}
    call GetMedianSubset{input: 
      medians=RDTestGenotypeOver5kb.copy_states,
      batch=batch,
      sv_pipeline_docker=sv_pipeline_docker}
    call MedianIntersect{input: 
      median=GetMedianSubset.regeno_median,
      regeno_list=GetRegenotype.regeno_bed,
      sv_pipeline_docker=sv_pipeline_docker}
  }
  output {
    File genotyped_vcf = ConcatGenotypedVcfs.genotyped_vcf
    File genotyped_vcf_index = ConcatGenotypedVcfs.genotyped_vcf_index
    File? regeno_list = MedianIntersect.regeno_bed
  }
}

task IntegrateDepthGq {
  input {
    File vcf
    File RD_melted_genotypes
    File RD_vargq
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

  output {
    File genotypes = "genotype.indiv.depth.txt.gz"
    File varGQ = "genotype.variant.depth.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_depthonly.sh \
      ~{vcf} \
      ~{RD_melted_genotypes} \
      ~{RD_vargq}
  
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
task GetRegenotype{
  input{
    File depth_genotyped_vcf
    File regeno_sample_counts_lookup
    File regeno_raw_combined_depth
    Int n_samples_cohort
    Float regeno_max_allele_freq # default = 0.01 set in Module04.wdl
    Int regeno_allele_count_threshold # default = 3 set in Module04.wdl
    String Batch
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
    }
  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 0,
    max_retries: 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command<<<
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
    output{
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
  input{
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
  command{
    cat ~{sep=' ' medians} |fgrep -v "start"> medians.txt
    R -e 'd<-read.table("medians.txt",header=T,check.names=F);result<-d[((apply(as.matrix(d[,5:length(d)]),FUN=median,1)<0.97 & apply(as.matrix(d[,5:length(d)]),FUN=median,1)>0.85)|(apply(as.matrix(d[,5:length(d)]),FUN=median,1)>1.03 & apply(as.matrix(d[,5:length(d)]),FUN=median,1)<1.65)),1:4];write.table(result,file="~{batch}_to_regeno.bed",sep="\t",quote=F,row.names=F,col.names=F)'    }
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
  input{
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
  command<<<
    set -euo pipefail
    cut -f 4 ~{median} > medmissingID.txt
    { fgrep -w -f medmissingID.txt ~{regeno_list} || true; }|sort -k1,1V -k2,2n -k3,3n -u > median_missing_sorted.bed
    >>>
  output{
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
