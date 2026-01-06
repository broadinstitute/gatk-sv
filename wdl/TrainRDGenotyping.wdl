version 1.0

import "Structs.wdl"
import "TasksGenotypeBatch.wdl" as tasksgenotypebatch

workflow TrainRDGenotyping {
  input {
    File bin_exclude
    File vcf                 # VCF to genotype
    File coveragefile        # batch coverage file
    File? coveragefile_index # batch coverage file
    File medianfile          # batch median file
    File rf_cutoffs          # Random forest cutoffs
    File seed_cutoffs
    Array[String] samples    # List of samples in batch
    String prefix            # prefix use in output files
    Int n_bins               # number of RdTest bins
    Int n_per_split          # number of variants per RdTest split
    String reference_build   #hg19 or hg38
    File ref_dict

    String sv_base_mini_docker
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_training_bed
    RuntimeAttr? runtime_attr_genotype_train
    RuntimeAttr? runtime_attr_generate_cutoff
    RuntimeAttr? runtime_attr_update_cutoff
    RuntimeAttr? runtime_attr_split_variants
    RuntimeAttr? runtime_attr_rdtest_genotype
    RuntimeAttr? runtime_attr_merge_genotypes
  }

  File bin_exclude_idx = bin_exclude + ".tbi"

  call MakeTrainingBed {
    input:
      sample_ID = samples[0],
      reference_build = reference_build,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_training_bed
  }

  call tasksgenotypebatch.RDTestGenotype as GenotypeTrain {
    input:
      bin_exclude=bin_exclude,
      bin_exclude_idx=bin_exclude_idx,
      bed = MakeTrainingBed.bed,
      coveragefile = coveragefile,
      coveragefile_index = coveragefile_index,
      medianfile = medianfile,
      samples = samples,
      gt_cutoffs = seed_cutoffs,
      n_bins = n_bins,
      prefix = prefix,
      generate_melted_genotypes = false,
      ref_dict = ref_dict,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_genotype_train
  }

  call GenerateCutoff {
    input:
      copy_states = GenotypeTrain.copy_states,
      max_copystate = 4,
      prefix = prefix,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_generate_cutoff
  }

  call UpdateCutoff {
    input:
      rf_cutoffs = rf_cutoffs,
      gt_cutoffs = GenerateCutoff.cutoffs,
      prefix = prefix,
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_update_cutoff
  }

  call tasksgenotypebatch.SplitVariants as SplitVariants {
    input:
      vcf = vcf,
      n_per_split = n_per_split,
      generate_bca = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_split_variants
  }

  scatter (pesr_bed in SplitVariants.lt5kb_beds) {
    call tasksgenotypebatch.RDTestGenotype as GenotypePESR {
      input:
        bin_exclude=bin_exclude,
        bin_exclude_idx=bin_exclude_idx,
        bed = pesr_bed,
        coveragefile = coveragefile,
        coveragefile_index = coveragefile_index,
        medianfile = medianfile,
        samples = samples,
        gt_cutoffs = UpdateCutoff.pesr_sepcutoff,
        n_bins = n_bins,
        prefix = basename(pesr_bed),
        generate_melted_genotypes = false,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_rdtest_genotype
    }
  }

  scatter (gt5kb_bed in SplitVariants.gt5kb_beds) {
    call tasksgenotypebatch.RDTestGenotype as GenotypeOver5kb {
      input:
        bin_exclude=bin_exclude,
        bin_exclude_idx=bin_exclude_idx,
        bed = gt5kb_bed,
        coveragefile = coveragefile,
        coveragefile_index = coveragefile_index,
        medianfile = medianfile,
        samples = samples,
        gt_cutoffs = UpdateCutoff.depth_sepcutoff,
        n_bins = n_bins,
        prefix = basename(gt5kb_bed),
        generate_melted_genotypes = false,
        ref_dict = ref_dict,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_rdtest_genotype
    }
  }

  call MergeGenotypeResults {
    input:
      pesr_genotypes = GenotypePESR.genotypes,
      gt5kb_genotypes = GenotypeOver5kb.genotypes,
      pesr_GQ = GenotypePESR.gq,
      gt5kb_GQ = GenotypeOver5kb.gq,
      pesr_varGQ = GenotypePESR.varGQ,
      gt5kb_varGQ = GenotypeOver5kb.varGQ,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_merge_genotypes
  }

  output {
    File genotypes = MergeGenotypeResults.genotypes
    File melted_genotypes = MergeGenotypeResults.melted_genotypes
    File varGQ = MergeGenotypeResults.varGQ
    File GQ = MergeGenotypeResults.GQ
    File pesr_sepcutoff = UpdateCutoff.pesr_sepcutoff
    File depth_sepcutoff = UpdateCutoff.depth_sepcutoff
  }
}

task MakeTrainingBed {
  input {
    String sample_ID
    String reference_build
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
    File bed = "train.bed"
  }
  command <<<

    set -euo pipefail
    if [ ~{reference_build} == "hg19" ]; then
      awk -v OFS="\t" -v sample="~{sample_ID}" '{$5=sample; print $1, $2, $3, $4, $5, $6}' /opt/RdTest/1kg.train.loci.bed > train.bed
    else
      awk -v OFS="\t" -v sample="~{sample_ID}" '{$5=sample; print $1, $2, $3, $4, $5, $6}' /opt/RdTest/train_hg38_reviewed_final.bed > train.bed
    fi
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

task UpdateCutoff {
  input {
    File rf_cutoffs
    File gt_cutoffs
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

  output {
    File pesr_sepcutoff = "~{prefix}.pesr_sepcutoff.txt"
    File depth_sepcutoff = "~{prefix}.depth_sepcutoff.txt"
  }
  command <<<

    set -euo pipefail
    sep=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
                      {if ($col["algtype"]=="PESR" && $col["min_svsize"]==1000 && $col["metric"]=="RD_Median_Separation")
                          print $col["cutoff"]}' \
               ~{rf_cutoffs})
    awk -v var=$sep \
      'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; print; next}
      {
        if ($col["copy_state"]=="1" && $col["cutoffs"]>1-var)
          $col["cutoffs"]=1-var;
        else if ($col["copy_state"]=="2" && $col["cutoffs"]<1+var)
          $col["cutoffs"]=1+var;
        print
      }' \
      ~{gt_cutoffs} \
      | tr ' ' '\t' \
      > ~{prefix}.pesr_sepcutoff.txt;
    sep=$(awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; next}
                      {if ($col["algtype"]=="Depth" && $col["metric"]=="RD_Median_Separation")
                         print $col["cutoff"]}' \
              ~{rf_cutoffs} \
          | sort -nr | head -n 1)
    awk -v var=$sep \
      'NR==1 {for(i=1;i<=NF;i++) col[$i]=i; print; next}
      {
        if ($col["copy_state"]=="1" && $col["cutoffs"]>1-var)
          $col["cutoffs"]=1-var;
        else if ($col["copy_state"]=="2" && $col["cutoffs"]<1+var)
          $col["cutoffs"]=1+var;
        print
      }' \
      ~{gt_cutoffs} \
    | tr ' ' '\t' \
    > ~{prefix}.depth_sepcutoff.txt;
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_base_mini_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

task MergeGenotypeResults {
  input {
    Array[File] pesr_genotypes
    Array[File] gt5kb_genotypes
    Array[File] pesr_GQ
    Array[File] gt5kb_GQ
    Array[File] pesr_varGQ
    Array[File] gt5kb_varGQ
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
    File genotypes = "rd.geno.all"
    File GQ = "rd.GQ.all"
    File melted_genotypes = "rd.geno.cnv.bed.gz"
    File varGQ = "rd.varGQ.all"
  }
  command <<<

    set -euo pipefail
    cat ~{sep=" "  pesr_genotypes} ~{sep=" "  gt5kb_genotypes} | awk '!_[$0]++' > rd.geno.all;
    cat ~{sep=" "  pesr_GQ} ~{sep=" "  gt5kb_GQ} | awk '!_[$0]++' > rd.GQ.all;
    cat ~{sep=" "  pesr_varGQ} ~{sep=" "  gt5kb_varGQ} | awk '!_[$0]++' > rd.varGQ.all;
    
    /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py rd.geno.all rd.GQ.all rd.geno.cnv.bed;
    sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}

task GenerateCutoff {
  input {
    File copy_states
    Int max_copystate
    String prefix
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
    File cutoffs = "${prefix}.cutoffs"
  }
  command <<<

    Rscript /opt/RdTest/generate_cutoff.R ~{copy_states} ~{max_copystate} ~{prefix}.cutoffs
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    noAddress: true
  }
}
