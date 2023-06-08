version 1.0
import "TasksGenotypeBatch.wdl" as tasksgenotypebatch

workflow Regenotype {
  input {
    File depth_vcf
    File regeno_bed
    File cohort_depth_vcf
    File batch_depth_vcf
    File coveragefile
    File coveragefile_idx
    File medianfile
    File? famfile
    File RD_depth_sepcutoff
    Int n_per_split
    Int n_RdTest_bins
    String batch
    File samples_list
    String sv_base_mini_docker
    String sv_pipeline_docker
    String sv_pipeline_rdtest_docker
    Array[String] samples = read_lines(samples_list)

    RuntimeAttr? runtime_attr_add_batch_samples
    RuntimeAttr? runtime_attr_get_regeno_g2
    RuntimeAttr? runtime_attr_split_beds
    RuntimeAttr? runtime_attr_make_subset_vcf
    RuntimeAttr? runtime_attr_rd_test_gt_regeno
    RuntimeAttr? runtime_attr_integrate_depth_gq
    RuntimeAttr? runtime_attr_add_genotypes
    RuntimeAttr? runtime_attr_concat_regenotyped_vcfs_g2
  }
  call tasksgenotypebatch.AddBatchSamples as AddBatchSamplesDepth {
    input:
      batch_vcf=batch_depth_vcf,
      cohort_vcf=cohort_depth_vcf,
      prefix="~{batch}.depth",
      runtime_attr_override = runtime_attr_add_batch_samples,
      sv_pipeline_docker=sv_pipeline_docker
  }
  call GetRegenotype {
    input:
      Batch=batch,
      master_regeno=regeno_bed,
      depth_genotyped_vcf=depth_vcf,
      runtime_attr_override = runtime_attr_get_regeno_g2,
      sv_pipeline_docker=sv_pipeline_docker
  }
  call SplitBeds as SplitBeds_regeno {
    input: 
      bed=GetRegenotype.regeno_bed,
      n_per_split=n_per_split,
      runtime_attr_override = runtime_attr_split_beds,
      sv_pipeline_docker=sv_pipeline_docker
  }
  scatter (regeno in SplitBeds_regeno.regeno_beds) {
    call tasksgenotypebatch.MakeSubsetVcf as make_subset_vcf_regeno {
      input:
        vcf=AddBatchSamplesDepth.updated_vcf,
        bed=regeno,
        runtime_attr_override = runtime_attr_make_subset_vcf,
        sv_base_mini_docker=sv_base_mini_docker
    }
    call RdTestGenotypeRegeno {
      input:
        bed=regeno,
        coveragefile=coveragefile,
        coveragefile_idx=coveragefile_idx,
        generate_melted_genotypes = true,
        medianfile=medianfile,
        famfile=famfile,
        samples=samples,
        gt_cutoffs=RD_depth_sepcutoff,
        n_bins=n_RdTest_bins,
        prefix=basename(regeno, ".bed"),
        runtime_attr_override = runtime_attr_rd_test_gt_regeno,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker
    }
    call tasksgenotypebatch.IntegrateDepthGq as IntegrateGQRegeno {
      input:
        vcf=make_subset_vcf_regeno.subset_vcf,
        RD_melted_genotypes=RdTestGenotypeRegeno.melted_genotypes,
        RD_vargq=RdTestGenotypeRegeno.varGQ,
        runtime_attr_override = runtime_attr_integrate_depth_gq,
        sv_pipeline_docker=sv_pipeline_docker
    }
    call tasksgenotypebatch.AddGenotypes as AddGenotypesRegeno {
      input:
        vcf=make_subset_vcf_regeno.subset_vcf,
        genotypes=IntegrateGQRegeno.genotypes,
        varGQ=IntegrateGQRegeno.varGQ,
        prefix=basename(regeno, ".bed"),
        runtime_attr_override = runtime_attr_add_genotypes,
        sv_pipeline_docker=sv_pipeline_docker
    }
  }
  call ConcatReGenotypedVcfs as ConcatRegenotypedVcfs {
    input:
      batch=batch,
      depth_vcf=depth_vcf,
      regeno_vcfs=AddGenotypesRegeno.genotyped_vcf,
      bed=regeno_bed,
      runtime_attr_override = runtime_attr_concat_regenotyped_vcfs_g2,
      sv_base_mini_docker=sv_base_mini_docker
  }
  output {
    File genotyped_vcf = ConcatRegenotypedVcfs.genotyped_vcf
    File genotyped_vcf_idx = ConcatRegenotypedVcfs.genotyped_vcf_idx
    File regeno_portion=ConcatRegenotypedVcfs.regenotyped_vcf
  }
}
task SplitBeds {
  input {
    File bed
    Int n_per_split
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
    cat ~{bed} \
      | split -l ~{n_per_split} -a 6 - regeno.
  >>>

  output {
    Array[File] regeno_beds = glob("regeno.*")
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

task RdTestGenotypeRegeno {
  input {
    File bed
    File coveragefile
    File coveragefile_idx
    File medianfile
    File? famfile
    Array[String] samples
    File gt_cutoffs
    Int n_bins
    String prefix
    Boolean generate_melted_genotypes
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    coveragefile: {
      localization_optional: true
    }
    coveragefile_idx: {
      localization_optional: true
    }
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
    File genotypes = "~{prefix}.geno"
    File copy_states = "~{prefix}.median_geno"
    File metrics = "~{prefix}.metrics"
    File gq = "~{prefix}.gq"
    File varGQ = "~{prefix}.vargq"
    File melted_genotypes = "rd.geno.cnv.bed.gz"
  }
  command <<<

    set -euo pipefail
    /opt/RdTest/localize_bincov.sh ~{bed} ~{coveragefile}
    Rscript /opt/RdTest/RdTest.R \
      -b ~{bed} \
      -c local_coverage.bed.gz \
      -m ~{medianfile} \
      ~{"-f " + famfile} \
      -n ~{prefix} \
      -v TRUE \
      -w ~{write_lines(samples)} \
      -i ~{n_bins} \
      -r ~{gt_cutoffs} \
      -y /opt/RdTest/bin_exclude.bed.gz \
      -g TRUE;
    if [ ~{generate_melted_genotypes} == "true" ]; then
      /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py ~{prefix}.geno ~{prefix}.gq rd.geno.cnv.bed;
      sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz
    else
      echo "" | bgzip -c > rd.geno.cnv.bed.gz
    fi
  
  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_rdtest_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task GetRegenotype{
  input {
    File depth_genotyped_vcf
    File master_regeno
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
    set -e
    svtk vcf2bed ~{depth_genotyped_vcf} ~{Batch}.bed
    sample=$(fgrep -v "#" ~{Batch}.bed | awk '{if($6!="" )print $6}' | head -n1 | cut -d"," -f1) # a sample in batch is needed for regenotyping
    while read chr start end variant s type; do
      printf "$chr\t$start\t$end\t$variant\t$sample\t$type\n" >> ~{Batch}.to_regeno.bed # modify this to suit input for Rdtest
    done < ~{master_regeno}
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
task ConcatReGenotypedVcfs {
  input {
    String batch
    File depth_vcf
    Array[File] regeno_vcfs
    File bed
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
    File genotyped_vcf = "~{batch}.depth.regeno.vcf.gz"
    File genotyped_vcf_idx = "~{batch}.depth.regeno.vcf.gz.tbi"
    File regenotyped_vcf="~{batch}.regeno.vcf.gz"
  }
  command <<<
    set -euo pipefail
    vcf-concat ~{sep=' ' regeno_vcfs}  \
      | vcf-sort -c \
      | bgzip -c > ~{batch}.regeno.vcf.gz
    cut -f 4 ~{bed} |awk '{print $0"\t"}'> varlist.txt
    zcat ~{depth_vcf} |fgrep -f varlist.txt -v |bgzip -c > no_variant.vcf.gz
    vcf-concat ~{batch}.regeno.vcf.gz no_variant.vcf.gz \
      | vcf-sort -c \
      | bgzip -c > ~{batch}.depth.regeno.vcf.gz
    tabix ~{batch}.depth.regeno.vcf.gz
  
  >>>
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
