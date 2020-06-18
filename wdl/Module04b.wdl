version 1.0

import "Genotype_2.wdl" as g2
import "CombineReassess.wdl" as creassess
import "Genotype_3.wdl" as g3

workflow Module04b{
    input{
        Array[File] regeno_beds
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
        File cluster_combined
        File cohort_samplelist
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_merge_list
    }
    call MergeList{
        input:
            regeno_beds = regeno_beds,
            prefix="master_regeno",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_merge_list
    }
    scatter (i in range(length(batches))){
        call g2.Regenotype as Genotype_2{
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

    call creassess.CombineReassess as CombineReassess{
        input:
            samplelist=cohort_samplelist,
            regeno_file=MergeList.master_regeno,
            cohort_cluster=cluster_combined,
            vcfs=Genotype_2.genotyped_vcf,
            sv_pipeline_docker=sv_pipeline_docker,
            sv_pipeline_base_docker=sv_pipeline_base_docker,
            runtime_attr_vcf2bed = runtime_attr_vcf2bed
    }
    
    scatter (i in range(length(Genotype_2.genotyped_vcf))){
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

task MergeList{
    input{
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

    command{
        cat ~{sep=' ' regeno_beds} |sort -k1,1V -k2,2n -k3,3n > ~{prefix}.bed
    }
    output{
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
    input{
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
