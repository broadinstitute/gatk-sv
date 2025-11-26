version 1.0

import "Structs.wdl"
import "CalculateFstGnomad.wdl" as calculate_fst

workflow CalculateFstLargeVcf {
    input{
        File vcf
        File vcf_idx
        File samp_pop
        String variant_type 
        String sv_fst_docker
        String sv_base_mini_docker
        File? region_bed
        RuntimeAttr? runtime_tabix_vcf
        RuntimeAttr? runtime_attr_fst
    }

    if (defined(region_bed)){

        Array[Array[String]] regions = read_tsv(select_first([region_bed]))

        scatter(region in regions){
            call TabixVcf{
                input:
                    vcf = vcf,
                    vcf_idx = vcf_idx, 
                    region = region[0],
                    sv_base_mini_docker = sv_base_mini_docker,
                    runtime_attr_override = runtime_tabix_vcf
            }

            call calculate_fst.CalculateFstGnomad as Calculate_fst_per_shard{
                input:
                    vcf = TabixVcf.split_vcf,
                    vcf_idx = TabixVcf.split_vcf_idx,
                    samp_pop = samp_pop,
                    variant_type = variant_type,
                    sv_fst_docker = sv_fst_docker,
                    runtime_attr_fst = runtime_attr_fst
            }
        }

        String prefix = basename(vcf,".vcf.gz")

        call ConcatBeds{
            input:
                shard_bed_files = Calculate_fst_per_shard.output_fst_sites,
                prefix = prefix,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call CalcuFstPop{
            input:
                Fst_sites = ConcatBeds.merged_file,
                sv_fst_docker = sv_fst_docker,
                runtime_attr_override = runtime_attr_fst_pop_from_sites
        }
    }

    if (!defined(region_bed)){
        call calculate_fst.CalculateFstGnomad as Calculate_fst{
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                samp_pop = samp_pop, 
                variant_type = variant_type,
                sv_fst_docker = sv_fst_docker,
                runtime_attr_fst = runtime_attr_fst
        }
    }


    File Fst_sites = select_first([ConcatBeds.merged_file, Calculate_fst.output_fst_sites])
    File Fst_pop = select_first([CalcuFstPop.fst_pop, Calculate_fst.output_fst_pop ])

    output{
        File output_Fst_sites = Fst_sites
        File output_Fst_pop = Fst_pop
    }
}


task CalcuFstPop{
    input{
        File Fst_sites
        String sv_fst_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File fst_pop = "~{filebase}.pop"
    }

    String filebase = basename(Fst_sites,".sites")

    command <<<
        set -Eeuo pipefail

        python /src/Calcu_Fst_pop_from_sites.py -i ~{Fst_sites} -p ~{filebase}.pop
   >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_fst_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }    
}

task TabixVcf{
    input{
        File vcf
        File vcf_idx
        String region
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 5, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File split_vcf = "~{filebase}.tmp.vcf.gz"
        File split_vcf_idx = "~{filebase}.tmp.vcf.gz.tbi"
    }

    String filebase = basename(vcf,".vcf.gz")

    command <<<
        set -Eeuo pipefail

        bcftools view ~{vcf} ~{region} | bgzip > ~{filebase}.tmp.vcf.gz
        tabix -p vcf ~{filebase}.tmp.vcf.gz
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

task ConcatBeds {
  input {
    Array[File] shard_bed_files
    String prefix
    Boolean? index_output
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Boolean call_tabix = select_first([index_output, true])
  String output_file="~{prefix}.sites"

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(shard_bed_files, "GB")
  RuntimeAttr runtime_default = object {
    mem_gb: 2.0,
    disk_gb: ceil(10.0 + input_size * 7.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eux

    # note head -n1 stops reading early and sends SIGPIPE to zcat,
    # so setting pipefail here would result in early termination
    head -n1 ~{shard_bed_files[0]} > header.txt

    # no more early stopping
    set -o pipefail

    while read SPLIT; do
      tail -n+2 $SPLIT 
    done < ~{write_lines(shard_bed_files)} \
      | cat header.txt - \
      > ~{output_file}

  >>>

  output {
    File merged_file = output_file
  }
}
