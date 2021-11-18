version 1.0

# Author: Xuefang Zhao <xzhao12@mgh.harvard.edu>
import "Structs.wdl"

workflow ArrayVali{
    input{
        File input_vcf
        File input_vcf_index
        File array_metrics
        File samplelist   
        File contiglist
        String ref_version

        String sv_pipeline_qc_docker
        String sv_base_mini_docker
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_split_vcf
        RuntimeAttr? runtime_attr_split_samples_from_vcf
        RuntimeAttr? runtime_attr_revise_vcf
        RuntimeAttr? runtime_override_array_validation
    }

    Array[Array[String]] contigs=read_tsv(contiglist)
    scatter (contig in contigs) {
        call SplitVcf{
            input:
                vcf = input_vcf,
                vcf_index = input_vcf_index,
                contig = contig[0],
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_split_vcf
        }

        call SplitSamplesFromVcf{
            input:
                vcf = SplitVcf.out_vcf,
                vcf_index = SplitVcf.out_idx,
                samplelist = samplelist,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_split_samples_from_vcf
        }

        call ReviseVcfFile{
            input:
                vcf = SplitSamplesFromVcf.out_vcf,
                vcf_index = SplitSamplesFromVcf.out_idx,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_revise_vcf
        }

        call runIRS{
            input:
                vcf = ReviseVcfFile.out_vcf,
                vcf_index = ReviseVcfFile.out_idx,
                array_metrics = array_metrics,
                samplelist = samplelist,
                ref_version = ref_version,
                sv_pipeline_qc_docker = sv_pipeline_qc_docker,
                runtime_attr_override = runtime_override_array_validation
        }
    }
    output{
        Array[File] vali_vcf = runIRS.out_vcf
        Array[File] vali_vcf_idx = runIRS.out_idx
        Array[File] vali_metrics = runIRS.out_metrics
    }
}

task SplitSamplesFromVcf{
    input{
        File vcf
        File vcf_index
        File samplelist   
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GiB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 1.0
    RuntimeAttr runtime_default = object {
      mem_gb: 5.0,
      disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
      cpu_cores: 1,
      preemptible_tries: 3,
      max_retries: 1,
      boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
      memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
      disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
      cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
      preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
      maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
      docker: sv_pipeline_docker
      bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String filename = basename(vcf, ".vcf.gz")

    command <<<
        set -euo pipefail
        bcftools view -S ~{samplelist   } ~{vcf} | bgzip > ~{filename}.target_samples.vcf.gz
        tabix -p vcf  ~{filename}.target_samples.vcf.gz        
    >>>

    output{
        File out_vcf = "~{filename}.target_samples.vcf.gz"
        File out_idx = "~{filename}.target_samples.vcf.gz.tbi"
    }
}

task ReviseVcfFile{
    input{
        File vcf
        File vcf_index
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size = size(vcf, "GiB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 1.0
    RuntimeAttr runtime_default = object {
        mem_gb: 5.0,
        disk_gb: ceil(base_disk_gb + input_size * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String filename = basename(vcf, ".vcf.gz")

    command <<<
        svtk vcf2bed ~{vcf} ~{filename}.bed
        zcat ~{vcf} | grep "#" > ~{filename}.revised.vcf
         paste \
             <(paste -d ';' \
             <(zcat ~{vcf} | grep -v "#" | cut -f1-8) \
             <(cut -f6 ~{filename}.bed | sed -e "s/^/SAMPLES=/")) \
             <(zcat ~{vcf} | grep -v "#" | cut -f9-) \
             >> ~{filename}.revised.vcf
         bgzip ~{filename}.revised.vcf
        tabix -p vcf ~{filename}.revised.vcf.gz
    >>>

    output{
        File out_vcf = "~{filename}.revised.vcf.gz"
        File out_idx = "~{filename}.revised.vcf.gz.tbi"
    }
}

task runIRS{
    input{
        File vcf
        File vcf_index
        File array_metrics
        File samplelist   
        String ref_version

        String sv_pipeline_qc_docker
        RuntimeAttr? runtime_attr_override
    }

    Float input_size_1 = size(vcf, "GiB")
    Float input_size_2 = size(array_metrics, "GiB")
    Float base_disk_gb = 10.0
    Float input_disk_scale = 1.0
    RuntimeAttr runtime_default = object {
        mem_gb: 50.0,
        disk_gb: ceil(base_disk_gb + input_size_1 * input_disk_scale + input_size_2 * input_disk_scale),
        cpu_cores: 1,
        preemptible_tries: 3,
        max_retries: 1,
        boot_disk_gb: 10
    }
    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_qc_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String filename = basename(vcf, ".vcf.gz")

    command <<<
        bash runIRS.sh \
            -s ~{vcf} \
            -o ~{filename}.output.vcf \
            -r ~{filename}.report.dat  \
            -d ~{samplelist   } \
            -g 38 \
            -a ~{array_metrics}"
        bgzip ~{filename}.output.vcf
        tabix -p vcf ~{filename}.output.vcf.gz
    >>>

    output{
        File out_vcf = "~{filename}.output.vcf.gz"
        File out_idx = "~{filename}.output.vcf.gz.tbi"
        File out_metrics = "~{filename}.report.dat"
    }
}

task SplitVcf {
    input {
        File vcf
        File vcf_index
        String contig
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
        File out_vcf = "~{contig}.vcf.gz"
        File out_idx = "~{contig}.vcf.gz.tbi"
    }
    command <<<

        set -euxo pipefail
        zcat ~{vcf} | sed -n -e '/^#/p' > ~{contig}.vcf
        zcat ~{vcf} | sed -e '/^#/d' | awk '{if ($1=="~{contig}") print}' >> ~{contig}.vcf
        bgzip ~{contig}.vcf
        tabix -p vcf ~{contig}.vcf.gz

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



