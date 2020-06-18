version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "igv_generate_plots.wdl" as igv
import "Structs.wdl"

workflow IGV_all_samples {
    input {
        Array[String] samples
        Array[String] crams
        Array[String] crams_idx
        File varfile
        File Fasta
        File Fasta_dict
        File Fasta_idx
        String prefix
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    scatter (i in range(length(samples))){
        call GeneratePerSampleBed{
            input:
                varfile = varfile,
                sample_id = samples[i],
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_override
        }

        call igv.IGVdenovo as IGVdenovo {
            input:
                varfile=GeneratePerSampleBed.per_sample_varfile,
                sample = samples[i],
                Cram_file = crams[i],
                Cram_file_idx = crams_idx[i],
                Fasta = Fasta,
                Fasta_idx = Fasta_idx,
                Fasta_dict = Fasta_dict, 
                prefix = prefix
        }

    }
    call IntegrateFigure{
        input:
            pe_tar_gz = IGVdenovo.tar_gz_pe,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker
   }

    output{
        File tar_gz_pe = IntegrateFigure.tar_gz_pe
    }
    }

task GeneratePerSampleBed{
    input {
        File varfile
        String sample_id
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }
    RuntimeAttr default_attr=object {
        cpu_cores: 1,
        mem_gb: 1,
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    String filename = basename(varfile, ".bed")
    command <<<
        set -euo pipefail
        grep -w ~{sample_id} ~{varfile} | cut -f1-5 > ~{filename}.~{sample_id}.bed
        >>>

    output{
        File per_sample_varfile= "~{filename}.~{sample_id}.bed"
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
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

task IntegrateFigure{
    input{
        Array[File] pe_tar_gz
        String prefix
        String sv_base_mini_docker
    }
    command <<<
        mkdir ~{prefix}_igv_screenshot/
        while read file; do
            tar -zxvf ${file}
            mv pe_screenshot/* ~{prefix}_igv_screenshot/
        done < ~{write_lines(pe_tar_gz)};

        tar -czf ~{prefix}_igv_screenshot.tar.gz ~{prefix}_igv_screenshot
    >>>

    runtime{
        docker: sv_base_mini_docker
        preemptible: 3
        memory: "10 GB"
        disks: "local-disk 50 SSD"    
    }
    output{
        File tar_gz_pe = "~{prefix}_igv_screenshot.tar.gz"
    }
}
