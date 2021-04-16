version 1.0

import "IGVTrioPlots.wdl" as igv
import "Structs.wdl"

workflow IGV_all_samples {
    input {
        Array[String] pb_list
        Array[String] fa_list
        Array[String] mo_list
        Array[File] pb_cram_list
        Array[File] pb_crai_list
        Array[File] fa_cram_list
        Array[File] fa_crai_list
        Array[File] mo_cram_list
        Array[File] mo_crai_list
        File varfile
        File Fasta
        File Fasta_dict
        File Fasta_idx
        String prefix
        String sv_base_mini_docker
        String igv_docker
        RuntimeAttr? runtime_attr_override
    }

    scatter (i in range(length(pb_list))){
        call generate_per_sample_bed{
            input:
                varfile = varfile,
                sample_id = pb_list[i],
                sv_base_mini_docker=sv_base_mini_docker,
                runtime_attr_override=runtime_attr_override
        }

        call igv.IGV_trio as IGV_trio {
            input:
                varfile=generate_per_sample_bed.per_sample_varfile,
                Fasta = Fasta,
                Fasta_idx = Fasta_idx,
                Fasta_dict = Fasta_dict,
                pb=pb_list[i],
                fa=fa_list[i],
                mo=mo_list[i],
                pb_cram=pb_cram_list[i],
                fa_cram=fa_cram_list[i],
                mo_cram=mo_cram_list[i],
                pb_crai=pb_crai_list[i],
                fa_crai=fa_crai_list[i],
                mo_crai=mo_crai_list[i],
                igv_docker = igv_docker
                }
        }
    call integrate_igv_plots{
        input:
            igv_tar = IGV_trio.tar_gz_pe,
            prefix = prefix, 
            sv_base_mini_docker = sv_base_mini_docker
    }

    output{
        File tar_gz_pe = integrate_igv_plots.plot_tar
    }
    }


task generate_per_sample_bed{
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
        grep -w ~{sample_id} ~{varfile} | cut -f1-5 | awk '{print $1,$2,$3,$5,$4}' | sed -e 's/ /\t/g' > ~{filename}.~{sample_id}.bed
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

task integrate_igv_plots{
    input {
        Array[File] igv_tar
        String prefix
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

    command <<<
        mkdir ~{prefix}_igv_plots
        while read file; do
            tar -zxf ${file}
            mv pe_igv_plots/*  ~{prefix}_igv_plots/
        done < ~{write_lines(igv_tar)};
        tar -czf ~{prefix}_igv_plots.tar.gz ~{prefix}_igv_plots
    >>>

    output{
        File plot_tar = "~{prefix}_igv_plots.tar.gz"
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
