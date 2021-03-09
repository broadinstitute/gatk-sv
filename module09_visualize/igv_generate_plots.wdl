version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs.wdl"

workflow IGV_denovo {
    input{
        File varfile
        File Fasta
        File Fasta_idx
        File Fasta_dict
        File Cram_file
        File Cram_file_idx
        Array[String] contigs
        String sample
        String prefix
        String igv_docker
        String sv_base_mini_docker
    }
    scatter (contig in contigs){
        call runIGV{
            input:
                varfile=varfile,
                fasta=Fasta,
                fasta_idx=Fasta_idx,
                sample =sample,
                prefix=prefix,
                local_cram=Cram_file,
                local_crai=Cram_file_idx,
                contig = contig,
                igv_docker = igv_docker
        }
       
    }
    call integrate_figure{
        input:
            pe_tar_gz = runIGV.pe_plots,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker
   }

     output{
        File tar_gz_pe = integrate_figure.tar_gz_pe
   }
}

task runIGV{
    input{
        File varfile
        File fasta
        File fasta_idx
        String sample 
        String prefix
        String contig
        File local_cram
        File local_crai
        String igv_docker
    }
    command <<<
              set -euo pipefail
              awk '{print $1,$2,$3,$5,$4}' ~{varfile} | sed -e 's/ /\t/g' | awk '{if ($1=="~{contig}") print}' > var_file
              python /src/makeigvsplit_cram.py var_file 500 ~{fasta} ~{local_cram} ~{sample} all
              bash pe.sh
              xvfb-run --server-args="-screen 0, 1920x3000x24" bash /IGV_2.4.14/igv.sh -b pe.txt

              tar -czf ~{prefix}.pe_screenshots.tar.gz pe_screenshot
        >>>
    runtime {
        docker: igv_docker
        preemptible: 3
        memory: "10 GB"
        disks: "local-disk 50 HDD"
        }
    output{
        File pe_plots="~{prefix}.pe_screenshots.tar.gz"
        File igv="pe.txt"
        }
    }

task integrate_figure{
    input{
        Array[File] pe_tar_gz
        String prefix
        String sv_base_mini_docker
    }
    command <<<
        mkdir ~{prefix}_igv_plots/
        while read file; do
            tar -zxvf ${file}
            mv pe_screenshot/* ~{prefix}_igv_plots/
        done < ~{write_lines(pe_tar_gz)};

        tar -czf ~{prefix}_igv_plots.tar.gz ~{prefix}_igv_plots
    >>>

    runtime{
        docker: sv_base_mini_docker
        preemptible: 3
        memory: "10 GB"
        disks: "local-disk 50 HDD"    
    }
    output{
        File tar_gz_pe = "~{prefix}_igv_plots.tar.gz"
    }
}
