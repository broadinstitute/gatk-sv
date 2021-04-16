version 1.0

import "Structs.wdl"

workflow IGV_denovo {
    input{
        File varfile
        File Fasta
        File Fasta_idx
        File Fasta_dict
        File Cram_file
        File Cram_file_idx
        String sample
        String prefix
        String igv_docker
        String sv_base_mini_docker
    }
    call runIGV_whole_genome{
            input:
                varfile=varfile,
                fasta=Fasta,
                fasta_idx=Fasta_idx,
                sample =sample,
                prefix=prefix,
                local_cram=Cram_file,
                local_crai=Cram_file_idx,
                var_file = varfile,
                igv_docker = igv_docker
        }

     output{
        File tar_gz_pe = runIGV_whole_genome.pe_plots
   }
}

task runIGV_whole_genome{
    input{
        File varfile
        File fasta
        File fasta_idx
        String sample 
        String prefix
        File local_cram
        File local_crai
        File var_file
        String igv_docker
    }
    command <<<
              set -euo pipefail
              python /src/makeigvsplit_cram.py ~{var_file} 500 ~{fasta} ~{local_cram} ~{sample} all
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

