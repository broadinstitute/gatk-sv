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
        String sample
        String prefix
        String igv_docker
    }
    call runIGV{
        input:
            varfile=varfile,
            fasta=Fasta,
            fasta_idx=Fasta_idx,
            sample =sample,
            prefix=prefix,
            local_cram=Cram_file,
            local_crai=Cram_file_idx,
            igv_docker = igv_docker
    }
    output{
        File tar_gz_pe = runIGV.pe_plots
        File igv_rec = runIGV.igv
    }
}

task runIGV{
    input{
        File varfile
        File fasta
        File fasta_idx
        String sample 
        String prefix
        String igv_docker
        File local_cram
        File local_crai
    }
    command <<<
              set -euo pipefail
              awk '{print $1,$2,$3,$5,$4}' ~{varfile} | sed -e 's/ /\t/g' > var_file
              python /makeigvpesr_cram.py var_file 500 ~{fasta} ~{local_cram} ~{sample}
              bash pe.sh
              xvfb-run --server-args="-screen 0, 1920x3000x24" bash /IGV_2.4.14/igv.sh -b pe.txt

              tar -czf ~{prefix}.pe_screenshots.tar.gz pe_screenshot
            fi
        >>>
    runtime {
        docker: igv_docker
        preemptible: 3
        memory: "10 GB"
        disks: "local-disk 50 SSD"
        }
    output{
        File pe_plots="~{prefix}.pe_screenshots.tar.gz"
        File igv="pe.txt"
        }
    }



