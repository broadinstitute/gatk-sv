version 1.0

import "Structs.wdl"

workflow IGV_trio {
    input{
        File varfile
        File Fasta
        File Fasta_idx
        File Fasta_dict
        String pb
        String fa
        String mo
        File pb_cram
        File pb_crai
        File fa_cram
        File fa_crai
        File mo_cram
        File mo_crai
        String igv_docker
    }

    call runIGV_whole_genome{
        input:
            varfile = varfile,
            fasta = Fasta,
            fasta_dict = Fasta_dict,
            fasta_idx = Fasta_idx,
            fa = fa,
            mo = mo,
            pb = pb,
            pb_cram =  pb_cram,
            pb_crai =  pb_crai,
            fa_cram =  fa_cram,
            fa_crai =  fa_crai,
            mo_cram =  mo_cram,
            mo_crai =  mo_crai,
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
        File fasta_dict
        String pb
        String fa
        String mo
        File pb_cram
        File pb_crai
        File fa_cram
        File fa_crai
        File mo_cram
        File mo_crai
        String igv_docker
    }
    command <<<
            set -euo pipefail
            #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
            python /src/makeigvpesr_trio.py ~{varfile} ~{fasta} ~{pb} ~{pb_cram},~{fa_cram},~{mo_cram} pe_igv_plots -b 500
            bash pe.sh
            xvfb-run --server-args="-screen 0, 1920x540x24" bash /IGV_2.4.14/igv.sh -b pe.txt
            tar -czf ~{pb}_pe_igv_plots.tar.gz pe_igv_plots

        >>>
    runtime {
        docker: igv_docker
        preemptible: 1
        memory: "15 GB"
        disks: "local-disk 100 HDD"
        }
    output{
        File pe_plots="~{pb}_pe_igv_plots.tar.gz"
        File pe_txt = "pe.txt"
        }
    }

