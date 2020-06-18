version 1.0

##########################################################################################

## Github commit: talkowski-lab/gatk-sv-v1:<ENTER HASH HERE IN FIRECLOUD>

##########################################################################################

import "Structs.wdl"

workflow IGVTrio {
    input{
        File varfile
        File Fasta
        File Fasta_idx
        File Fasta_dict
        File ped_file
        File sample_cram
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

    call RunIGVWholeGeome{
        input:
            varfile = varfile,
            fasta = Fasta,
            fasta_dict = Fasta_dict,
            fasta_idx = Fasta_idx,
            ped_file = ped_file,
            cram_list = sample_cram,
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
        File tar_gz_pe = RunIGVWholeGeome.pe_plots
    }
}

task RunIGVWholeGeome{
    input{
        File varfile
        File fasta
        File fasta_idx
        File fasta_dict
        File ped_file
        File cram_list
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

