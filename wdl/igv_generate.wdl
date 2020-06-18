workflow IGV_denovo {
    File varfile
    File Fasta
    File Fasta_idx
    call runIGV{
        input:
            varfile=varfile,
            fasta=Fasta,
            fasta_idx=Fasta_idx
        }
    output{
        runIGV.*
    }
}
task runIGV{
    File varfile
    File fasta
    File fasta_idx
    command <<<
        set -euo pipefail
        export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        python /makeigvpesr_cram.py ${varfile} 500 ${fasta}
        cat pe.sh
        bash pe.sh
        xvfb-run --server-args="-screen 0, 1920x1080x24" bash /IGV_2.4.14/igv.sh -b pe.txt
        tar -czf pe_screenshots.tar.gz pe_screenshot
        #tar -czf pe_bam.tar.gz pe_bam
        python /makeigvsplit_cram.py ${varfile} 500 ${fasta}
        cat sr.sh
        bash sr.sh
        xvfb-run --server-args="-screen 0, 1920x3000x24" bash /IGV_2.4.14/igv.sh -b sr.txt
        tar -czf sr_screenshots.tar.gz sr_screenshot
        #tar -czf pe_bam.tar.gz pe_bam
    >>>
  runtime {
      docker: "talkowski/igv:latest"
      preemptible: 3
      memory: "10 GB"
      disks: "local-disk 50 SSD"
  }
  output{
      File pe_plots="pe_screenshots.tar.gz"
      File sr_plots="sr_screenshots.tar.gz"
      #File bams="pe_bam.tar.gz"
      #File igvfile="pe.txt"
      #File igvbam="pe.sh"
  }
}
