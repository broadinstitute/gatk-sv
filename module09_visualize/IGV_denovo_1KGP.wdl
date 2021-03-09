workflow IGV_denovo {
    File varfile
    File Fasta
    File Fasta_idx
    File Fasta_dict
    File Cram_file
    File Cram_file_idx
    String sample
    String chromosome

    call localize_cram{
      input:
        Cram_file = Cram_file,
        Cram_file_idx = Cram_file_idx,
        chromosome = chromosome,
        sample = sample,
        fasta = Fasta,
        fasta_idx = Fasta_idx,
        fasta_dict = Fasta_dict
    }

    call runIGV{
        input:
            varfile=varfile,
            fasta=Fasta,
            fasta_idx=Fasta_idx, 
            chromosome = chromosome,
            sample = sample,
            local_cram = localize_cram.local_cram,
            local_crai = localize_cram.local_crai
        }
    output{
        runIGV.*
    }
}

task localize_cram{
  File Cram_file
  File Cram_file_idx
  File fasta
  File fasta_idx
  File fasta_dict
  String chromosome
  String sample

  command <<<
    /gatk/gatk PrintReads \
       -R ${fasta} \
       -I ${Cram_file} \
       -L ${chromosome} \
       -O ${sample}_${chromosome}.bam
    samtools index ${sample}_${chromosome}.bam
  >>>
  runtime{
    docker: "broadinstitute/gatk:4.1.3.0"
    preemptible: 3
    memory: "10 GB"
    disks: "local-disk 50 HDD"
  }
  output{
    File local_cram = '${sample}_${chromosome}.bam'
    File local_crai = '${sample}_${chromosome}.bam.bai'
  }
}

task runIGV{
    File varfile
    File fasta
    File fasta_idx
    String sample 
    String chromosome
    File local_cram
    File local_crai
    command <<<
    	set -e
        #export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
        python /makeigvpesr_cram.py ${varfile} 500 ${fasta} ${local_cram} ${sample} ${chromosome}
        #cat pe.sh
        bash pe.sh
        xvfb-run --server-args="-screen 0, 1920x1080x24" bash /IGV_2.4.14/igv.sh -b pe.txt
        tar -czf pe_screenshots.tar.gz pe_screenshot
        #tar -czf pe_bam.tar.gz pe_bam
        python /makeigvsplit_cram.py ${varfile} 500 ${fasta} ${local_cram} ${sample} ${chromosome}
        #cat sr.sh
        bash sr.sh
        xvfb-run --server-args="-screen 0, 1920x3000x24" bash /IGV_2.4.14/igv.sh -b sr.txt
        tar -czf sr_screenshots.tar.gz sr_screenshot
        #tar -czf pe_bam.tar.gz pe_bam
    >>>
  runtime {
      docker: "talkowski/igv_gatk:latest"
      preemptible: 3
      memory: "10 GB"
      disks: "local-disk 50 HDD"
  }
  output{
      File pe_plots="pe_screenshots.tar.gz"
      File sr_plots="sr_screenshots.tar.gz"
      #File bams="pe_bam.tar.gz"
      #File igvfile="pe.txt"
      #File igvbam="pe.sh"
  }
}
