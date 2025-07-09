version 1.0

import "MergeVcfsByChromosome.wdl" as MergeVcfsByChromosome

# Runs a truvari inter-sample merge
# adapted to do whole-genome merge using the data table but otherwise not changing the parameters at all
workflow TruvariIntersample {
    input {
        Array[File] input_vcfs
        Array[File] input_vcfs_idx
        Array[String] sample_ids
        String chrom

        String sv_base_mini_docker
        
        File reference_fai
        File density_counter_py
        Int max_records_per_chunk = 10000
    }

    parameter_meta {
        max_records_per_chunk: "Discards chunks that contain more than this many records. Setting it to 10k keeps 99.9% of all chunks in AoU Phase 1 (1027 samples) on CHM13."
    }

    scatter (idx in range(length(input_vcfs))) {

        call MergeVcfsByChromosome.ExtractChromosomeVcf {
            input:
              input_vcf = input_vcfs[idx],
              input_vcf_idx = input_vcfs_idx[idx],
              chromosome = chrom,
              sv_base_mini_docker = sv_base_mini_docker
          }
    }
    
    call TruvariIntersampleImpl {
            input:
                vcfs = ExtractChromosomeVcf.output_vcf,
                sample_ids = sample_ids,
                reference_fai = reference_fai,
                density_counter_py = density_counter_py,
                max_records_per_chunk = max_records_per_chunk
    }
    
    output {
        File merged_vcf = TruvariIntersampleImpl.merged_vcf
        File merged_vcf_tbi = TruvariIntersampleImpl.merged_vcf_tbi
        File bcftools_merge_vcf = TruvariIntersampleImpl.bcftools_merge_vcf
        File bcftools_merge_vcf_tbi = TruvariIntersampleImpl.bcftools_merge_vcf_tbi
    }
}


task TruvariIntersampleImpl {
    input {
        Array[File] vcfs
        Array[String] sample_ids
        File reference_fai
        File density_counter_py
        Int max_records_per_chunk
    }
    
    Int ram_size_gb = 64  # Arbitrary
    
    command <<<
        set -euxo pipefail
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        echo ~{sep=' ' sample_ids} | sed 's/ /\n/g' > sample_ids.txt
        echo ~{sep=' ' vcfs} | sed 's/ /\n/g' > vcfs.txt
        paste sample_ids.txt vcfs.txt > samples_vcfs.txt
        while read sample vcf; do
            if [[ ! "$vcf" == *"$sample"* ]]; then echo "ERROR: possible sample/vcf mismatch for sample $sample, vcf $vcf" >&2; fi
            echo $sample > tmp_sample_file.txt
            bcftools reheader -s tmp_sample_file.txt $vcf -o renamed_${sample}.vcf.gz
            tabix renamed_${sample}.vcf.gz
        done < samples_vcfs.txt 

        # BCFTOOLS
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none renamed_*vcf.gz --output-type z > merged.vcf.gz
        tabix -f merged.vcf.gz
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z merged.vcf.gz > normed.vcf.gz
        tabix -f normed.vcf.gz

        # Discarding chunks with too many records
        python ~{density_counter_py} normed.vcf.gz > chunks.bed
        awk '$4 >= ~{max_records_per_chunk}' chunks.bed > excluded.bed
        bedtools complement -i excluded.bed -g ~{reference_fai} > included.bed
        bedtools sort -faidx ~{reference_fai} -i included.bed > included.sorted.bed
        
        # TRUVARI
        truvari collapse --input normed.vcf.gz --collapsed-output removed.vcf.gz --sizemin 0 --sizemax 1000000 --keep common --bed included.sorted.bed --gt all \
            | bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z > collapsed.vcf.gz
        tabix -f collapsed.vcf.gz

        mv normed.vcf.gz bcftools_merge.vcf.gz
        tabix -f bcftools_merge.vcf.gz
    >>>
    
    output {
        File merged_vcf = "collapsed.vcf.gz"
        File merged_vcf_tbi = "collapsed.vcf.gz.tbi"
        File sample_vcf_correspondence = "samples_vcfs.txt" # just to check that the sample:vcf correspondence is correct
        File bcftools_merge_vcf = "bcftools_merge.vcf.gz"
        File bcftools_merge_vcf_tbi = "bcftools_merge.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 8
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}