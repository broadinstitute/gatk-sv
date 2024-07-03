version 1.0

import "Structs.wdl"

workflow SVCodingConstraint {
    input{
        Int permutation_rounds
        File src_tar
        File gene_anno_tars
        File permutated_genes_tars
        File SV_sites_file
        File contig_file
        String sv_base_mini_docker

    }


    scatter(i in range(permutation_rounds)){

        call SVsVsGenesPart1{
            input:
                permu = i,
                src_tar = src_tar,
                gene_tars = permutated_genes_tars,
                SV_sites_file = SV_sites_file,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVsVsGenesPart2{
            input:
                permu = i,
                src_tar = src_tar,
                gene_anno_tars = gene_anno_tars,
                gene_tars = permutated_genes_tars,
                SV_sites_file = SV_sites_file,
                sv_base_mini_docker = sv_base_mini_docker,
                vs_3_prime_utr  = SVsVsGenesPart1.vs_3_prime_utr,
                vs_5_prime_utr  = SVsVsGenesPart1.vs_5_prime_utr,
                vs_promoter     = SVsVsGenesPart1.vs_promoter,
                vs_intact_exon_overlap  = SVsVsGenesPart1.vs_intact_exon_overlap,
                vs_partial_exon_overlap  = SVsVsGenesPart1.vs_partial_exon_overlap,
                vs_tss_transcripts_overlap  = SVsVsGenesPart1.vs_tss_transcripts_overlap,
                vs_partial_transcripts_overlap  = SVsVsGenesPart1.vs_partial_transcripts_overlap,
                vs_inside_exons = SVsVsGenesPart1.vs_inside_exons,
                vs_inside_introns = SVsVsGenesPart1.vs_inside_introns,
                vs_whole_transcript_overlap  = SVsVsGenesPart1.vs_whole_transcript_overlap
        }
    }

    output{
    	Array[File] gene_SV_rdata_list = SVsVsGenesPart2.gene_SV_rdata
    }
}



task SVsVsGenesPart1{
    input{
        Int permu
        File gene_tars
        File src_tar
        File SV_sites_file
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 15, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vs_whole_transcript_overlap = "~{filebase}.whole_transcript_overlap"
        File vs_3_prime_utr = "~{filebase}.3_prime_utr"
        File vs_5_prime_utr = "~{filebase}.5_prime_utr"
        File vs_intact_exon_overlap = "~{filebase}.intact_exon_overlap"
        File vs_partial_exon_overlap = "~{filebase}.partial_exon_overlap"
        File vs_tss_transcripts_overlap = "~{filebase}.tss_transcripts_overlap"
        File vs_partial_transcripts_overlap = "~{filebase}.partial_transcripts_overlap"
        File vs_inside_exons = "~{filebase}.inside_exons"
        File vs_inside_introns = "~{filebase}.inside_introns"
        File vs_promoter = "~{filebase}.promoter"
    }

    String filebase = basename(SV_sites_file,".gz")

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{gene_tars} ./
        tar zxvf gene_permu.tar.gz 

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 


        bedtools intersect -wo -a <(zcat ~{SV_sites_file} | cut -f1-5) -b gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz | bgzip >  ~{filebase}.transcript.bed.gz
        bedtools intersect -wo -a <(zcat ~{SV_sites_file} | cut -f1-5) -b gene_permu/r3.gencode.v39.ensembl.105.CDS.permu_~{permu}.bed.gz        |  bgzip > ~{filebase}.CDS.bed.gz
        bedtools intersect -wo -a <(zcat ~{SV_sites_file} | cut -f1-5) -b gene_permu/r3.gencode.v39.ensembl.105.intron.permu_~{permu}.bed.gz     | bgzip >  ~{filebase}.intron.bed.gz
        bedtools intersect -wo -a <(zcat ~{SV_sites_file} | cut -f1-5) -b gene_permu/r3.gencode.v39.ensembl.105.utr_3.permu_~{permu}.bed.gz      | bgzip >  ~{filebase}.utr_3.bed.gz
        bedtools intersect -wo -a <(zcat ~{SV_sites_file} | cut -f1-5) -b gene_permu/r3.gencode.v39.ensembl.105.utr_5.permu_~{permu}.bed.gz      | bgzip >  ~{filebase}.utr_5.bed.gz
        bedtools intersect -wo -a <(zcat ~{SV_sites_file} | cut -f1-5) -b gene_permu/r3.gencode.v39.ensembl.105.promoter.permu_~{permu}.bed.gz | bgzip >    ~{filebase}.promoter.bed.gz


        zcat ~{filebase}.transcript.bed.gz | awk '{if ($2<$7 && $3>$8) print}' | cut -f4,5,10,11  > ~{filebase}.whole_transcript_overlap 

        zcat ~{filebase}.transcript.bed.gz | awk '{if ($9=="+" && $2<$7+1 && $3<$8+1) print}' | cut -f4,5,10,11 > ~{filebase}.tss_transcripts_overlap
        zcat ~{filebase}.transcript.bed.gz | awk '{if ($9=="-" && $2>$7-1 && $3>$8-1) print}' | cut -f4,5,10,11 >> ~{filebase}.tss_transcripts_overlap
        zcat ~{filebase}.transcript.bed.gz | awk '{if ($9=="-" && $2<$7+1 && $3<$8+1) print}' | cut -f4,5,10,11 > ~{filebase}.partial_transcripts_overlap
        zcat ~{filebase}.transcript.bed.gz | awk '{if ($9=="+" && $2>$7-1 && $3>$8-1) print}' | cut -f4,5,10,11 >> ~{filebase}.partial_transcripts_overlap

        zcat ~{filebase}.utr_5.bed.gz | awk '{if ($9=="+" && $3<$8+1 && $3>$7-1) print}' | cut -f4,5,10,11 > ~{filebase}.5_prime_utr
        zcat ~{filebase}.utr_5.bed.gz | awk '{if ($9=="-" && $2<$8+1 && $2>$7-1) print}' | cut -f4,5,10,11 >> ~{filebase}.5_prime_utr
        zcat ~{filebase}.utr_3.bed.gz | awk '{if ($9=="+" && $2<$8+1 && $2>$7-1) print}' | cut -f4,5,10,11 > ~{filebase}.3_prime_utr
        zcat ~{filebase}.utr_3.bed.gz | awk '{if ($9=="-" && $3<$8+1 && $3>$7-1) print}' | cut -f4,5,10,11 >> ~{filebase}.3_prime_utr

        zcat ~{filebase}.CDS.bed.gz | awk '{if ($2>$7-1 && $3<$8+1) print}' | cut -f4,5,10,11  > ~{filebase}.inside_exons 
        zcat ~{filebase}.intron.bed.gz | awk '{if ($2>$7-1 && $3<$8+1) print}' | cut -f4,5,10,11  > ~{filebase}.inside_introns
        zcat ~{filebase}.promoter.bed.gz | awk '{if ($9=="+" && $3>$7-1 && $3<$8+1) print}' | cut -f4,5,10,11  > ~{filebase}.promoter 
        zcat ~{filebase}.promoter.bed.gz | awk '{if ($9=="-" && $2>$7-1 && $2<$8+1) print}' | cut -f4,5,10,11  >> ~{filebase}.promoter

        zcat ~{filebase}.transcript.bed.gz | awk '{if ($2>$7 && $3<$8) print}' | cut -f4,5,10,11  > ~{filebase}.SVs_inside_transcripts
        cut -f2  ~{filebase}.SVs_inside_transcripts | sort | uniq > ~{filebase}.SVs_inside_transcripts.gene_id

        Rscript ./src/categorize_intact_vs_partial_exon_overlap.R \
                -c ~{filebase}.CDS.bed.gz \
                -g ~{filebase}.SVs_inside_transcripts \
                -p ~{filebase}

    >>>

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

task SVsVsGenesPart2{
    input{
        File vs_whole_transcript_overlap 
        File vs_3_prime_utr  
        File vs_5_prime_utr 
        File vs_intact_exon_overlap 
        File vs_partial_exon_overlap 
        File vs_tss_transcripts_overlap 
        File vs_partial_transcripts_overlap 
        File vs_inside_exons
        File vs_inside_introns
        File vs_promoter

        Int permu
        File gene_tars
        File gene_anno_tars
        File src_tar
        File SV_sites_file
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 15, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File gene_SV_rdata = "~{filebase}.gene_SV_data.rData"
    }

    String vs_promoter.base = basename(vs_promoter,"")
    String vs_3_prime_utr.base = basename(vs_3_prime_utr,"")
    String vs_5_prime_utr.base = basename(vs_5_prime_utr,"")
    String vs_inside_exons.base = basename(vs_inside_exons,"")
    String vs_inside_introns.base = basename(vs_inside_introns,"")
    String vs_intact_exon_overlap.base = basename(vs_intact_exon_overlap,"")
    String vs_partial_exon_overlap.base = basename(vs_partial_exon_overlap,"")
    String vs_tss_transcripts_overlap.base = basename(vs_tss_transcripts_overlap,"")
    String vs_partial_transcripts_overlap.base = basename(vs_partial_transcripts_overlap,"")
    String vs_whole_transcript_overlap.base = basename(vs_whole_transcript_overlap,"")

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{gene_tars} ./
        tar zxvf gene_permu.tar.gz 

        gsutil cp ~{gene_anno_tars} ./
        tar zxvf gene_annotation.tar.gz

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 

        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_whole_transcript_overlap} -o ~{vs_whole_transcript_overlap.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_3_prime_utr}  -o ~{vs_3_prime_utr.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_5_prime_utr}  -o ~{vs_5_prime_utr.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_intact_exon_overlap}   -o ~{vs_intact_exon_overlap.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_partial_exon_overlap}  -o ~{vs_partial_exon_overlap.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_tss_transcripts_overlap}  -o ~{vs_tss_transcripts_overlap.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_partial_transcripts_overlap}  -o ~{vs_partial_transcripts_overlap.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_inside_exons} -o ~{vs_inside_exons.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_inside_introns} -o ~{vs_inside_introns.base}.reorganized
        Rscript ./src/reorganize_SVID_vs_gene.R -g gene_permu/r3.gencode.v39.ensembl.105.transcript.permu_~{permu}.bed.gz -i ~{vs_promoter} -o ~{vs_promoter.base}.reorganized

        Rscript ./src/integrate_SVID_vs_genes.across_different_overlaps.R -p ~{filebase}
        bgzip ~{filebase}.integrated

        Rscript src/calcu.gene.data.reaano.R \
            --sv_file_real ~{SV_sites_file} \
            --sv_vs_gene ~{filebase}.integrated.gz \
            --output ~{filebase}.gene_SV_data.rData \
            --gene_feature gene_annotation/gene_features.bed.gz \
            --gene_anno gene_annotation/genes_grch38_annotated_4_mapped_gencode_v39.CDS.UTR.tsv.gz \
            --gene_anno_phaplo_ptriplo gene_annotation/gene_information_loeuf_pHaplo_pTriplo_pLI_4.5.txt.gz \
            --gene_loeuf gene_annotation/loeuf.csv 
    >>>

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




