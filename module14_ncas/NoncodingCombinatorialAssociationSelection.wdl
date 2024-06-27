version 1.0

import "Structs.wdl"
import "SVvsConservative.wdl" as SVvsConservative

workflow NoncodingCombinatorialAssociationSelection {
    input{
        Int permutation_rounds
        File aps
        File SV_sites_vcf
        File contig_file
    }

    scatter(i in range(length(permutation_rounds))){
        call GeneratePermutatedSVs{
            input:
                permu = i, 
                SV_sites_vcf = SV_sites_vcf,
                ncas_docker = ncas_docker
        }

        call SVvsConservative.SVvsConservative{
            input:
                SV_file = GeneratePermutatedSVs.permutated_SV,
                contig_file = contig_file,
                ncas_docker = ncas_docker                
        }

        call SVvsGencode{
            input: 
                SV_file = GeneratePermutatedSVs.permutated_SV,
                ncas_docker = ncas_docker
        }

        call SVvsNoncoding{
            input:
                SV_file = GeneratePermutatedSVs.permutated_SV,
                ncas_docker = ncas_docker
        }

        call SVvsGene{
            input:
                SV_file = GeneratePermutatedSVs.permutated_SV,
                ncas_docker = ncas_docker

        }

        call GenerateNcasMetrics{
            input:
                aps = aps,
                sv_file_real = SV_sites_vcf,
                sv_file_permu = GeneratePermutatedSVs.permutated_SV,
                sv_vs_gencode = SVvsGencode.SV_vs_gencode
                sv_vs_conserve    = SVvsConservative.SV_vs_conserved_elements
                sv_vs_noncoding = SVvsNoncoding.SV_vs_noncoding,
                sv_vs_gene = SVvsGene.SV_vs_trans,
                sv_vs_coding = SVvsGene.SV_vs_cds,
                ncas_docker = ncas_docker
        }
    }
    output{
        Array[File] ncas_data_list = GenerateNcasMetrics.ncas_rdata
    }
}

task GeneratePermutatedSVs{
    input{
        Int permu
        File SV_sites_vcf
        String ncas_docker
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
        File permutated_SV = "~{filebase}.permu_~{permu}.corrected.gz"
    }

    String filebase = basename(SV_sites_vcf,".gz")

    command <<<
        set -Eeuo pipefail


        Rscript src/generate_SV_permutations.R -p ~{permu} -i ~{SV_sites_vcf} -o ~{filebase}.permu_~{permu}
        bgzip ~{filebase}.permu_~{permu}

        Rscript src/correct_permutated_SVs.R -i ~{filebase}.permu_~{permu}.gz
        bgzip ~{filebase}.permu_~{permu}.corrected

    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task GenerateNcasMetrics{
    input{
        File aps
        File sv_file_real
        File sv_file_permu
        File sv_vs_gencode
        File sv_vs_conserve    
        File sv_vs_noncoding
        File sv_vs_gene
        File sv_vs_coding

        String ncas_docker
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
        File ncas_rdata = "~{filebase}.rData"
    }

    String filebase = basename(sv_file_permu,".gz")

    command <<<
        set -Eeuo pipefail

             Rscript src/generate_cwas_metrics.sh \
             --sv_file_real ~{sv_file_real} \
             --sv_file_permu ~{sv_file_permu} \
             --sv_vs_gencode ~{sv_vs_gencode} \
             --sv_vs_conserve ~{sv_vs_conserve} \
             --sv_vs_noncoding ~{sv_vs_noncoding} \
             --sv_vs_gene ~{sv_vs_gene} \
             --sv_vs_coding ~{sv_vs_coding} \
             --aps ~{aps} \
             --output ~{filebase}.rData
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SVvsGene{
    input{
        File SV_file
        String ncas_docker
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
        File SV_vs_cds = "~{filebase}.coding_interruptive.bed"
        File SV_vs_trans = "~{filebase}.genic.bed"
    }

    String filebase = basename(SV_file,".gz")

    command <<<
        set -Eeuo pipefail
        bedtools coverage -wo -a ~{SV_file} -b gene/r3.gencode.v39.ensembl.105.CDS.sorted.bed.gz | awk '{if ($NF>0) print}' > ~{filebase}.coding_interruptive.bed
        bedtools coverage -wo -a ~{SV_file} -b gene/r3.gencode.v39.ensembl.105.transcript.sorted.bed.gz | awk '{if ($NF>0) print}' > ~{filebase}.genic.bed
        bgzip ~{filebase}.coding_interruptive.bed
        bgzip ~{filebase}.genic.bed

    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SVvsGencode{
    input{
        File SV_file
        String ncas_docker
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
        File SV_vs_gencode = "~{filebase}.vs.gencode.integrated.gz"
    }

    String filebase = basename(SV_file,".gz")

    command <<<
        set -Eeuo pipefail

        zcat ~{SV_file} | cut -f1-4 > input_SVs.bed
        zcat gencode/gencode.v39.integrated.bed.gz | cut -f1-4,8 > gencode.bed
        bedtools intersect -wo -a input_SVs.bed -b gencode.bed > ~{filebase}.vs.gencode.bed
        cat gencode/header ~{filebase}.vs.gencode.bed | bgzip >  ~{filebase}.vs.gencode.bed.gz
        Rscript src/integrate_SV_vs_gencode.R -i ~{filebase}.vs.gencode.bed.gz  -o ~{filebase}.vs.gencode.integrated
        bgzip ~{filebase}.vs.gencode.integrated
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SVvsNoncoding{
    input{
        File SV_file
        String ncas_docker
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
        File SV_vs_noncoding = "~{filebase}.vs.nc_elements.integrated"
    }

    String filebase = basename(SV_file,".gz")

    command <<<
        set -Eeuo pipefail
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/vista_cores.neg.vs_zscore.uniq  | sed -e  's/$/\tvista_cores.neg/' | bgzip > ~{filebase}.vs.vista_cores.neg.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/vista_cores.pos.vs_zscore.uniq| sed -e  's/$/\tvista_cores.pos/' | bgzip > ~{filebase}.vs.vista_cores.pos.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.general.vs_zscore.uniq| sed -e  's/$/\tDCR2.general/' | bgzip > ~{filebase}.vs.DCR2.general.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.adrenal_glands.vs_zscore.uniq| sed -e  's/$/\tDCR2.adrenal_glands/' | bgzip > ~{filebase}.vs.DCR2.adrenal_glands.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.brain.vs_zscore.uniq| sed -e  's/$/\tDCR2.brain/' | bgzip > ~{filebase}.vs.DCR2.brain.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.eye.vs_zscore.uniq| sed -e  's/$/\tDCR2.eye/' | bgzip > ~{filebase}.vs.DCR2.eye.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.face.vs_zscore.uniq| sed -e  's/$/\tDCR2.face/' | bgzip > ~{filebase}.vs.DCR2.face.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.heart.vs_zscore.uniq| sed -e  's/$/\tDCR2.heart/' | bgzip > ~{filebase}.vs.DCR2.heart.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.intestine.vs_zscore.uniq| sed -e  's/$/\tDCR2.intestine/' | bgzip > ~{filebase}.vs.DCR2.intestine.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.kidney.vs_zscore.uniq| sed -e  's/$/\tDCR2.kidney/' | bgzip > ~{filebase}.vs.DCR2.kidney.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.limb.vs_zscore.uniq| sed -e  's/$/\tDCR2.limb/' | bgzip > ~{filebase}.vs.DCR2.limb.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.liver.vs_zscore.uniq| sed -e  's/$/\tDCR2.liver/' | bgzip > ~{filebase}.vs.DCR2.liver.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.lung.vs_zscore.uniq| sed -e  's/$/\tDCR2.lung/' | bgzip > ~{filebase}.vs.DCR2.lung.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.muscle.vs_zscore.uniq| sed -e  's/$/\tDCR2.muscle/' | bgzip > ~{filebase}.vs.DCR2.muscle.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.neural_tube.vs_zscore.uniq| sed -e  's/$/\tDCR2.neural_tube/' | bgzip > ~{filebase}.vs.DCR2.neural_tube.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.placenta.vs_zscore.uniq| sed -e  's/$/\tDCR2.placenta/' | bgzip > ~{filebase}.vs.DCR2.placenta.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.skin.vs_zscore.uniq| sed -e  's/$/\tDCR2.skin/' | bgzip > ~{filebase}.vs.DCR2.skin.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DCR2.stomach.vs_zscore.uniq| sed -e  's/$/\tDCR2.stomach/' | bgzip > ~{filebase}.vs.DCR2.stomach.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.loeuf.phi.pts.tall_v0.1.vs_zscore.uniq | sed -e  's/$/\tabc.cognate.loeuf.phi.pts.tall_v0.1/' | bgzip > ~{filebase}.vs.abc.cognate.loeuf.phi.pts.tall_v0.1.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.loeuf_constraint.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.loeuf_constraint/' | bgzip > ~{filebase}.vs.abc.cognate.loeuf_constraint.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.loeuf_unconstraint.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.loeuf_unconstraint/' | bgzip > ~{filebase}.vs.abc.cognate.loeuf_unconstraint.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.phaplo_constraint.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.phaplo_constraint/' | bgzip > ~{filebase}.vs.abc.cognate.phaplo_constraint.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.phaplo_unconstraint.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.phaplo_unconstraint/' | bgzip > ~{filebase}.vs.abc.cognate.phaplo_unconstraint.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.ptriplo_constraint.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.ptriplo_constraint/' | bgzip > ~{filebase}.vs.abc.cognate.ptriplo_constraint.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.ptriplo_unconstraint.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.ptriplo_unconstraint/' | bgzip > ~{filebase}.vs.abc.cognate.ptriplo_unconstraint.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.loeuf_middle.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.loeuf_middle/' | bgzip > ~{filebase}.vs.abc.cognate.loeuf_middle.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.phaplo_middle.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.phaplo_middle/' | bgzip > ~{filebase}.vs.abc.cognate.phaplo_middle.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/abc.cognate.ptriplo_middle.vs_zscore.uniq| sed -e  's/$/\tabc.cognate.ptriplo_middle/' | bgzip > ~{filebase}.vs.abc.cognate.ptriplo_middle.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/chromHMM.Enh.txEnh.hg38.vs_zscore.uniq| sed -e  's/$/\tchromHMM.Enh.txEnh/' | bgzip > ~{filebase}.vs.chromHMM.Enh.txEnh.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/chromHMM.insulator.hg38.vs_zscore.uniq| sed -e  's/$/\tchromHMM.insulator/' | bgzip > ~{filebase}.vs.chromHMM.insulator.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/chromHMM.polycomb.hg38.vs_zscore.uniq| sed -e  's/$/\tchromHMM.polycomb/' | bgzip > ~{filebase}.vs.chromHMM.polycomb.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/chromHMM.promoter.TSS.hg38.vs_zscore.uniq| sed -e  's/$/\tchromHMM.promoter.TSS/' | bgzip > ~{filebase}.vs.chromHMM.promoter.TSS.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/enhancerAtlasV2.hg38.vs_zscore.uniq| sed -e  's/$/\tenhancerAtlasV2/' | bgzip > ~{filebase}.vs.enhancerAtlasV2.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/atac-seq_OCRs_2020.hg38.autosomes.reformat.vs_zscore.uniq| sed -e  's/$/\tatac-seq_OCRs_2020.autosomes/' | bgzip > ~{filebase}.vs.atac-seq_OCRs_2020.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/atac-seq_pREs_enhancers_2020.hg38.autosomes.reformat.vs_zscore.uniq| sed -e  's/$/\tatac-seq_pREs_enhancers_2020.autosomes/' | bgzip > ~{filebase}.vs.atac-seq_pREs_enhancers_2020.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Constrained_2493.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tConstrained_2493.autosomes/' | bgzip > ~{filebase}.vs.Constrained_2493.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/DNaseHS.UCSC.hg38.vs_zscore.uniq| sed -e  's/$/\tDNaseHS.UCSC/' | bgzip > ~{filebase}.vs.DNaseHS.UCSC.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.CTCF-only,CTCF-bound.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.CTCF-only,CTCF-bound/' | bgzip > ~{filebase}.vs.Encode3_ccREs.CTCF-only,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.dELS.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.dELS/' | bgzip > ~{filebase}.vs.Encode3_ccREs.dELS.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.dELS,CTCF-bound.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.dELS,CTCF-bound/' | bgzip > ~{filebase}.vs.Encode3_ccREs.dELS,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.DNase-H3K4me3.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.DNase-H3K4me3/' | bgzip > ~{filebase}.vs.Encode3_ccREs.DNase-H3K4me3.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.DNase-H3K4me3,CTCF-bound.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.DNase-H3K4me3,CTCF-bound/' | bgzip > ~{filebase}.vs.Encode3_ccREs.DNase-H3K4me3,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.pELS.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.pELS/' | bgzip > ~{filebase}.vs.Encode3_ccREs.pELS.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.pELS,CTCF-bound.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.pELS,CTCF-bound/' | bgzip > ~{filebase}.vs.Encode3_ccREs.pELS,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.PLS.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.PLS/' | bgzip > ~{filebase}.vs.Encode3_ccREs.PLS.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Encode3_ccREs.hg38.PLS,CTCF-bound.vs_zscore.uniq| sed -e  's/$/\tEncode3_ccREs.PLS,CTCF-bound/' | bgzip > ~{filebase}.vs.Encode3_ccREs.PLS,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.Bivalent.vs_zscore.uniq| sed -e  's/$/\tencode4.Bivalent/' | bgzip > ~{filebase}.vs.encode4.Bivalent.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.ConstitutiveHet.vs_zscore.uniq| sed -e  's/$/\tencode4.ConstitutiveHet/' | bgzip > ~{filebase}.vs.encode4.ConstitutiveHet.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.CTCF.vs_zscore.uniq| sed -e  's/$/\tencode4.CTCF/' | bgzip > ~{filebase}.vs.encode4.CTCF.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.Enhancer.vs_zscore.uniq| sed -e  's/$/\tencode4.Enhancer/' | bgzip > ~{filebase}.vs.encode4.Enhancer.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.EnhancerLow.vs_zscore.uniq| sed -e  's/$/\tencode4.EnhancerLow/' | bgzip > ~{filebase}.vs.encode4.EnhancerLow.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.FacultativeHet.vs_zscore.uniq| sed -e  's/$/\tencode4.FacultativeHet/' | bgzip > ~{filebase}.vs.encode4.FacultativeHet.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.K9K36.vs_zscore.uniq| sed -e  's/$/\tencode4.K9K36/' | bgzip > ~{filebase}.vs.encode4.K9K36.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.Promoter.vs_zscore.uniq| sed -e  's/$/\tencode4.Promoter/' | bgzip > ~{filebase}.vs.encode4.Promoter.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.PromoterFlanking.vs_zscore.uniq| sed -e  's/$/\tencode4.PromoterFlanking/' | bgzip > ~{filebase}.vs.encode4.PromoterFlanking.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/encode4.Transcribed.vs_zscore.uniq| sed -e  's/$/\tencode4.Transcribed/' | bgzip > ~{filebase}.vs.encode4.Transcribed.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/enhancer_DCR2v3.hg38.autosomes.reformat.vs_zscore.uniq| sed -e  's/$/\tenhancer_DCR2v3.autosomes/' | bgzip > ~{filebase}.vs.enhancer_DCR2v3.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.all_cells.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.all_cells/' | bgzip > ~{filebase}.vs.fantom5_enhancers.all_cells.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.all_organs.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.all_organs/' | bgzip > ~{filebase}.vs.fantom5_enhancers.all_organs.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.enhancer_clusters.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.enhancer_clusters/' | bgzip > ~{filebase}.vs.fantom5_enhancers.enhancer_clusters.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.permissive_enhancers.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.permissive_enhancers/' | bgzip > ~{filebase}.vs.fantom5_enhancers.permissive_enhancers.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.robust_enhancers.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.robust_enhancers/' | bgzip > ~{filebase}.vs.fantom5_enhancers.robust_enhancers.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.specific_enhancers_cells.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.specific_enhancers_cells/' | bgzip > ~{filebase}.vs.fantom5_enhancers.specific_enhancers_cells.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.specific_enhancers_organs.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.specific_enhancers_organs/' | bgzip > ~{filebase}.vs.fantom5_enhancers.specific_enhancers_organs.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.ubiquitous_enhancers_cells.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.ubiquitous_enhancers_cells/' | bgzip > ~{filebase}.vs.fantom5_enhancers.ubiquitous_enhancers_cells.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom5_enhancers.hg38.ubiquitous_enhancers_organs.vs_zscore.uniq| sed -e  's/$/\tfantom5_enhancers.ubiquitous_enhancers_organs/' | bgzip > ~{filebase}.vs.fantom5_enhancers.ubiquitous_enhancers_organs.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/fantom_enhancers.hg38.enhancers.siwei_paper.vs_zscore.uniq| sed -e  's/$/\tfantom_enhancers.enhancers.siwei_paper/' | bgzip > ~{filebase}.vs.fantom_enhancers.enhancers.siwei_paper.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/FetalBrain_2685.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tFetalBrain_2685.autosomes/' | bgzip > ~{filebase}.vs.FetalBrain_2685.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.all.noncoding.exon.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.all.noncoding.exon.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.all.noncoding.exon.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.all.noncoding.promoter.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.all.noncoding.promoter.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.all.noncoding.promoter.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.others.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.others.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.others.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.proteincoding.cds.extend10bp.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.proteincoding.cds.extend10bp.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.proteincoding.cds.extend10bp.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.proteincoding.promoter.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.proteincoding.promoter.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.proteincoding.promoter.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.proteincoding.utr3.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.proteincoding.utr3.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.proteincoding.utr3.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.proteincoding.utr5.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.proteincoding.utr5.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.proteincoding.utr5.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.pseudogenes.exon.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.pseudogenes.exon.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.pseudogenes.exon.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/gencode.v38.pseudogenes.promoter.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tgencode.v38.pseudogenes.promoter.autosomes/' | bgzip > ~{filebase}.vs.gencode.v38.pseudogenes.promoter.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/HAR.hg38.2023.vs_zscore.uniq| sed -e  's/$/\tHAR.2023/' | bgzip > ~{filebase}.vs.HAR.2023.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/human_fetal_brain_TF_footprints.hg38.autosomes.reformat.vs_zscore.uniq| sed -e  's/$/\thuman_fetal_brain_TF_footprints.autosomes/' | bgzip > ~{filebase}.vs.human_fetal_brain_TF_footprints.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/masked.regions.merged.LCRs.gaps.hg38.autosomes.reformat.sorted.vs_zscore.uniq| sed -e  's/$/\tmasked.regions.LCRs.gaps.autosomes/' | bgzip > ~{filebase}.vs.masked.regions.LCRs.gaps.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/NDD_664.autosomes.merged.sorted.vs_zscore.uniq| sed -e  's/$/\tNDD_664.autosomes/' | bgzip > ~{filebase}.vs.NDD_664.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Neuron_PROcap_enhancers.hg38.autosomes.reformat.vs_zscore.uniq| sed -e  's/$/\tNeuron_PROcap_enhancers.autosomes/' | bgzip > ~{filebase}.vs.Neuron_PROcap_enhancers.autosomes.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/SE_package.lft38.vs_zscore.uniq| sed -e  's/$/\tSE_package.lft38/' | bgzip > ~{filebase}.vs.SE_package.lft38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/Shi_Lab.TAD_Boundry.hg38.vs_zscore.uniq| sed -e  's/$/\tShi_Lab.TAD_Boundry/' | bgzip > ~{filebase}.vs.Shi_Lab.TAD_Boundry.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/TFChip.UCSC.hg38.vs_zscore.uniq| sed -e  's/$/\tTFChip.UCSC/' | bgzip > ~{filebase}.vs.TFChip.UCSC.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a <(zcat ~{SV_file} ) -b noncoding/UCNE_coord.hg38.vs_zscore.uniq| sed -e  's/$/\tUCNE_coord/' | bgzip > ~{filebase}.vs.UCNE_coord.over_50perc_ovr.bed.gz

            zcat ~{filebase}.vs.abc.cognate.loeuf_constraint.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.loeuf_constraint.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.loeuf_middle.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.loeuf_middle.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.loeuf.phi.pts.tall_v0.1.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.loeuf.phi.pts.tall_v0.1.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.loeuf_unconstraint.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.loeuf_unconstraint.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.phaplo_constraint.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.phaplo_constraint.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.phaplo_middle.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.phaplo_middle.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.phaplo_unconstraint.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.phaplo_unconstraint.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.ptriplo_constraint.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.ptriplo_constraint.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.ptriplo_middle.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.ptriplo_middle.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.cognate.ptriplo_unconstraint.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.abc.cognate.ptriplo_unconstraint.over_50perc_ovr.stat
            zcat ~{filebase}.vs.atac-seq_OCRs_2020.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.atac-seq_OCRs_2020.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.atac-seq_pREs_enhancers_2020.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.atac-seq_pREs_enhancers_2020.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.Enh.txEnh.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.chromHMM.Enh.txEnh.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.insulator.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.chromHMM.insulator.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.polycomb.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.chromHMM.polycomb.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.promoter.TSS.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.chromHMM.promoter.TSS.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Constrained_2493.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Constrained_2493.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.adrenal_glands.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.adrenal_glands.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.brain.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.brain.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.eye.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.eye.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.face.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.face.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.general.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.general.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.heart.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.heart.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.intestine.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.intestine.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.kidney.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.kidney.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.limb.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.limb.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.liver.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.liver.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.lung.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.lung.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.muscle.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.muscle.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.neural_tube.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.neural_tube.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.placenta.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.placenta.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.skin.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.skin.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.stomach.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DCR2.stomach.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DNaseHS.UCSC.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.DNaseHS.UCSC.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.CTCF-only,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.CTCF-only,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.dELS,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.dELS,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.dELS.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.dELS.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.DNase-H3K4me3,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.DNase-H3K4me3,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.DNase-H3K4me3.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.DNase-H3K4me3.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.pELS,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.pELS,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.pELS.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.pELS.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.PLS,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.PLS,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.PLS.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Encode3_ccREs.PLS.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Bivalent.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.Bivalent.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.ConstitutiveHet.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.ConstitutiveHet.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.CTCF.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.CTCF.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.EnhancerLow.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.EnhancerLow.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Enhancer.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.Enhancer.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.FacultativeHet.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.FacultativeHet.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.K9K36.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.K9K36.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.PromoterFlanking.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.PromoterFlanking.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Promoter.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.Promoter.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Transcribed.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.encode4.Transcribed.over_50perc_ovr.stat
            zcat ~{filebase}.vs.enhancerAtlasV2.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.enhancerAtlasV2.over_50perc_ovr.stat
            zcat ~{filebase}.vs.enhancer_DCR2v3.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.enhancer_DCR2v3.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.all_cells.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.all_cells.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.all_organs.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.all_organs.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.enhancer_clusters.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.enhancer_clusters.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.permissive_enhancers.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.permissive_enhancers.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.robust_enhancers.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.robust_enhancers.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.specific_enhancers_cells.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.specific_enhancers_cells.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.specific_enhancers_organs.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.specific_enhancers_organs.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.ubiquitous_enhancers_cells.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.ubiquitous_enhancers_cells.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.ubiquitous_enhancers_organs.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom5_enhancers.ubiquitous_enhancers_organs.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom_enhancers.enhancers.siwei_paper.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.fantom_enhancers.enhancers.siwei_paper.over_50perc_ovr.stat
            zcat ~{filebase}.vs.FetalBrain_2685.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.FetalBrain_2685.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.all.noncoding.exon.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.all.noncoding.exon.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.all.noncoding.promoter.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.all.noncoding.promoter.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.others.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.others.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.proteincoding.cds.extend10bp.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.proteincoding.cds.extend10bp.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.proteincoding.promoter.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.proteincoding.promoter.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.proteincoding.utr3.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.proteincoding.utr3.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.proteincoding.utr5.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.proteincoding.utr5.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.pseudogenes.exon.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.pseudogenes.exon.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.gencode.v38.pseudogenes.promoter.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.gencode.v38.pseudogenes.promoter.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.HAR.2023.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.HAR.2023.over_50perc_ovr.stat
            zcat ~{filebase}.vs.human_fetal_brain_TF_footprints.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.human_fetal_brain_TF_footprints.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.masked.regions.LCRs.gaps.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.masked.regions.LCRs.gaps.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.NDD_664.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.NDD_664.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Neuron_PROcap_enhancers.autosomes.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Neuron_PROcap_enhancers.autosomes.over_50perc_ovr.stat
            zcat ~{filebase}.vs.SE_package.lft38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.SE_package.lft38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Shi_Lab.TAD_Boundry.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.Shi_Lab.TAD_Boundry.over_50perc_ovr.stat
            zcat ~{filebase}.vs.TFChip.UCSC.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.TFChip.UCSC.over_50perc_ovr.stat
            zcat ~{filebase}.vs.UCNE_coord.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.UCNE_coord.over_50perc_ovr.stat
            zcat ~{filebase}.vs.vista_cores.neg.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.vista_cores.neg.over_50perc_ovr.stat
            zcat ~{filebase}.vs.vista_cores.pos.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c > ~{filebase}.vs.vista_cores.pos.over_50perc_ovr.stat

            ls *.over_50perc_ovr.stat > SV_vs_nc_stat_list.tsv
            Rscript ../ncas_docker/src/integrate_stats_across_nc_elements.R -i SV_vs_nc_stat_list.tsv -o ~{filebase}.vs.nc_elements.integrated
            bgzip ~{filebase}.vs.nc_elements.integrated

    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: rdpesr_benchmark_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}



