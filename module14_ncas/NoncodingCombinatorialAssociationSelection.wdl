version 1.0

import "Structs.wdl"
import "SVvsConservative.wdl" as SVvsConservative

workflow NoncodingCombinatorialAssociationSelection {
    input{
        Int permutation_rounds
        File SV_sites_file
        File SV_function_file
        File SVID_genomic_context
        File src_tar
        File ref_tar
        File gene_tar
        File gencode_tar
        File conserve_tar
        File noncoding_tar
        String prefix
        File contig_file
        String sv_base_mini_docker

    }

    call CalculateAPS{
        input:
            ref_tar = ref_tar,
            src_tar = src_tar,
            SV_sites_file = SV_sites_file,
            SV_function_file = SV_function_file,
            SVID_genomic_context = SVID_genomic_context,
            prefix = prefix,
            sv_base_mini_docker = sv_base_mini_docker
    }

    scatter(i in range(permutation_rounds)){
        call GeneratePermutatedSVs{
            input:
                permu = i, 
                ref_tar = ref_tar,
                src_tar = src_tar,
                SV_sites_file = SV_sites_file,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVvsConservative.SV_vs_Conservative as SV_vs_Conservative{
            input:
                src_tar = src_tar,
                SV_file = GeneratePermutatedSVs.permutated_SV,
                conserve_tar = conserve_tar,
                contig_file = contig_file,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVvsGencode{
            input: 
                src_tar = src_tar,
                gencode_tar = gencode_tar,
                SV_file = GeneratePermutatedSVs.permutated_SV,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVvsNoncoding{
            input:
                src_tar = src_tar,
                noncoding_tar = noncoding_tar,
                SV_file = GeneratePermutatedSVs.permutated_SV,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call SVvsGene{
            input:
                src_tar = src_tar,
                gene_tar = gene_tar,
                SV_file = GeneratePermutatedSVs.permutated_SV,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call GenerateNcasMetrics{
            input:
                permu = i,
                src_tar = src_tar,
                svid_aps = CalculateAPS.svid_aps,
                svid_genomic_context = SVID_genomic_context,
                sv_file_real = SV_sites_file,
                sv_file_permu = GeneratePermutatedSVs.permutated_SV,
                sv_vs_gencode = SVvsGencode.SV_vs_gencode,
                sv_vs_conserve = SV_vs_Conservative.SV_vs_conserved,
                sv_vs_noncoding = SVvsNoncoding.SV_vs_noncoding,
                sv_vs_gene = SVvsGene.SV_vs_trans,
                sv_vs_coding = SVvsGene.SV_vs_cds,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call CalcuNcasStat{
            input:
                permu = i,
                src_tar = src_tar,
                ncas_rdata = GenerateNcasMetrics.ncas_rdata,
                sv_base_mini_docker = sv_base_mini_docker
        }
    }
    output{
        Array[File] ncas_data_list = GenerateNcasMetrics.ncas_rdata
        Array[File] ncas_stat = CalcuNcasStat.ncas_stat
    }
}

task CalculateAPS{
    input{
        File src_tar
        File ref_tar
        File SV_sites_file
        File SV_function_file
        File SVID_genomic_context
        String prefix
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
        File svid_aps = "~{prefix}.SVID_aps"
        File aps_deviation = "~{prefix}.aps_deviation.tsv"
        File aps_figures = "~{prefix}.aps.tar.gz"
    }

    String filebase = basename(SV_sites_file,".gz")

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{ref_tar} ./
        tar zxvf ref.tar.gz 
        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 

        Rscript ./src/calculate_APS.R \
        --SV_color ref/SV_colors.tsv \
        --genomic_context ~{SVID_genomic_context} \
        --SV_function_prediction ~{SV_function_file} \
        --SV_info ~{SV_sites_file} \
        --prefix ~{prefix} \
        -f "~{prefix}.aps"

        tar czvf ~{prefix}.aps.tar.gz ~{prefix}.aps/
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

task CalcuNcasStat{
    input{
        String permu
        File src_tar
        File ncas_rdata
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 10, 
        disk_gb: 20,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File ncas_stat = "ncas_stat.permu_~{permu}.tar.gz"
    }

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 

        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t ALL -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t DEL -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t DUP -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INV -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t CPX -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS:ME:ALU -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS:ME:SVA -g noncoding
        Rscript ./src/calculate_cwas_statistics.R -d ~{ncas_rdata} -a permu_~{permu}.stat -t INS:ME:LINE1 -g noncoding
        mkdir ncas_stat.permu_~{permu}
        mv *.stat ncas_stat.permu_~{permu}/
        tar czvf ncas_stat.permu_~{permu}.tar.gz ncas_stat.permu_~{permu}/
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

task GeneratePermutatedSVs{
    input{
        Int permu
        File src_tar
        File ref_tar
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
        File permutated_SV = "~{filebase}.permu_~{permu}.corrected.gz"
    }

    String filebase = basename(SV_sites_file,".gz")

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{ref_tar} ./
        tar zxvf ref.tar.gz 

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz 

        Rscript ./src/generate_SV_permutations.R -p ~{permu} -i ~{SV_sites_file} -o ~{filebase}.permu_~{permu} -g ref/hg38.genome.tsv
        bgzip ~{filebase}.permu_~{permu}

        Rscript ./src/correct_permutated_SVs.R -i ~{filebase}.permu_~{permu}.gz -r ~{SV_sites_file}
        bgzip ~{filebase}.permu_~{permu}.corrected

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

task GenerateNcasMetrics{
    input{
        Int permu
        File src_tar
        File sv_file_real
        File sv_file_permu
        File sv_vs_gencode
        File sv_vs_conserve    
        File sv_vs_noncoding
        File sv_vs_gene
        File sv_vs_coding
        File svid_aps
        File svid_genomic_context

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
        File ncas_rdata = "~{filebase}.rData"
    }

    String filebase = basename(sv_file_permu,".gz")

    command <<<
        set -Eeuo pipefail

            gsutil cp ~{src_tar} ./
            tar zxvf src.tar.gz

            Rscript ./src/generate_cwas_metrics.R \
                    --permu "permu_~{permu}" \
                    --sv_file_real ~{sv_file_real} \
                    --sv_file_permu ~{sv_file_permu} \
                    --sv_vs_gencode ~{sv_vs_gencode} \
                    --sv_vs_conserve ~{sv_vs_conserve} \
                    --sv_vs_noncoding ~{sv_vs_noncoding} \
                    --sv_vs_gene ~{sv_vs_gene} \
                    --sv_vs_coding ~{sv_vs_coding} \
                    --SVID_aps  ~{svid_aps} \
                    --SVID_genomic_context  ~{svid_genomic_context} \
                    --output ~{filebase}.rData
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

task SVvsGene{
    input{
        File src_tar
        File SV_file
        File gene_tar
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
        File SV_vs_cds = "~{filebase}.coding_interruptive.bed.gz"
        File SV_vs_trans = "~{filebase}.genic.bed.gz"
    }

    String filebase = basename(SV_file,".gz")

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{gene_tar} ./
        tar zxvf gene.tar.gz

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz

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
        docker: sv_base_mini_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task SVvsGencode{
    input{
        File SV_file
        File src_tar
        File gencode_tar
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
        File SV_vs_gencode = "~{filebase}.vs.gencode.integrated.gz"
    }

    String filebase = basename(SV_file,".gz")

    command <<<
        set -Eeuo pipefail

        gsutil cp ~{gencode_tar} ./
        tar zxvf gencode.tar.gz

        gsutil cp ~{src_tar} ./
        tar zxvf src.tar.gz

        zcat ~{SV_file} | cut -f1-4 > input_SVs.bed
        zcat gencode/gencode.v39.integrated.bed.gz | cut -f1-4,8 > gencode.bed
        bedtools intersect -wo -a input_SVs.bed -b gencode.bed > ~{filebase}.vs.gencode.bed
        cat gencode/header ~{filebase}.vs.gencode.bed | bgzip >  ~{filebase}.vs.gencode.bed.gz
        Rscript ./src/integrate_SV_vs_gencode.R -i ~{filebase}.vs.gencode.bed.gz  -o ~{filebase}.vs.gencode.integrated
        bgzip ~{filebase}.vs.gencode.integrated
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

task SVvsNoncoding{
    input{
        File SV_file
        File src_tar
        File noncoding_tar
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
        File SV_vs_noncoding = "~{filebase}.vs.nc_elements.integrated.gz"
    }

    String filebase = basename(SV_file,".gz")

    command <<<
        set -Eeuo pipefail

            gsutil cp ~{noncoding_tar} ./
            tar zxvf noncoding.tar.gz

            gsutil cp ~{src_tar} ./
            tar zxvf src.tar.gz

            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/atac-seq_OCRs_2020.hg38.autosomes.reformat.vs_zscore.uniq  | sed -e  's/$/\tatac-seq_OCRs_2020.hg38.autosomes.reformat/' | bgzip >  ~{filebase}.vs.atac-seq_OCRs_2020.hg38.autosomes.reformat.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/atac-seq_pREs_enhancers_2020.hg38.autosomes.reformat.vs_zscore.uniq  | sed -e  's/$/\tatac-seq_pREs_enhancers_2020.hg38.autosomes.reformat/' | bgzip >  ~{filebase}.vs.atac-seq_pREs_enhancers_2020.hg38.autosomes.reformat.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DNaseHS.UCSC.hg38.vs_zscore.uniq  | sed -e  's/$/\tDNaseHS.UCSC.hg38/' | bgzip >  ~{filebase}.vs.DNaseHS.UCSC.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/enhancer_DCR2v3.hg38.autosomes.reformat.vs_zscore.uniq  | sed -e  's/$/\tenhancer_DCR2v3.hg38.autosomes.reformat/' | bgzip >  ~{filebase}.vs.enhancer_DCR2v3.hg38.autosomes.reformat.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/GW18_cortex_enh.vs_zscore.uniq  | sed -e  's/$/\tGW18_cortex_enh/' | bgzip >  ~{filebase}.vs.GW18_cortex_enh.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/GW18_cortex_sil.vs_zscore.uniq  | sed -e  's/$/\tGW18_cortex_sil/' | bgzip >  ~{filebase}.vs.GW18_cortex_sil.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/hg38_gene_deserts.vs_zscore.uniq  | sed -e  's/$/\thg38_gene_deserts/' | bgzip >  ~{filebase}.vs.hg38_gene_deserts.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/human_fetal_brain_TF_footprints.hg38.autosomes.reformat.vs_zscore.uniq  | sed -e  's/$/\thuman_fetal_brain_TF_footprints.hg38.autosomes.reformat/' | bgzip >  ~{filebase}.vs.human_fetal_brain_TF_footprints.hg38.autosomes.reformat.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/masked.regions.merged.LCRs.gaps.hg38.autosomes.reformat.sorted.vs_zscore.uniq  | sed -e  's/$/\tmasked.regions.merged.LCRs.gaps.hg38.autosomes.reformat.sorted/' | bgzip >  ~{filebase}.vs.masked.regions.merged.LCRs.gaps.hg38.autosomes.reformat.sorted.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/Neuron_PROcap_enhancers.hg38.autosomes.reformat.vs_zscore.uniq  | sed -e  's/$/\tNeuron_PROcap_enhancers.hg38.autosomes.reformat/' | bgzip >  ~{filebase}.vs.Neuron_PROcap_enhancers.hg38.autosomes.reformat.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/Shi_Lab.TAD_Boundry.hg38.vs_zscore.uniq  | sed -e  's/$/\tShi_Lab.TAD_Boundry.hg38/' | bgzip >  ~{filebase}.vs.Shi_Lab.TAD_Boundry.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/TFChip.UCSC.hg38.vs_zscore.uniq  | sed -e  's/$/\tTFChip.UCSC.hg38/' | bgzip >  ~{filebase}.vs.TFChip.UCSC.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/tri1_brain_atac.vs_zscore.uniq  | sed -e  's/$/\ttri1_brain_atac/' | bgzip >  ~{filebase}.vs.tri1_brain_atac.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/vista_cores.neg.vs_zscore.uniq  | sed -e  's/$/\tvista_cores.neg/' | bgzip >  ~{filebase}.vs.vista_cores.neg.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/vista_cores.pos.vs_zscore.uniq  | sed -e  's/$/\tvista_cores.pos/' | bgzip >  ~{filebase}.vs.vista_cores.pos.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_A549_treated_with_ethanol_0_02_percent_for_1_hour_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_A549_treated_with_ethanol_0_02_percent_for_1_hour_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_A549_treated_with_ethanol_0_02_percent_for_1_hour_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_A673_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_A673_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_A673_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_adipose_tissue_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_adipose_tissue_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_adipose_tissue_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_adrenal_gland_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_adrenal_gland_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_adrenal_gland_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_adrenal_gland_fetal_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_adrenal_gland_fetal_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_adrenal_gland_fetal_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_astrocyte_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_astrocyte_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_astrocyte_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_B_cell_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_B_cell_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_B_cell_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_bipolar_neuron_from_iPSC_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_bipolar_neuron_from_iPSC_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_bipolar_neuron_from_iPSC_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_BJAB_anti_IgM_anti_CD40_4hr_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_BJAB_anti_IgM_anti_CD40_4hr_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_BJAB_anti_IgM_anti_CD40_4hr_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_BJAB_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_BJAB_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_BJAB_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_body_of_pancreas_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_body_of_pancreas_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_body_of_pancreas_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_breast_epithelium_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_breast_epithelium_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_breast_epithelium_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_brite_adipose_Loft2014.vs_zscore.uniq  | sed -e  's/$/\tABC_brite_adipose_Loft2014/' | bgzip >  ~{filebase}.vs.ABC_brite_adipose_Loft2014.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_cardiac_muscle_cell_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_cardiac_muscle_cell_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_cardiac_muscle_cell_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocytes_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocytes_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocytes_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_BG_1h_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_BG_1h_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_1h_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_BG_4h_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_BG_4h_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_4h_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_BG_d1_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_BG_d1_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_d1_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_BG_d6_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_BG_d6_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_d6_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_LPS_1h_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_LPS_1h_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_1h_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_LPS_d1_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_LPS_d1_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_d1_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_LPS_d6_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_LPS_d6_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_d6_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_RPMI_1h_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_RPMI_1h_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_1h_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_RPMI_4h_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_RPMI_4h_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_4h_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_RPMI_d1_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_RPMI_d1_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_d1_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD14_positive_monocyte_treated_with_RPMI_d6_Novakovic2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD14_positive_monocyte_treated_with_RPMI_d6_Novakovic2016/' | bgzip >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_d6_Novakovic2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD19_positive_B_cell_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_CD19_positive_B_cell_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_CD19_positive_B_cell_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD34_positive_mobilized_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_CD34_positive_mobilized_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_CD34_positive_mobilized_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD3_positive_T_cell_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_CD3_positive_T_cell_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_CD3_positive_T_cell_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD4_positive_helper_T_cell_Corces2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD4_positive_helper_T_cell_Corces2016/' | bgzip >  ~{filebase}.vs.ABC_CD4_positive_helper_T_cell_Corces2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD4_positive_helper_T_cell_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_CD4_positive_helper_T_cell_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_CD4_positive_helper_T_cell_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD56_positive_natural_killer_cells_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_CD56_positive_natural_killer_cells_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_CD56_positive_natural_killer_cells_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD8_positive_alpha_beta_T_cell_Corces2016.vs_zscore.uniq  | sed -e  's/$/\tABC_CD8_positive_alpha_beta_T_cell_Corces2016/' | bgzip >  ~{filebase}.vs.ABC_CD8_positive_alpha_beta_T_cell_Corces2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_CD8_positive_alpha_beta_T_cell_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_CD8_positive_alpha_beta_T_cell_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_CD8_positive_alpha_beta_T_cell_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_coronary_artery_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_coronary_artery_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_coronary_artery_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_coronary_artery_smooth_muscle_cell_Miller2016.vs_zscore.uniq  | sed -e  's/$/\tABC_coronary_artery_smooth_muscle_cell_Miller2016/' | bgzip >  ~{filebase}.vs.ABC_coronary_artery_smooth_muscle_cell_Miller2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_dendritic_cell_treated_with_Lipopolysaccharide_0_ng_mL_for_0_hour_Garber2017.vs_zscore.uniq  | sed -e  's/$/\tABC_dendritic_cell_treated_with_Lipopolysaccharide_0_ng_mL_for_0_hour_Garber2017/' | bgzip >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_0_ng_mL_for_0_hour_Garber2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_1_hour_Garber2017.vs_zscore.uniq  | sed -e  's/$/\tABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_1_hour_Garber2017/' | bgzip >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_1_hour_Garber2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_2_hour_Garber2017.vs_zscore.uniq  | sed -e  's/$/\tABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_2_hour_Garber2017/' | bgzip >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_2_hour_Garber2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_30_minute_Garber2017.vs_zscore.uniq  | sed -e  's/$/\tABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_30_minute_Garber2017/' | bgzip >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_30_minute_Garber2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_4_hour_Garber2017.vs_zscore.uniq  | sed -e  's/$/\tABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_4_hour_Garber2017/' | bgzip >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_4_hour_Garber2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_6_hour_Garber2017.vs_zscore.uniq  | sed -e  's/$/\tABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_6_hour_Garber2017/' | bgzip >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_6_hour_Garber2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_endothelial_cell_of_umbilical_vein_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_endothelial_cell_of_umbilical_vein_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_endothelial_cell_of_umbilical_vein_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_epithelial_cell_of_prostate_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_epithelial_cell_of_prostate_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_epithelial_cell_of_prostate_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_erythroblast_Corces2016.vs_zscore.uniq  | sed -e  's/$/\tABC_erythroblast_Corces2016/' | bgzip >  ~{filebase}.vs.ABC_erythroblast_Corces2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_fibroblast_of_arm_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_fibroblast_of_arm_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_fibroblast_of_arm_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_fibroblast_of_dermis_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_fibroblast_of_dermis_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_fibroblast_of_dermis_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_fibroblast_of_lung_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_fibroblast_of_lung_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_fibroblast_of_lung_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_foreskin_fibroblast_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_foreskin_fibroblast_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_foreskin_fibroblast_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_gastrocnemius_medialis_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_gastrocnemius_medialis_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_gastrocnemius_medialis_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_GM12878_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_GM12878_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_GM12878_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_H1_BMP4_Derived_Mesendoderm_Cultured_Cells_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_H1_BMP4_Derived_Mesendoderm_Cultured_Cells_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_H1_BMP4_Derived_Mesendoderm_Cultured_Cells_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_H1_BMP4_Derived_Trophoblast_Cultured_Cells_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_H1_BMP4_Derived_Trophoblast_Cultured_Cells_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_H1_BMP4_Derived_Trophoblast_Cultured_Cells_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_H1_Derived_Mesenchymal_Stem_Cells_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_H1_Derived_Mesenchymal_Stem_Cells_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_H1_Derived_Mesenchymal_Stem_Cells_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_H1_Derived_Neuronal_Progenitor_Cultured_Cells_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_H1_Derived_Neuronal_Progenitor_Cultured_Cells_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_H1_Derived_Neuronal_Progenitor_Cultured_Cells_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_H1_hESC_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_H1_hESC_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_H1_hESC_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_H7.vs_zscore.uniq  | sed -e  's/$/\tABC_H7/' | bgzip >  ~{filebase}.vs.ABC_H7.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_H9_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_H9_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_H9_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_HAP1.vs_zscore.uniq  | sed -e  's/$/\tABC_HAP1/' | bgzip >  ~{filebase}.vs.ABC_HAP1.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_HCT116_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_HCT116_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_HCT116_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_heart_ventricle_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_heart_ventricle_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_heart_ventricle_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_HeLa_S3_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_HeLa_S3_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_HeLa_S3_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_hepatocyte_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_hepatocyte_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_hepatocyte_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_HepG2_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_HepG2_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_HepG2_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_HT29.vs_zscore.uniq  | sed -e  's/$/\tABC_HT29/' | bgzip >  ~{filebase}.vs.ABC_HT29.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_IMR90_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_IMR90_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_IMR90_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_induced_pluripotent_stem_cell_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_induced_pluripotent_stem_cell_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_induced_pluripotent_stem_cell_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_iPS_DF_19_11_Cell_Line_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_iPS_DF_19_11_Cell_Line_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_iPS_DF_19_11_Cell_Line_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_Jurkat_anti_CD3_PMA_4hr_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_Jurkat_anti_CD3_PMA_4hr_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_Jurkat_anti_CD3_PMA_4hr_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_Jurkat_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_Jurkat_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_Jurkat_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_K562_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_K562_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_K562_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_Karpas_422_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_Karpas_422_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_Karpas_422_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_keratinocyte_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_keratinocyte_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_keratinocyte_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_large_intestine_fetal_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_large_intestine_fetal_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_large_intestine_fetal_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_liver_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_liver_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_liver_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_LNCAP.vs_zscore.uniq  | sed -e  's/$/\tABC_LNCAP/' | bgzip >  ~{filebase}.vs.ABC_LNCAP.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/abc.loeuf.0.15.hg38.vs_zscore.uniq  | sed -e  's/$/\tabc.loeuf.0.15.hg38/' | bgzip >  ~{filebase}.vs.abc.loeuf.0.15.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_LoVo.vs_zscore.uniq  | sed -e  's/$/\tABC_LoVo/' | bgzip >  ~{filebase}.vs.ABC_LoVo.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_mammary_epithelial_cell_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_mammary_epithelial_cell_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_mammary_epithelial_cell_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_MCF10A_Ji2017.vs_zscore.uniq  | sed -e  's/$/\tABC_MCF10A_Ji2017/' | bgzip >  ~{filebase}.vs.ABC_MCF10A_Ji2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_MCF10A_treated_with_TAM24hr_Ji2017.vs_zscore.uniq  | sed -e  's/$/\tABC_MCF10A_treated_with_TAM24hr_Ji2017/' | bgzip >  ~{filebase}.vs.ABC_MCF10A_treated_with_TAM24hr_Ji2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_MCF_7_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_MCF_7_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_MCF_7_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_MDA_MB_231.vs_zscore.uniq  | sed -e  's/$/\tABC_MDA_MB_231/' | bgzip >  ~{filebase}.vs.ABC_MDA_MB_231.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_megakaryocyte_erythroid_progenitor_Corces2016.vs_zscore.uniq  | sed -e  's/$/\tABC_megakaryocyte_erythroid_progenitor_Corces2016/' | bgzip >  ~{filebase}.vs.ABC_megakaryocyte_erythroid_progenitor_Corces2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/abc.merge.hg38.vs_zscore.uniq  | sed -e  's/$/\tabc.merge.hg38/' | bgzip >  ~{filebase}.vs.abc.merge.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_MM_1S_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_MM_1S_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_MM_1S_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_muscle_of_leg_fetal_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_muscle_of_leg_fetal_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_muscle_of_leg_fetal_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_muscle_of_trunk_fetal_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_muscle_of_trunk_fetal_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_muscle_of_trunk_fetal_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_myotube_originated_from_skeletal_muscle_myoblast_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_myotube_originated_from_skeletal_muscle_myoblast_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_myotube_originated_from_skeletal_muscle_myoblast_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_natural_killer_cell_Corces2016.vs_zscore.uniq  | sed -e  's/$/\tABC_natural_killer_cell_Corces2016/' | bgzip >  ~{filebase}.vs.ABC_natural_killer_cell_Corces2016.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_NCCIT.vs_zscore.uniq  | sed -e  's/$/\tABC_NCCIT/' | bgzip >  ~{filebase}.vs.ABC_NCCIT.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_OCI_LY7_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_OCI_LY7_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_OCI_LY7_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_osteoblast_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_osteoblast_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_osteoblast_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_ovary_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_ovary_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_ovary_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_Panc1_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_Panc1_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_Panc1_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_pancreas_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_pancreas_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_pancreas_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_PC_9_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_PC_9_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_PC_9_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/abc.pHaplo.0.9.hg38.vs_zscore.uniq  | sed -e  's/$/\tabc.pHaplo.0.9.hg38/' | bgzip >  ~{filebase}.vs.abc.pHaplo.0.9.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_placenta_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_placenta_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_placenta_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_psoas_muscle_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_psoas_muscle_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_psoas_muscle_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/abc.pTriplo.0.84.hg38.vs_zscore.uniq  | sed -e  's/$/\tabc.pTriplo.0.84.hg38/' | bgzip >  ~{filebase}.vs.abc.pTriplo.0.84.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_sigmoid_colon_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_sigmoid_colon_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_sigmoid_colon_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_skeletal_muscle_myoblast_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_skeletal_muscle_myoblast_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_skeletal_muscle_myoblast_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_SK_N_SH_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_SK_N_SH_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_SK_N_SH_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_small_intestine_fetal_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_small_intestine_fetal_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_small_intestine_fetal_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_spinal_cord_fetal_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_spinal_cord_fetal_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_spinal_cord_fetal_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_spleen_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_spleen_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_spleen_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_stomach_fetal_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_stomach_fetal_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_stomach_fetal_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_stomach_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_stomach_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_stomach_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_T_cell_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_T_cell_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_T_cell_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP1_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_THP1_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_THP1_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP1_LPS_4hr_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_THP1_LPS_4hr_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_THP1_LPS_4hr_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_1_macrophage_VanBortle2017.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_1_macrophage_VanBortle2017/' | bgzip >  ~{filebase}.vs.ABC_THP_1_macrophage_VanBortle2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_1_monocyte_VanBortle2017.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_1_monocyte_VanBortle2017/' | bgzip >  ~{filebase}.vs.ABC_THP_1_monocyte_VanBortle2017.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_0h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_0h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_0h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_120h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_120h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_120h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_12h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_12h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_12h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_1h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_1h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_1h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_24h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_24h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_24h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_2h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_2h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_2h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_48h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_48h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_48h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_6h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_6h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_6h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_72h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_72h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_72h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_THP_pmaLPS_ATAC_96h.vs_zscore.uniq  | sed -e  's/$/\tABC_THP_pmaLPS_ATAC_96h/' | bgzip >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_96h.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_thymus_fetal_Roadmap.vs_zscore.uniq  | sed -e  's/$/\tABC_thymus_fetal_Roadmap/' | bgzip >  ~{filebase}.vs.ABC_thymus_fetal_Roadmap.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_thyroid_gland_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_thyroid_gland_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_thyroid_gland_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_transverse_colon_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_transverse_colon_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_transverse_colon_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_trophoblast_cell_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_trophoblast_cell_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_trophoblast_cell_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_U937_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_U937_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_U937_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_U937_LPS_4hr_Engreitz.vs_zscore.uniq  | sed -e  's/$/\tABC_U937_LPS_4hr_Engreitz/' | bgzip >  ~{filebase}.vs.ABC_U937_LPS_4hr_Engreitz.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_uterus_ENCODE.vs_zscore.uniq  | sed -e  's/$/\tABC_uterus_ENCODE/' | bgzip >  ~{filebase}.vs.ABC_uterus_ENCODE.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/abc/ABC_white_adipose_Loft2014.vs_zscore.uniq  | sed -e  's/$/\tABC_white_adipose_Loft2014/' | bgzip >  ~{filebase}.vs.ABC_white_adipose_Loft2014.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/chromHMM/chromHMM.Enh.txEnh.hg38.vs_zscore.uniq  | sed -e  's/$/\tchromHMM.Enh.txEnh.hg38/' | bgzip >  ~{filebase}.vs.chromHMM.Enh.txEnh.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/chromHMM/chromHMM.insulator.hg38.vs_zscore.uniq  | sed -e  's/$/\tchromHMM.insulator.hg38/' | bgzip >  ~{filebase}.vs.chromHMM.insulator.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/chromHMM/chromHMM.polycomb.hg38.vs_zscore.uniq  | sed -e  's/$/\tchromHMM.polycomb.hg38/' | bgzip >  ~{filebase}.vs.chromHMM.polycomb.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/chromHMM/chromHMM.promoter.TSS.hg38.vs_zscore.uniq  | sed -e  's/$/\tchromHMM.promoter.TSS.hg38/' | bgzip >  ~{filebase}.vs.chromHMM.promoter.TSS.hg38.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.adrenal_glands.vs_zscore.uniq  | sed -e  's/$/\tDCR2.adrenal_glands/' | bgzip >  ~{filebase}.vs.DCR2.adrenal_glands.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.brain.vs_zscore.uniq  | sed -e  's/$/\tDCR2.brain/' | bgzip >  ~{filebase}.vs.DCR2.brain.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.eye.vs_zscore.uniq  | sed -e  's/$/\tDCR2.eye/' | bgzip >  ~{filebase}.vs.DCR2.eye.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.face.vs_zscore.uniq  | sed -e  's/$/\tDCR2.face/' | bgzip >  ~{filebase}.vs.DCR2.face.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.general.vs_zscore.uniq  | sed -e  's/$/\tDCR2.general/' | bgzip >  ~{filebase}.vs.DCR2.general.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.heart.vs_zscore.uniq  | sed -e  's/$/\tDCR2.heart/' | bgzip >  ~{filebase}.vs.DCR2.heart.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.intestine.vs_zscore.uniq  | sed -e  's/$/\tDCR2.intestine/' | bgzip >  ~{filebase}.vs.DCR2.intestine.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.kidney.vs_zscore.uniq  | sed -e  's/$/\tDCR2.kidney/' | bgzip >  ~{filebase}.vs.DCR2.kidney.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.limb.vs_zscore.uniq  | sed -e  's/$/\tDCR2.limb/' | bgzip >  ~{filebase}.vs.DCR2.limb.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.liver.vs_zscore.uniq  | sed -e  's/$/\tDCR2.liver/' | bgzip >  ~{filebase}.vs.DCR2.liver.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.lung.vs_zscore.uniq  | sed -e  's/$/\tDCR2.lung/' | bgzip >  ~{filebase}.vs.DCR2.lung.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.muscle.vs_zscore.uniq  | sed -e  's/$/\tDCR2.muscle/' | bgzip >  ~{filebase}.vs.DCR2.muscle.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.neural_tube.vs_zscore.uniq  | sed -e  's/$/\tDCR2.neural_tube/' | bgzip >  ~{filebase}.vs.DCR2.neural_tube.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.placenta.vs_zscore.uniq  | sed -e  's/$/\tDCR2.placenta/' | bgzip >  ~{filebase}.vs.DCR2.placenta.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.skin.vs_zscore.uniq  | sed -e  's/$/\tDCR2.skin/' | bgzip >  ~{filebase}.vs.DCR2.skin.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/DCR/DCR2.stomach.vs_zscore.uniq  | sed -e  's/$/\tDCR2.stomach/' | bgzip >  ~{filebase}.vs.DCR2.stomach.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.CTCF-only,CTCF-bound.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.CTCF-only,CTCF-bound/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.CTCF-only,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.dELS,CTCF-bound.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.dELS,CTCF-bound/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.dELS,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.dELS.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.dELS/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.dELS.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.DNase-H3K4me3,CTCF-bound.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.DNase-H3K4me3,CTCF-bound/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.DNase-H3K4me3,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.DNase-H3K4me3.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.DNase-H3K4me3/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.DNase-H3K4me3.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.pELS,CTCF-bound.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.pELS,CTCF-bound/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.pELS,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.pELS.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.pELS/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.pELS.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.PLS,CTCF-bound.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.PLS,CTCF-bound/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.PLS,CTCF-bound.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode3/Encode3_ccREs.hg38.PLS.vs_zscore.uniq  | sed -e  's/$/\tEncode3_ccREs.hg38.PLS/' | bgzip >  ~{filebase}.vs.Encode3_ccREs.hg38.PLS.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.Bivalent.vs_zscore.uniq  | sed -e  's/$/\tencode4.Bivalent/' | bgzip >  ~{filebase}.vs.encode4.Bivalent.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.ConstitutiveHet.vs_zscore.uniq  | sed -e  's/$/\tencode4.ConstitutiveHet/' | bgzip >  ~{filebase}.vs.encode4.ConstitutiveHet.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.CTCF.vs_zscore.uniq  | sed -e  's/$/\tencode4.CTCF/' | bgzip >  ~{filebase}.vs.encode4.CTCF.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.EnhancerLow.vs_zscore.uniq  | sed -e  's/$/\tencode4.EnhancerLow/' | bgzip >  ~{filebase}.vs.encode4.EnhancerLow.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.Enhancer.vs_zscore.uniq  | sed -e  's/$/\tencode4.Enhancer/' | bgzip >  ~{filebase}.vs.encode4.Enhancer.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.FacultativeHet.vs_zscore.uniq  | sed -e  's/$/\tencode4.FacultativeHet/' | bgzip >  ~{filebase}.vs.encode4.FacultativeHet.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.K9K36.vs_zscore.uniq  | sed -e  's/$/\tencode4.K9K36/' | bgzip >  ~{filebase}.vs.encode4.K9K36.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.PromoterFlanking.vs_zscore.uniq  | sed -e  's/$/\tencode4.PromoterFlanking/' | bgzip >  ~{filebase}.vs.encode4.PromoterFlanking.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.Promoter.vs_zscore.uniq  | sed -e  's/$/\tencode4.Promoter/' | bgzip >  ~{filebase}.vs.encode4.Promoter.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/encode4/encode4.Transcribed.vs_zscore.uniq  | sed -e  's/$/\tencode4.Transcribed/' | bgzip >  ~{filebase}.vs.encode4.Transcribed.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.all_cells.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.all_cells/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.all_cells.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.all_organs.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.all_organs/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.all_organs.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.enhancer_clusters.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.enhancer_clusters/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.enhancer_clusters.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.permissive_enhancers.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.permissive_enhancers/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.permissive_enhancers.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.robust_enhancers.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.robust_enhancers/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.robust_enhancers.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.specific_enhancers_cells.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.specific_enhancers_cells/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.specific_enhancers_cells.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.specific_enhancers_organs.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.specific_enhancers_organs/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.specific_enhancers_organs.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.ubiquitous_enhancers_cells.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.ubiquitous_enhancers_cells/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.ubiquitous_enhancers_cells.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom5_enhancers.hg38.ubiquitous_enhancers_organs.vs_zscore.uniq  | sed -e  's/$/\tfantom5_enhancers.hg38.ubiquitous_enhancers_organs/' | bgzip >  ~{filebase}.vs.fantom5_enhancers.hg38.ubiquitous_enhancers_organs.over_50perc_ovr.bed.gz
            bedtools intersect -wo -f .5 -a ~{SV_file}  -b noncoding/fantom5/fantom_enhancers.hg38.enhancers.siwei_paper.vs_zscore.uniq  | sed -e  's/$/\tfantom_enhancers.hg38.enhancers.siwei_paper/' | bgzip >  ~{filebase}.vs.fantom_enhancers.hg38.enhancers.siwei_paper.over_50perc_ovr.bed.gz



            zcat ~{filebase}.vs.ABC_A549_treated_with_ethanol_0_02_percent_for_1_hour_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_A549_treated_with_ethanol_0_02_percent_for_1_hour_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_A673_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_A673_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_adipose_tissue_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_adipose_tissue_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_adrenal_gland_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_adrenal_gland_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_adrenal_gland_fetal_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_adrenal_gland_fetal_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_astrocyte_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_astrocyte_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_B_cell_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_B_cell_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_bipolar_neuron_from_iPSC_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_bipolar_neuron_from_iPSC_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_BJAB_anti_IgM_anti_CD40_4hr_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_BJAB_anti_IgM_anti_CD40_4hr_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_BJAB_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_BJAB_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_body_of_pancreas_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_body_of_pancreas_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_breast_epithelium_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_breast_epithelium_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_brite_adipose_Loft2014.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_brite_adipose_Loft2014.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_cardiac_muscle_cell_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_cardiac_muscle_cell_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocytes_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocytes_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_1h_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_1h_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_4h_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_4h_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_d1_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_d1_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_d6_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_BG_d6_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_1h_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_1h_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_d1_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_d1_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_d6_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_LPS_d6_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_1h_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_1h_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_4h_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_4h_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_d1_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_d1_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_d6_Novakovic2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD14_positive_monocyte_treated_with_RPMI_d6_Novakovic2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD19_positive_B_cell_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD19_positive_B_cell_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD34_positive_mobilized_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD34_positive_mobilized_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD3_positive_T_cell_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD3_positive_T_cell_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD4_positive_helper_T_cell_Corces2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD4_positive_helper_T_cell_Corces2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD4_positive_helper_T_cell_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD4_positive_helper_T_cell_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD56_positive_natural_killer_cells_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD56_positive_natural_killer_cells_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD8_positive_alpha_beta_T_cell_Corces2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD8_positive_alpha_beta_T_cell_Corces2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_CD8_positive_alpha_beta_T_cell_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_CD8_positive_alpha_beta_T_cell_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_coronary_artery_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_coronary_artery_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_coronary_artery_smooth_muscle_cell_Miller2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_coronary_artery_smooth_muscle_cell_Miller2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_0_ng_mL_for_0_hour_Garber2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_0_ng_mL_for_0_hour_Garber2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_1_hour_Garber2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_1_hour_Garber2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_2_hour_Garber2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_2_hour_Garber2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_30_minute_Garber2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_30_minute_Garber2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_4_hour_Garber2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_4_hour_Garber2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_6_hour_Garber2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_dendritic_cell_treated_with_Lipopolysaccharide_100_ng_mL_for_6_hour_Garber2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_endothelial_cell_of_umbilical_vein_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_endothelial_cell_of_umbilical_vein_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_epithelial_cell_of_prostate_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_epithelial_cell_of_prostate_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_erythroblast_Corces2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_erythroblast_Corces2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_fibroblast_of_arm_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_fibroblast_of_arm_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_fibroblast_of_dermis_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_fibroblast_of_dermis_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_fibroblast_of_lung_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_fibroblast_of_lung_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_foreskin_fibroblast_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_foreskin_fibroblast_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_gastrocnemius_medialis_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_gastrocnemius_medialis_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_GM12878_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_GM12878_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_H1_BMP4_Derived_Mesendoderm_Cultured_Cells_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_H1_BMP4_Derived_Mesendoderm_Cultured_Cells_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_H1_BMP4_Derived_Trophoblast_Cultured_Cells_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_H1_BMP4_Derived_Trophoblast_Cultured_Cells_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_H1_Derived_Mesenchymal_Stem_Cells_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_H1_Derived_Mesenchymal_Stem_Cells_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_H1_Derived_Neuronal_Progenitor_Cultured_Cells_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_H1_Derived_Neuronal_Progenitor_Cultured_Cells_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_H1_hESC_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_H1_hESC_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_H7.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_H7.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_H9_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_H9_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_HAP1.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_HAP1.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_HCT116_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_HCT116_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_heart_ventricle_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_heart_ventricle_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_HeLa_S3_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_HeLa_S3_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_hepatocyte_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_hepatocyte_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_HepG2_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_HepG2_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_HT29.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_HT29.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_IMR90_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_IMR90_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_induced_pluripotent_stem_cell_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_induced_pluripotent_stem_cell_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_iPS_DF_19_11_Cell_Line_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_iPS_DF_19_11_Cell_Line_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_Jurkat_anti_CD3_PMA_4hr_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_Jurkat_anti_CD3_PMA_4hr_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_Jurkat_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_Jurkat_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_K562_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_K562_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_Karpas_422_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_Karpas_422_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_keratinocyte_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_keratinocyte_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_large_intestine_fetal_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_large_intestine_fetal_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_liver_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_liver_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_LNCAP.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_LNCAP.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.loeuf.0.15.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.abc.loeuf.0.15.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_LoVo.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_LoVo.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_mammary_epithelial_cell_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_mammary_epithelial_cell_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_MCF10A_Ji2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_MCF10A_Ji2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_MCF10A_treated_with_TAM24hr_Ji2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_MCF10A_treated_with_TAM24hr_Ji2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_MCF_7_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_MCF_7_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_MDA_MB_231.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_MDA_MB_231.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_megakaryocyte_erythroid_progenitor_Corces2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_megakaryocyte_erythroid_progenitor_Corces2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.merge.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.abc.merge.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_MM_1S_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_MM_1S_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_muscle_of_leg_fetal_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_muscle_of_leg_fetal_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_muscle_of_trunk_fetal_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_muscle_of_trunk_fetal_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_myotube_originated_from_skeletal_muscle_myoblast_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_myotube_originated_from_skeletal_muscle_myoblast_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_natural_killer_cell_Corces2016.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_natural_killer_cell_Corces2016.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_NCCIT.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_NCCIT.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_OCI_LY7_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_OCI_LY7_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_osteoblast_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_osteoblast_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_ovary_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_ovary_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_Panc1_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_Panc1_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_pancreas_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_pancreas_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_PC_9_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_PC_9_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.pHaplo.0.9.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.abc.pHaplo.0.9.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_placenta_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_placenta_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_psoas_muscle_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_psoas_muscle_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.abc.pTriplo.0.84.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.abc.pTriplo.0.84.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_sigmoid_colon_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_sigmoid_colon_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_skeletal_muscle_myoblast_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_skeletal_muscle_myoblast_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_SK_N_SH_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_SK_N_SH_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_small_intestine_fetal_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_small_intestine_fetal_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_spinal_cord_fetal_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_spinal_cord_fetal_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_spleen_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_spleen_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_stomach_fetal_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_stomach_fetal_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_stomach_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_stomach_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_T_cell_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_T_cell_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP1_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP1_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP1_LPS_4hr_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP1_LPS_4hr_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_1_macrophage_VanBortle2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_1_macrophage_VanBortle2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_1_monocyte_VanBortle2017.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_1_monocyte_VanBortle2017.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_0h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_0h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_120h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_120h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_12h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_12h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_1h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_1h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_24h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_24h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_2h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_2h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_48h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_48h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_6h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_6h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_72h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_72h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_96h.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_THP_pmaLPS_ATAC_96h.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_thymus_fetal_Roadmap.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_thymus_fetal_Roadmap.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_thyroid_gland_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_thyroid_gland_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_transverse_colon_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_transverse_colon_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_trophoblast_cell_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_trophoblast_cell_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_U937_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_U937_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_U937_LPS_4hr_Engreitz.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_U937_LPS_4hr_Engreitz.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_uterus_ENCODE.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_uterus_ENCODE.over_50perc_ovr.stat
            zcat ~{filebase}.vs.ABC_white_adipose_Loft2014.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.ABC_white_adipose_Loft2014.over_50perc_ovr.stat
            zcat ~{filebase}.vs.atac-seq_OCRs_2020.hg38.autosomes.reformat.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.atac-seq_OCRs_2020.hg38.autosomes.reformat.over_50perc_ovr.stat
            zcat ~{filebase}.vs.atac-seq_pREs_enhancers_2020.hg38.autosomes.reformat.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.atac-seq_pREs_enhancers_2020.hg38.autosomes.reformat.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.Enh.txEnh.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.chromHMM.Enh.txEnh.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.insulator.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.chromHMM.insulator.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.polycomb.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.chromHMM.polycomb.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.chromHMM.promoter.TSS.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.chromHMM.promoter.TSS.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.adrenal_glands.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.adrenal_glands.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.brain.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.brain.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.eye.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.eye.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.face.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.face.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.general.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.general.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.heart.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.heart.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.intestine.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.intestine.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.kidney.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.kidney.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.limb.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.limb.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.liver.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.liver.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.lung.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.lung.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.muscle.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.muscle.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.neural_tube.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.neural_tube.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.placenta.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.placenta.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.skin.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.skin.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DCR2.stomach.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DCR2.stomach.over_50perc_ovr.stat
            zcat ~{filebase}.vs.DNaseHS.UCSC.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.DNaseHS.UCSC.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.CTCF-only,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.CTCF-only,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.dELS,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.dELS,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.dELS.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.dELS.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.DNase-H3K4me3,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.DNase-H3K4me3,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.DNase-H3K4me3.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.DNase-H3K4me3.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.pELS,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.pELS,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.pELS.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.pELS.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.PLS,CTCF-bound.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.PLS,CTCF-bound.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Encode3_ccREs.hg38.PLS.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Encode3_ccREs.hg38.PLS.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Bivalent.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.Bivalent.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.ConstitutiveHet.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.ConstitutiveHet.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.CTCF.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.CTCF.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.EnhancerLow.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.EnhancerLow.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Enhancer.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.Enhancer.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.FacultativeHet.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.FacultativeHet.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.K9K36.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.K9K36.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.PromoterFlanking.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.PromoterFlanking.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Promoter.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.Promoter.over_50perc_ovr.stat
            zcat ~{filebase}.vs.encode4.Transcribed.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.encode4.Transcribed.over_50perc_ovr.stat
            zcat ~{filebase}.vs.enhancer_DCR2v3.hg38.autosomes.reformat.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.enhancer_DCR2v3.hg38.autosomes.reformat.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.all_cells.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.all_cells.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.all_organs.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.all_organs.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.enhancer_clusters.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.enhancer_clusters.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.permissive_enhancers.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.permissive_enhancers.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.robust_enhancers.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.robust_enhancers.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.specific_enhancers_cells.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.specific_enhancers_cells.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.specific_enhancers_organs.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.specific_enhancers_organs.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.ubiquitous_enhancers_cells.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.ubiquitous_enhancers_cells.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom5_enhancers.hg38.ubiquitous_enhancers_organs.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom5_enhancers.hg38.ubiquitous_enhancers_organs.over_50perc_ovr.stat
            zcat ~{filebase}.vs.fantom_enhancers.hg38.enhancers.siwei_paper.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.fantom_enhancers.hg38.enhancers.siwei_paper.over_50perc_ovr.stat
            zcat ~{filebase}.vs.GW18_cortex_enh.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.GW18_cortex_enh.over_50perc_ovr.stat
            zcat ~{filebase}.vs.GW18_cortex_sil.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.GW18_cortex_sil.over_50perc_ovr.stat
            zcat ~{filebase}.vs.hg38_gene_deserts.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.hg38_gene_deserts.over_50perc_ovr.stat
            zcat ~{filebase}.vs.human_fetal_brain_TF_footprints.hg38.autosomes.reformat.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.human_fetal_brain_TF_footprints.hg38.autosomes.reformat.over_50perc_ovr.stat
            zcat ~{filebase}.vs.masked.regions.merged.LCRs.gaps.hg38.autosomes.reformat.sorted.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.masked.regions.merged.LCRs.gaps.hg38.autosomes.reformat.sorted.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Neuron_PROcap_enhancers.hg38.autosomes.reformat.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Neuron_PROcap_enhancers.hg38.autosomes.reformat.over_50perc_ovr.stat
            zcat ~{filebase}.vs.Shi_Lab.TAD_Boundry.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.Shi_Lab.TAD_Boundry.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.TFChip.UCSC.hg38.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.TFChip.UCSC.hg38.over_50perc_ovr.stat
            zcat ~{filebase}.vs.tri1_brain_atac.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.tri1_brain_atac.over_50perc_ovr.stat
            zcat ~{filebase}.vs.vista_cores.neg.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.vista_cores.neg.over_50perc_ovr.stat
            zcat ~{filebase}.vs.vista_cores.pos.over_50perc_ovr.bed.gz | cut -f4 | sort | uniq -c >  ~{filebase}.vs.vista_cores.pos.over_50perc_ovr.stat


            ls *.over_50perc_ovr.stat > SV_vs_nc_stat_list.tsv
            Rscript ./src/integrate_stats_across_nc_elements.R -i SV_vs_nc_stat_list.tsv -o ~{filebase}.vs.nc_elements.integrated -p ~{filebase}
            bgzip ~{filebase}.vs.nc_elements.integrated

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



