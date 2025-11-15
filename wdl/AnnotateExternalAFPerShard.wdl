version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks

workflow AnnotateExternalAFPerShard {
    input {
        File vcf
        File vcf_idx
        String prefix
        File split_ref_bed_del
        File split_ref_bed_dup
        File split_ref_bed_ins
        File split_ref_bed_inv
        File split_ref_bed_bnd

        Array[String] population
        String ref_prefix

        String sv_base_mini_docker
        String sv_pipeline_docker

        # overrides for local tasks
        RuntimeAttr? runtime_attr_modify_vcf
        RuntimeAttr? runtime_attr_split_query_vcf
        RuntimeAttr? runtime_attr_bedtools_closest
        RuntimeAttr? runtime_attr_select_matched_svs
    }

    call SplitQueryVcf {
        input:
            vcf = vcf,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_split_query_vcf
    }

    call BedtoolsClosest as compare_del {
        input:
            bed_a = SplitQueryVcf.del,
            bed_b = split_ref_bed_del,
            svtype = "del",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    call BedtoolsClosest as compare_dup {
        input:
            bed_a = SplitQueryVcf.dup,
            bed_b = split_ref_bed_dup,
            svtype = "dup",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    call BedtoolsClosest as compare_ins {
        input:
            bed_a = SplitQueryVcf.ins,
            bed_b = split_ref_bed_ins,
            svtype = "ins",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    call BedtoolsClosest as compare_inv {
        input:
            bed_a = SplitQueryVcf.inv,
            bed_b = split_ref_bed_inv,
            svtype = "inv",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    call BedtoolsClosest as compare_bnd {
        input:
            bed_a = SplitQueryVcf.bnd,
            bed_b = split_ref_bed_bnd,
            svtype = "bnd",
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_bedtools_closest
    }

    call SelectMatchedSVs as calcu_del {
        input:
            svtype = "del",
            input_bed = compare_del.output_bed,
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call SelectMatchedSVs as calcu_dup {
        input:
            input_bed=compare_dup.output_bed,
            svtype="dup",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call SelectMatchedINSs as calcu_ins {
        input:
            input_bed = compare_ins.output_bed,
            svtype = "ins",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call SelectMatchedSVs as calcu_inv {
        input:
            input_bed = compare_inv.output_bed,
            svtype = "inv",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }

    call SelectMatchedINSs as calcu_bnd {
        input:
            input_bed = compare_bnd.output_bed,
            svtype = "bnd",
            population = population,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_select_matched_svs
    }
 
    call ModifyVcf {
        input:
            labeled_del = calcu_del.output_comp,
            labeled_dup = calcu_dup.output_comp,
            labeled_ins = calcu_ins.output_comp,
            labeled_inv = calcu_inv.output_comp,
            labeled_bnd = calcu_bnd.output_comp,
            vcf = vcf,
            prefix = prefix,
            ref_prefix = ref_prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_modify_vcf       
    }

    output {
        File annotated_vcf = ModifyVcf.annotated_vcf
        File annotated_vcf_tbi = ModifyVcf.annotated_vcf_tbi
    }

}

task SplitRefBed {
    input {
        File bed
        String contig
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 3,
        disk_gb: 5,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(bed, ".bed.gz")
    
    command <<<
        set -eu
        zcat ~{bed} | head -1 > header
        set -o pipefail
        cat header <(zcat ~{bed} | awk '{if ($1=="~{contig}" && $6=="DEL") print}') > ~{prefix}.~{contig}.DEL.bed
        cat header <(zcat ~{bed} | awk '{if ($1=="~{contig}" && $6=="DUP") print}') > ~{prefix}.~{contig}.DUP.bed
        cat header <(zcat ~{bed} | awk '{if ($1=="~{contig}" && $6=="INS" || $6=="INS:ME" || $6=="INS:ME:ALU" || $6=="INS:ME:LINE1" || $6=="INS:ME:SVA" || $6=="ALU" || $6=="LINE1" || $6=="SVA" || $6=="HERVK" ) print}') > ~{prefix}.~{contig}.INS.bed
        cat header <(zcat ~{bed} | awk '{if ($1=="~{contig}" && $6=="INV" || $6=="CPX") print}' ) > ~{prefix}.~{contig}.INV_CPX.bed
        cat header <(zcat ~{bed} | awk '{if ($1=="~{contig}" && $6=="BND" || $6=="CTX") print}' ) > ~{prefix}.~{contig}.BND_CTX.bed
    >>>

    output {
        File del = "~{prefix}.~{contig}.DEL.bed"
        File dup = "~{prefix}.~{contig}.DUP.bed"
        File ins = "~{prefix}.~{contig}.INS.bed"
        File inv = "~{prefix}.~{contig}.INV_CPX.bed"
        File bnd = "~{prefix}.~{contig}.BND_CTX.bed"
    }    
}

task SplitQueryVcf {
    input {
        File vcf
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 3,
        disk_gb: 10,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, ".vcf.gz")
    
    command <<<
        set -euo pipefail
        svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} tmp.bed
        cut -f1-4,7-8 tmp.bed > ~{prefix}.bed
        set +o pipefail
        head -1 ~{prefix}.bed > header
        set -o pipefail
        cat header <(awk '{if ($5=="DEL") print}' ~{prefix}.bed )> ~{prefix}.DEL.bed
        cat header <(awk '{if ($5=="DUP") print}' ~{prefix}.bed )> ~{prefix}.DUP.bed
        cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK" ) print}' ~{prefix}.bed )> ~{prefix}.INS.bed
        cat header <(awk '{if ($5=="INV" || $5=="CPX") print}' ~{prefix}.bed )> ~{prefix}.INV_CPX.bed
        cat header <(awk '{if ($5=="BND" || $5=="CTX") print}' ~{prefix}.bed )> ~{prefix}.BND_CTX.bed
    >>>

    output {
        File bed = "~{prefix}.bed"
        File del = "~{prefix}.DEL.bed"
        File dup = "~{prefix}.DUP.bed"
        File ins = "~{prefix}.INS.bed"
        File inv = "~{prefix}.INV_CPX.bed"
        File bnd = "~{prefix}.BND_CTX.bed"
    }
}

task BedtoolsClosest {
    input {
        File bed_a
        File bed_b
        String svtype
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 3,
        disk_gb: 5,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        set -eu
        paste <(head -1 ~{bed_a}) <(head -1 ~{bed_b}) | sed -e "s/#//g" > ~{svtype}.bed
        set -o pipefail
        bedtools closest -wo -a <(sort -k1,1 -k2,2n ~{bed_a}) -b <(sort -k1,1 -k2,2n ~{bed_b}) >> ~{svtype}.bed
    >>>

    output {
        File output_bed = "~{svtype}.bed"
    }
}

task SelectMatchedSVs {
    input {
        File input_bed
        String svtype
        Array[String] population
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 3,
        disk_gb: 5,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(input_bed, ".bed")
    File pop_list = write_lines(population)

    command <<<
        set -euo pipefail
        Rscript /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison \
            -p ~{pop_list}
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }    
}

task SelectMatchedINSs {
    input {
        File input_bed
        String svtype
        Array[String] population
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 3,
        disk_gb: 5,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(input_bed, ".bed")
    File pop_list = write_lines(population)

    command <<<
        Rscript /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison \
            -p ~{pop_list}
    >>>

    output {
        File output_comp = "~{prefix}.comparison"
    }
}

task ModifyVcf {
    input {
        File labeled_del
        File labeled_dup
        File labeled_ins
        File labeled_inv
        File labeled_bnd
        File vcf
        String prefix
        String ref_prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr runtime_default = object {
        mem_gb: 10,
        disk_gb: 15,
        cpu_cores: 1,
        preemptible_tries: 1,
        max_retries: 1,
        boot_disk_gb: 10
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_attr.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_attr.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    command <<<
        cat ~{labeled_del} > labeled.bed
        cat ~{labeled_dup} >> labeled.bed
        cat ~{labeled_ins} >> labeled.bed
        cat ~{labeled_inv} >> labeled.bed
        cat ~{labeled_bnd} >> labeled.bed

        python /opt/sv-pipeline/05_annotation/scripts/annotate_external_af_modify_vcf.py \
            --vcf ~{vcf} \
            --ref-prefix ~{ref_prefix} \
            --labeled-bed labeled.bed \
            --output-filename ~{prefix}.annotated.vcf

        bgzip ~{prefix}.annotated.vcf
        tabix ~{prefix}.annotated.vcf.gz
    >>>

    output {
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_tbi = "~{prefix}.annotated.vcf.gz.tbi"
    }        
}




