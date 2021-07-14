version 1.0

# Author: Xuefang Zhao <XZHAO12@mgh.harvard.edu>

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks
import "AnnotateExternalAFperContig.wdl" as AnnotateExternalAFperContig

workflow AnnotateExternalAF {
    input {
        File vcf
        File vcf_idx
        File ref_bed
        Array[String] population
        Array[String] contigs
        String ref_prefix
        String prefix

        Int max_shards_per_chrom_step1
        Int min_records_per_shard_step1

        String sv_base_mini_docker
        String sv_pipeline_docker

        # overrides for local tasks
        RuntimeAttr? runtime_attr_modify_vcf
        RuntimeAttr? runtime_override_combine_vcfs
        RuntimeAttr? runtime_override_split_vcf

    }
    call SplitBed as split_ref_bed{
        input:
            bed = ref_bed,
            sv_base_mini_docker = sv_base_mini_docker
    }
    call SplitVcf as split_query_vcf{
        input:
            vcf = vcf,
            sv_pipeline_docker = sv_pipeline_docker
    }

    Array[String] svtype_list = ["DEL","DUP","INS","INV_CPX","BND_CTX"]

    scatter ( contig in contigs ) {
        call AnnotateExternalAFperContig.AnnotateExternalAFperContig as AnnotateExternalAFperContig{
            input:
                vcf = vcf,
                vcf_idx = vcf_idx,
                ref_bed = ref_bed,
                split_query_vcf_del = split_query_vcf.del,
                split_query_vcf_dup = split_query_vcf.dup,
                split_query_vcf_ins = split_query_vcf.ins,
                split_query_vcf_inv = split_query_vcf.inv,
                split_query_vcf_bnd = split_query_vcf.bnd,
                split_ref_bed_del = split_ref_bed.del,
                split_ref_bed_dup = split_ref_bed.dup,
                split_ref_bed_ins = split_ref_bed.ins,
                split_ref_bed_inv = split_ref_bed.inv,
                split_ref_bed_bnd = split_ref_bed.bnd,
                population = population,
                contig = contig,
                ref_prefix = ref_prefix,

                max_shards_per_chrom_step1 = max_shards_per_chrom_step1,
                min_records_per_shard_step1 = min_records_per_shard_step1,
                sv_base_mini_docker = sv_base_mini_docker,
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_override_split_vcf = runtime_override_split_vcf,
                runtime_attr_modify_vcf = runtime_attr_modify_vcf
        }
    }

    call MiniTasks.ConcatVcfs as CombineVcfStep2 {
      input:
        vcfs = AnnotateExternalAFperContig.annotated_vcf,
        vcfs_idx = AnnotateExternalAFperContig.annotated_vcf_tbi,
        naive = true,
        outfile_prefix = "~{prefix}.annotated.vcf",
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_override_combine_vcfs
    }

     output{
        File annotated_vcf = CombineVcfStep2.concat_vcf
        File annotated_vcf_tbi = CombineVcfStep2.concat_vcf_idx
    }

}

task SplitBed{
    input{
        File bed
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

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_base_mini_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(bed, ".bed.gz")
    
    command <<<
        zcat ~{bed} | head -1 > header
        cat header <(zcat ~{bed} | awk '{if ($6=="DEL") print}') > ~{prefix}.DEL.bed
        cat header <(zcat ~{bed} | awk '{if ($6=="DUP") print}') > ~{prefix}.DUP.bed
        cat header <(zcat ~{bed} | awk '{if ($6=="INS" || $6=="INS:ME" || $6=="INS:ME:ALU" || $6=="INS:ME:LINE1" || $6=="INS:ME:SVA" || $6=="ALU" || $6=="LINE1" || $6=="SVA" || $6=="HERVK" ) print}') > ~{prefix}.INS.bed
        cat header <(zcat ~{bed} | awk '{if ($6=="INV" || $6=="CPX") print}' ) > ~{prefix}.INV_CPX.bed
        cat header <(zcat ~{bed} | awk '{if ($6=="BND" || $6=="CTX") print}' ) > ~{prefix}.BND_CTX.bed
    >>>

    output{
        File del = "~{prefix}.DEL.bed"
        File dup = "~{prefix}.DUP.bed"
        File ins = "~{prefix}.INS.bed"
        File inv = "~{prefix}.INV_CPX.bed"
        File bnd = "~{prefix}.BND_CTX.bed"
    }    
}

task SplitVcf{
    input{
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

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf, ".vcf.gz")
    
    command <<<
        svtk vcf2bed -i SVTYPE -i SVLEN ~{vcf} - |
            cut -f1-4,7-8 > tmp.bed
        head -1 tmp.bed > ~{prefix}.bed
        awk 'NR > 1' < tmp.bed \
            | sort -k1,1V -k2,2n -k3,3n >> ~{prefix}.bed
        rm tmp.bed
        head -1 ~{prefix}.bed > header
        cat header <(awk '{if ($5=="DEL") print}' ~{prefix}.bed )> ~{prefix}.DEL.bed
        cat header <(awk '{if ($5=="DUP") print}' ~{prefix}.bed )> ~{prefix}.DUP.bed
        cat header <(awk '{if ($5=="INS" || $5=="INS:ME" || $5=="INS:ME:ALU" || $5=="INS:ME:LINE1" || $5=="INS:ME:SVA" || $5=="ALU" || $5=="LINE1" || $5=="SVA" || $5=="HERVK" ) print}' ~{prefix}.bed )> ~{prefix}.INS.bed
        cat header <(awk '{if ($5=="INV" || $5=="CPX") print}' ~{prefix}.bed )> ~{prefix}.INV_CPX.bed
        cat header <(awk '{if ($5=="BND" || $5=="CTX") print}' ~{prefix}.bed )> ~{prefix}.BND_CTX.bed
    >>>

    output{
        File bed = "~{prefix}.bed"
        File del = "~{prefix}.DEL.bed"
        File dup = "~{prefix}.DUP.bed"
        File ins = "~{prefix}.INS.bed"
        File inv = "~{prefix}.INV_CPX.bed"
        File bnd = "~{prefix}.BND_CTX.bed"
    }
}

task BedtoolsClosest{
    input{
        File bed_a
        File bed_b
        String svtype
        String contig
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

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }
    
    command <<<
        awk '{if ($1=="~{contig}") print}' ~{bed_a} > filea.bed
        awk '{if ($1=="~{contig}") print}' ~{bed_b} > fileb.bed

        paste <(head -1 ~{bed_a}) <(head -1 ~{bed_b}) | sed -e "s/#//g" > ~{svtype}.bed

        bedtools closest -wo -a <(sort -k1,1 -k2,2n filea.bed) -b <(sort -k1,1 -k2,2n fileb.bed) >> ~{svtype}.bed
    >>>

    output{
        File output_bed = "~{svtype}.bed"
    }
}

task SelectMatchedSVs{
    input{
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

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(input_bed, ".bed")
    File pop_list = write_lines(population)

    command <<<
        Rscript /opt/sv-pipeline/05_annotation/scripts/R1.bedtools_closest_CNV.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison \
            -p ~{pop_list}
    >>>

    output{
        File output_comp = "~{prefix}.comparison"
    }    
}

task SelectMatchedINSs{
    input{
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

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(input_bed, ".bed")
    File pop_list = write_lines(population)

    command <<<
        Rscript /opt/sv-pipeline/05_annotation/scripts/R2.bedtools_closest_INS.R \
            -i ~{input_bed} \
            -o ~{prefix}.comparison \
            -p ~{pop_list}
    >>>

    output{
        File output_comp = "~{prefix}.comparison"
    }
}

task ModifyVcf{
    input{
        Array[File] labeled_del
        Array[File] labeled_dup
        Array[File] labeled_ins
        Array[File] labeled_inv
        Array[File] labeled_bnd
        File vcf
        String ref_prefix
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

    RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
    
    runtime {
        memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
        disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
        cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
        preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
        maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
        docker: sv_pipeline_docker
        bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
    }

    String prefix = basename(vcf,'.vcf.gz')
    command <<<
        cat ~{sep=" " labeled_del} > labeled.bed
        cat ~{sep=" " labeled_dup} >> labeled.bed
        cat ~{sep=" " labeled_ins} >> labeled.bed
        cat ~{sep=" " labeled_inv} >> labeled.bed
        cat ~{sep=" " labeled_bnd} >> labeled.bed

        python <<CODE
        import os
        fin=os.popen(r'''zcat %s'''%("~{vcf}"))
        header = []
        body = {}
        SVID_key = []
        for line in fin:
            pin=line.strip().split()
            if pin[0][:2]=='##':
                header.append(pin)
            else:
                body[pin[2]]=pin
                SVID_key.append(pin[2])
        header.append(['##INFO=<ID='+"~{ref_prefix}"+'_SVID'+',Number=1,Type=Float,Description="Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in gnomad.">'])

        fin.close()
        fin=open('labeled.bed')
        colname = fin.readline().strip().split()

        for j in range(len(colname)-1):
            if j>1:
                header.append(['##INFO=<ID='+"~{ref_prefix}"+'_'+colname[j]+',Number=1,Type=Float,Description="Allele frequency (for biallelic sites) or copy-state frequency (for multiallelic sites) of an overlapping event in gnomad.">'])

        for line in fin:
            pin=line.strip().split()
            if pin[0]=='query_svid': continue
            info_add = ["~{ref_prefix}"+'_SVID'+'='+pin[1]]
            for j in range(len(colname)-1):
                if j>1:
                    info_add.append("~{ref_prefix}"+'_'+colname[j]+'='+pin[j])
            body[pin[0]][7]+=';'+';'.join(info_add)
        fin.close()

        fo=open('~{prefix}.annotated.vcf','w')
        for i in header:
            print(' '.join(i), file=fo)
        for i in SVID_key:
            print('\t'.join(body[i]), file=fo)
        fo.close()
        CODE

        bgzip ~{prefix}.annotated.vcf
        tabix ~{prefix}.annotated.vcf.gz
    >>>

    output{
        File annotated_vcf = "~{prefix}.annotated.vcf.gz"
        File annotated_vcf_tbi = "~{prefix}.annotated.vcf.gz.tbi"
    }        
}




