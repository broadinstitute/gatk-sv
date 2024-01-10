version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "TasksBenchmark.wdl" as tasks10

import "Duphold.wdl" as duphold
import "RdPeSrAnno.wdl" as rdpesr

workflow AnnoILFeatures {
    input{
        String prefix
        String il_bam
        String il_bam_bai
        File vcf_file

        File pe_matrix
        File pe_index
        File sr_matrix
        File sr_index
        File rd_matrix
        File rd_index

        File ref_SegDup
        File ref_SimpRep
        File ref_RepMask

        File ref_fasta
        File ref_fai
        File ref_dict
        File contig_list

        Array[File] raw_vcfs
        Array[String] raw_algorithms

        String rdpesr_benchmark_docker
        String vapor_docker
        String duphold_docker
        String sv_base_mini_docker
        String sv_pipeline_docker

        RuntimeAttr? runtime_attr_Vapor 
        RuntimeAttr? runtime_attr_duphold
        RuntimeAttr? runtime_attr_rdpesr
        RuntimeAttr? runtime_attr_bcf2vcf
        RuntimeAttr? runtime_attr_LocalizeCram
        RuntimeAttr? runtime_attr_vcf2bed
        RuntimeAttr? runtime_attr_SplitVcf
        RuntimeAttr? runtime_attr_ConcatBeds
        RuntimeAttr? runtime_attr_ConcatVcfs
    }

    call vcf2bed{
        input:
            vcf = vcf_file,
            prefix = prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_vcf2bed
    }

    call RunGenomicContextAnnotation{
        input:
            bed = vcf2bed.bed,
            prefix = prefix,
            ref_SegDup = ref_SegDup,
            ref_SimpRep = ref_SimpRep,
            ref_RepMask = ref_RepMask,
            rdpesr_benchmark_docker = rdpesr_benchmark_docker,
            runtime_attr_override = runtime_attr_rdpesr
    }

    call Bed2QueryAndRef{
        input:
            bed = vcf2bed.bed,
            sv_base_mini_docker = sv_base_mini_docker
    }

    call ExtracGTGQ{
      input:
        prefix = prefix,
        vcf_file = vcf_file,
        sv_pipeline_docker = sv_pipeline_docker
    }

    scatter (i in range(length(raw_vcfs))){

        call vcf2bed as vcf2bed_raw{
            input:
                vcf = raw_vcfs[i],
                prefix = "${prefix}.${raw_algorithms[i]}",
                sv_pipeline_docker = sv_pipeline_docker,
                runtime_attr_override = runtime_attr_vcf2bed
        }


        call Bed2QueryAndRef as Bed2QueryAndRef_Raw{
            input:
                bed = vcf2bed_raw.bed,
                sv_base_mini_docker = sv_base_mini_docker
        }

        call BedComparison as BedComparison_vs_raw{
            input:
                query = Bed2QueryAndRef_Raw.query,
                ref = Bed2QueryAndRef.ref,
                prefix = "${prefix}.vs.${raw_algorithms[i]}",
                sv_pipeline_docker=sv_pipeline_docker
        }
    }

    call ExtracAlgorithmEvidenceFilter{
        input:
            prefix = prefix,
            vcf_file = vcf_file,
            sv_pipeline_docker=sv_pipeline_docker
    }

    call vcf2bed as vcf2bed_all{
        input:
            vcf = vcf_file,
            prefix = prefix,
            sv_pipeline_docker = sv_pipeline_docker,
            runtime_attr_override = runtime_attr_vcf2bed
    }

    call RunRdPeSrAnnotation{
        input:
            prefix = prefix,
            bed = vcf2bed_all.bed,
            pe_matrix = pe_matrix,
            pe_index = pe_index,
            sr_matrix = sr_matrix,
            sr_index = sr_index,
            rd_matrix = rd_matrix,
            rd_index = rd_index,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict=ref_dict,
            rdpesr_benchmark_docker = rdpesr_benchmark_docker,
            runtime_attr_override = runtime_attr_rdpesr
    }
    
    Array[String] contigs = transpose(read_tsv(contig_list))[0]
    scatter ( contig in contigs ) {
        call tasks10.SplitVcf as SplitVcf{
            input:
                contig = contig,
                vcf_file = vcf_file,
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_override=runtime_attr_SplitVcf
            }

        call tasks10.LocalizeCramRequestPay as LocalizeCramIL{
            input:
                contig = contig,
                ref_fasta=ref_fasta,
                ref_fai=ref_fai,
                ref_dict=ref_dict,
                project_id="talkowski-sv-gnomad",
                bam_or_cram_file=il_bam,
                bam_or_cram_index=il_bam_bai,
                sv_pipeline_docker=sv_pipeline_docker,
                runtime_attr_override=runtime_attr_LocalizeCram
            }

        call ShiftVcfForDuphold{
            input:
                prefix = contig,
                vcf_file = SplitVcf.contig_vcf,
                vcf_index = SplitVcf.contig_vcf_index,
                ref_fai = ref_fai,
                rdpesr_benchmark_docker = rdpesr_benchmark_docker,
                runtime_attr_override = runtime_attr_duphold
            } 

        call RunDupholdPerContig as RunDupholdPerContigIL{
            input:
                prefix = prefix,
                contig = contig,
                bam_or_cram_file=LocalizeCramIL.local_bam,
                bam_or_cram_index=LocalizeCramIL.local_bai,
                vcf_file = SplitVcf.contig_vcf,
                vcf_index = SplitVcf.contig_vcf_index,
                vcf_le_file = ShiftVcfForDuphold.le_flank,
                vcf_le_index = ShiftVcfForDuphold.le_flank_index,
                vcf_ri_file = ShiftVcfForDuphold.ri_flank,
                vcf_ri_index = ShiftVcfForDuphold.ri_flank_index,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                rdpesr_benchmark_docker = duphold_docker,
                runtime_attr_override = runtime_attr_duphold
            }

        call Bcf2Vcf as Bcf2VcfIL{
            input:
                prefix = prefix,
                contig = contig,
                bcf = RunDupholdPerContigIL.bcf,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_bcf2vcf
            }

        call Bcf2Vcf as Bcf2VcfIL_le_flank{
            input:
                prefix = prefix,
                contig = contig,
                bcf = RunDupholdPerContigIL.bcf_le,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_bcf2vcf
            }

        call Bcf2Vcf as Bcf2VcfIL_ri_flank{
            input:
                prefix = prefix,
                contig = contig,
                bcf = RunDupholdPerContigIL.bcf_ri,
                sv_base_mini_docker = sv_base_mini_docker,
                runtime_attr_override = runtime_attr_bcf2vcf
            }
   
    }

    call MiniTasks.ConcatVcfs as ConcatVcfsIL{
        input:
            vcfs=Bcf2VcfIL.vcf,
            outfile_prefix="~{prefix}.IL",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_ConcatVcfs
    }

    call MiniTasks.ConcatVcfs as ConcatVcfsIL_le_flank{
        input:
            vcfs=Bcf2VcfIL_le_flank.vcf,
            outfile_prefix="~{prefix}.IL_le_flank",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_ConcatVcfs
    }

    call MiniTasks.ConcatVcfs as ConcatVcfsIL_ri_flank{
        input:
            vcfs=Bcf2VcfIL_ri_flank.vcf,
            outfile_prefix="~{prefix}.IL_ri_flank",
            sv_base_mini_docker=sv_base_mini_docker,
            runtime_attr_override=runtime_attr_ConcatVcfs
    }

    output{
            File duphold_vcf_il = ConcatVcfsIL.concat_vcf
            File duphold_vcf_il_le = ConcatVcfsIL_le_flank.concat_vcf
            File duphold_vcf_il_ri = ConcatVcfsIL_ri_flank.concat_vcf

            File PesrAnno  = RunRdPeSrAnnotation.pesr_anno
            File RdAnno    = RunRdPeSrAnnotation.cov
            File RdAnno_le = RunRdPeSrAnnotation.cov_le_flank 
            File RdAnno_ri = RunRdPeSrAnnotation.cov_ri_flank

            File GCAnno = RunGenomicContextAnnotation.anno_bed
            File GTGQ = ExtracGTGQ.GQ_GT
            File vcf_info = ExtracAlgorithmEvidenceFilter.vcf_info
            Array[File] vs_raw = BedComparison_vs_raw.comparison
        }
    }
 
task RunDupholdPerContig{
    input{
        String prefix
        String contig
        File bam_or_cram_file
        File bam_or_cram_index
        File vcf_file
        File vcf_index
        File vcf_le_file
        File vcf_le_index
        File vcf_ri_file
        File vcf_ri_index
        File ref_fasta
        File ref_fai
        File ref_dict
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])


    output {
        File bcf = "~{prefix}.~{contig}.bcf"
        File bcf_le = "~{prefix}.~{contig}.le_flank.bcf"
        File bcf_ri = "~{prefix}.~{contig}.ri_flank.bcf"
    }
    command <<<

        set -Eeuo pipefail
        
        duphold -t 4 \
        -v ~{vcf_file} \
        -b ~{bam_or_cram_file} \
        -f ~{ref_fasta} \
        -o ~{prefix}.~{contig}.bcf

        duphold -t 4 \
        -v ~{vcf_le_file} \
        -b ~{bam_or_cram_file} \
        -f ~{ref_fasta} \
        -o ~{prefix}.~{contig}.le_flank.bcf

        duphold -t 4 \
        -v ~{vcf_ri_file} \
        -b ~{bam_or_cram_file} \
        -f ~{ref_fasta} \
        -o ~{prefix}.~{contig}.ri_flank.bcf

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

task vcf2bed{
    input{
        String prefix
        File vcf
        File? vcf_index
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 2,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    String filename = basename(vcf, ".vcf.gz")

    output {
        File bed = "~{prefix}.bed"
    }

    command <<<

        set -Eeuo pipefail
        
        gsutil cp ~{vcf} ./tmp.vcf.gz
        tabix -p vcf ./tmp.vcf.gz
        svtk vcf2bed -i SVTYPE -i SVLEN tmp.vcf.gz ~{prefix}.bed
        
    >>>
    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task Bed2QueryAndRef{
    input{
        File bed
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File query = "${filebase}.query.gz"
        File ref = "${filebase}.ref.gz"
    }

    String filebase=basename(bed,".bed")
    command <<<
        echo "#chroms tart end name SVTYPE SVLEN" | sed -e 's/ /\t/g' > ~{filebase}.query
        echo "#chrom start end VID svtype length AF samples" | sed -e 's/ /\t/g' > ~{filebase}.ref

        cut -f1-4,7,8 ~{bed} | grep -v "#" >> ~{filebase}.query
        cut -f1-4,7,8 ~{bed} | sed -e "s/$/\t0\t~{filebase}/" | grep -v "#" >> ~{filebase}.ref

        bgzip ~{filebase}.query
        bgzip ~{filebase}.ref

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

task BedComparison{
    input{
        File query
        File ref
        String prefix
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File comparison = "~{prefix}.bed"
    }

    command <<<
        bash /opt/sv-pipeline/scripts/vcf_qc/compare_callsets_V2.sh \
            -O ~{prefix}.bed -p ~{prefix} ~{query} ~{ref}
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task Bcf2Vcf{
    input{
        String prefix
        String contig
        File bcf
        String sv_base_mini_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vcf = "~{prefix}.~{contig}.duphold.vcf.gz"
    }

    command <<<
            set -Eeuo pipefail
            bcftools view ~{bcf} | bgzip > ~{prefix}.~{contig}.duphold.vcf.gz
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

task ShiftVcfForDuphold{
    input{
        String prefix
        File vcf_file
        File vcf_index
        File ref_fai
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 7.5, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File le_flank = "~{prefix}.le_flank.vcf.gz"
        File ri_flank = "~{prefix}.ri_flank.vcf.gz"
        File le_flank_index = "~{prefix}.le_flank.vcf.gz.tbi"
        File ri_flank_index = "~{prefix}.ri_flank.vcf.gz.tbi"
    }

    command <<<
        python3 /src/Modify_vcf_by_steps.py ~{vcf_file} ~{prefix}.ri_flank.vcf -s 1000 -c ~{ref_fai}
        python3 /src/Modify_vcf_by_steps.py ~{vcf_file} ~{prefix}.le_flank.vcf -s -1000 -c    ~{ref_fai}

        bgzip ~{prefix}.ri_flank.vcf
        bgzip ~{prefix}.le_flank.vcf

        tabix -p vcf ~{prefix}.ri_flank.vcf.gz
        tabix -p vcf ~{prefix}.le_flank.vcf.gz
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

task RunDuphold{
    input{
        String prefix
        File bam_or_cram_file
        File bam_or_cram_index
        File vcf_file
        File ref_fasta
        File ref_fai
        File ref_dict
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 2,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File bcf = "~{prefix}.bcf"
    }

    command <<<

        set -Eeuo pipefail
        
        duphold -t 4 \
        -v ~{vcf_file} \
        -b ~{bam_or_cram_file} \
        -f ~{ref_fasta} \
        -o ~{prefix}.bcf
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

task RunVaPoR{
    input{
        String prefix
        String contig
        File bam_or_cram_file
        File bam_or_cram_index
        File bed
        File ref_fasta
        File ref_fai
        File ref_dict
        String rdpesr_benchmark_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output {
        File vapor = "~{bed}.vapor"
        File vapor_plot = "~{prefix}.~{contig}.tar.gz"
    }

    command <<<

        set -Eeuo pipefail

        mkdir ~{prefix}.~{contig}
    
        vapor bed \
        --sv-input ~{bed} \
        --output-path ~{prefix}.~{contig} \
        --reference ~{ref_fasta} \
        --pacbio-input ~{bam_or_cram_index} \

        tar -czf ~{prefix}.~{contig}.tar.gz ~{prefix}.~{contig}
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

task RunRdPeSrAnnotation{
    input{
        String prefix
        File bed
        File pe_matrix
        File pe_index
        File sr_matrix
        File sr_index
        File rd_matrix
        File rd_index
        File ref_fasta
        File ref_fai
        File ref_dict
        String rdpesr_benchmark_docker
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

    output {
        File cov = "~{filebase}.bed.Rd.gz"
        File cov_ri_flank = "~{filebase}.ri_flank.Rd.gz"
        File cov_le_flank = "~{filebase}.le_flank.Rd.gz"
        File pesr_anno = "~{filebase}.bed.PeSr.gz"
    }

    String filebase = basename(bed,".bed")

    command <<<

        set -Eeuo pipefail

        Rscript /src/modify_bed_for_PE_SR_RD_labeling.R \
            -i ~{bed} \
            --le_bp ~{bed}.le_bp \
            --ri_bp ~{bed}.ri_bp \
            --le_flank ~{bed}.le_flank \
            --ri_flank ~{bed}.ri_flank

        zcat ~{rd_matrix} | grep -v '@' | grep -v CONTIG |bgzip >    bincov.tsv.gz
        Rscript /src/bincov_to_normCov.R -i bincov.tsv.gz
        bgzip normCov.tsv
        tabix -b 2 -e 2 normCov.tsv.gz

        python3 /src/add_RD_to_SVs.py ~{bed} normCov.tsv.gz ~{filebase}.bed.Rd
        python3 /src/add_RD_to_SVs.py ~{bed}.ri_flank normCov.tsv.gz ~{filebase}.ri_flank.Rd
        python3 /src/add_RD_to_SVs.py ~{bed}.le_flank normCov.tsv.gz ~{filebase}.le_flank.Rd
        python3 /src/add_SR_PE_to_PB_INS.V2.py ~{bed} ~{pe_matrix} ~{sr_matrix} ~{filebase}.bed.PeSr

        bgzip ~{filebase}.bed.Rd
        bgzip ~{filebase}.ri_flank.Rd
        bgzip ~{filebase}.le_flank.Rd
        bgzip ~{filebase}.bed.PeSr

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

task RunGenomicContextAnnotation{
    input{
        File bed
        File ref_SegDup
        File ref_SimpRep
        File ref_RepMask
        String prefix
        String rdpesr_benchmark_docker
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

    command <<<

        awk '{print $1,$2,$2,$4,$5}' ~{bed} | sed -e 's/ /\t/g' >    ~{prefix}.le_bp
        awk '{print $1,$3,$3,$4,$5}' ~{bed} | sed -e 's/ /\t/g' >    ~{prefix}.ri_bp
        bedtools coverage -a ~{prefix}.le_bp -b ~{ref_RepMask} | awk '{if ($9>0) print}'>    ~{prefix}.le_bp.vs.RM
        bedtools coverage -a ~{prefix}.le_bp -b ~{ref_SegDup}    | awk '{if ($9>0) print}'>    ~{prefix}.le_bp.vs.SD
        bedtools coverage -a ~{prefix}.le_bp -b ~{ref_SimpRep} | awk '{if ($9>0) print}'>    ~{prefix}.le_bp.vs.SR
        bedtools coverage -a ~{prefix}.ri_bp -b ~{ref_RepMask} | awk '{if ($9>0) print}'>    ~{prefix}.ri_bp.vs.RM
        bedtools coverage -a ~{prefix}.ri_bp -b ~{ref_SegDup}    | awk '{if ($9>0) print}'>    ~{prefix}.ri_bp.vs.SD
        bedtools coverage -a ~{prefix}.ri_bp -b ~{ref_SimpRep} | awk '{if ($9>0) print}'>    ~{prefix}.ri_bp.vs.SR


        Rscript /src/add_GC_anno_to_bed.R \
        -b ~{bed} \
        -o ~{prefix}.GC_anno.bed \
        --left_vs_SR    ~{prefix}.le_bp.vs.SR \
        --left_vs_SD    ~{prefix}.le_bp.vs.SD \
        --left_vs_RM    ~{prefix}.le_bp.vs.RM \
        --right_vs_SR ~{prefix}.ri_bp.vs.SR \
        --right_vs_SD ~{prefix}.ri_bp.vs.SD \
        --right_vs_RM ~{prefix}.ri_bp.vs.RM 
    >>>

    output{
        File anno_bed = "~{prefix}.GC_anno.bed"
    }

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

task ExtracGTGQ{
    input{
        String prefix
        File vcf_file
        String sv_pipeline_docker
        RuntimeAttr? runtime_attr_override
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
  
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File GQ_GT = "~{prefix}.SVID_gt.tsv"
    }

    command <<<
        zcat ~{vcf_file} | grep -v '#' > ~{prefix}.SVID_gt

        python <<CODE
        import os
        fin=open("~{prefix}.SVID_gt")
        svid_gt={}
        for line in fin:
          pin=line.strip().split()
          svid_gt[pin[2]]=[pin[9].split(':')[pin[8].split(':').index('GT')], 
                          pin[9].split(':')[pin[8].split(':').index('GQ')], 
                          pin[9].split(':')[pin[8].split(':').index('RD_CN')], 
                          pin[9].split(':')[pin[8].split(':').index('RD_GQ')], 
                          pin[9].split(':')[pin[8].split(':').index('PE_GT')], 
                          pin[9].split(':')[pin[8].split(':').index('PE_GQ')], 
                          pin[9].split(':')[pin[8].split(':').index('SR_GT')], 
                          pin[9].split(':')[pin[8].split(':').index('SR_GQ')]]
        fin.close()

        fo=open("~{prefix}.SVID_gt.tsv", 'w')
        print('\t'.join(['SVID','GT','GQ','RD_CN','RD_GQ','PE_GT','PE_GQ','SR_GT','SR_GQ']), file=fo)
        for i in svid_gt.keys():
          print('\t'.join([i]+svid_gt[i]), file=fo)
        fo.close()
        CODE

    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}

task ExtracAlgorithmEvidenceFilter{
  input{
    String prefix
    File vcf_file
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

    RuntimeAttr default_attr = object {
        cpu_cores: 1, 
        mem_gb: 3.75, 
        disk_gb: 10,
        boot_disk_gb: 10,
        preemptible_tries: 1,
        max_retries: 1
    }
  
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

    output{
        File vcf_info = "~{prefix}.info"
    }

    command <<<
        zcat ~{vcf_file} | grep -v '##' | cut -f3,7  > ~{prefix}.SVID_filter
        svtk vcf2bed -i SVTYPE -i SVLEN -i ALGORITHMS -i EVIDENCE ~{vcf_file} ~{prefix}.bed
        paste  <(cut -f4,7-10 ~{prefix}.bed) \
               <(cut -f2 ~{prefix}.SVID_filter) \
               > ~{prefix}.info
    >>>

    runtime {
        cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
        memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
        disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
        docker: sv_pipeline_docker
        preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    }
}
