version 1.0

workflow PhysicalPhasingPerContig {

    input {
        File all_chr_bam
        File all_chr_bai
        File reference_fasta
        File reference_fasta_fai
        File small_vcf
        File small_vcf_idx
        File sv_vcf
        File sv_vcf_idx
        File? trgt_vcf
        File? trgt_vcf_idx
        Int hiphase_memory
        String hiphase_extra_args
        String sample_id
        String region

        String sv_pipeline_base_docker
    }

    #call PreProcessVCF as PreProcessSVVCF { 
    #    input:
    #        vcf = sv_vcf,
    #        vcf_idx = sv_vcf_idx,
    #        prefix = sample_id + ".uppercased.unphased.subset.sv",
    #        locus = region
    #}

    call SubsetVCF as SubsetVcfSv{
        input:
            vcf_gz = sv_vcf,
            vcf_idx = sv_vcf_idx,
            ocus = region,
            docker_image = sv_pipeline_base_docker
    }
    call SubsetVCF as SubsetVcfShort { 
        input:
            vcf_gz = small_vcf,
            vcf_idx = small_vcf_idx,
            locus = region,
            docker_image = sv_pipeline_base_docker
    }

    if (defined(trgt_vcf) && defined(trgt_vcf_idx)) {
        call SubsetVCF as SubsetVcfTRGT { 
            input:
                vcf_gz = select_first([trgt_vcf]),
                vcf_idx = select_first([trgt_vcf_idx]),
                locus = region,
                docker_image = sv_pipeline_base_docker
        }
        call HiPhaseTRGT as HiphaseTrgt { 
            input:
                bam = all_chr_bam,
                bai = all_chr_bai,
                unphased_sv_vcf = SubsetVcfSv.subset_vcf,
                unphased_sv_idx = SubsetVcfSv.subset_idx,
                unphased_snp_vcf = SubsetVcfShort.subset_vcf,
                unphased_snp_idx = SubsetVcfShort.subset_idx,
                unphased_trgt_vcf = SubsetVcfTRGT.subset_vcf,
                unphased_trgt_idx = SubsetVcfTRGT.subset_idx,
                ref_fasta = reference_fasta,
                ref_fasta_fai = reference_fasta_fai,
                samplename = sample_id,
                memory = hiphase_memory,
                extra_args = hiphase_extra_args
        }
    }
    if (! defined(trgt_vcf)) {
        call HiPhase { 
            input:
                bam = all_chr_bam,
                bai = all_chr_bai,
                unphased_sv_vcf = SubsetVcfSv.subset_vcf,
                unphased_sv_idx = SubsetVcfSv.subset_idx,
                unphased_snp_vcf = SubsetVcfShort.subset_vcf,
                unphased_snp_idx = SubsetVcfShort.subset_idx,
                ref_fasta = reference_fasta,
                ref_fasta_fai = reference_fasta_fai,
                samplename = sample_id,
                memory = hiphase_memory,
                extra_args = hiphase_extra_args
        }
    }


    output {
        File hiphase_short_vcf = select_first([HiPhase.phased_snp_vcf, HiphaseTrgt.phased_snp_vcf])
        File hiphase_short_idx = select_first([HiPhase.phased_snp_vcf_idx, HiphaseTrgt.phased_snp_vcf_idx]) 
        File hiphase_sv_vcf = select_first([HiPhase.phased_sv_vcf, HiphaseTrgt.phased_sv_vcf])
        File hiphase_sv_idx = select_first([HiPhase.phased_sv_vcf_idx, HiphaseTrgt.phased_sv_vcf_idx])
        File? hiphase_trgt_vcf = HiphaseTrgt.phased_trgt_vcf
        File? hiphase_trgt_idx = HiphaseTrgt.phased_trgt_vcf_idx
    }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}

struct DataTypeParameters {
    Int num_shards
    String map_preset
}

task PreProcessVCF {
    input {
        File vcf
        File vcf_idx
        String prefix
        String locus
        RuntimeAttr? runtime_attr_override
    }
    String docker_dir = "/truvari_intrasample"

    command <<<
        set -euxo pipefail

        # convert lower case alleles to upper case
        cp ~{docker_dir}/convert_lower_case.py .

        python convert_lower_case.py -i ~{vcf} -o ~{prefix}.uppercase.vcf
        bgzip ~{prefix}.uppercase.vcf ~{prefix}.uppercase.vcf.gz
        tabix -p vcf ~{prefix}.uppercase.vcf.gz

        # unphase genotypes
        bcftools +setGT ~{prefix}.uppercase.vcf.gz --no-version -Oz -o ~{prefix}.unphased.vcf.gz -- --target-gt a --new-gt u
        bcftools index -t ~{prefix}.unphased.vcf.gz

        # subset VCF
        bcftools view --no-version ~{prefix}.unphased.vcf.gz --regions ~{locus} -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz

        ls -l
        
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_idx = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            50,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "hangsuunc/cleanvcf:v1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task SubsetVCF {

    meta {
        description: "Subset a VCF file to a given locus"
    }

    parameter_meta {
        vcf_gz: "VCF file to be subsetted"
        vcf_idx: "Tabix index for the VCF file"
        locus: "Locus to be subsetted"
        prefix: "Prefix for the output file"
        runtime_attr_override: "Override default runtime attributes"
    }

    input {
        File vcf_gz
        File? vcf_idx
        String locus
        String prefix = "subset"
        String docker_image

        RuntimeAttr? runtime_attr_override
    }

    command <<<
        set -euxo pipefail
        if ! ~{defined(vcf_idx)}; then
            bcftools index ~{vcf_gz}
        fi
        bcftools view --no-version ~{vcf_gz} --regions ~{locus} -Oz -o ~{prefix}.vcf.gz
        bcftools index -t ~{prefix}.vcf.gz
    >>>

    output {
        File subset_vcf = "~{prefix}.vcf.gz"
        File subset_idx = "~{prefix}.vcf.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            100,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             docker_image
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task HiPhase {

    meta {
        description: "Generates phased VCF. Note this runs fast so no need to parallize."
    }


    input {
        File bam
        File bai

        File unphased_snp_vcf
        File unphased_snp_idx
        File unphased_sv_vcf
        File unphased_sv_idx

        File ref_fasta
        File ref_fasta_fai
        String samplename

        Int memory
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        String extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = bam_sz + 200
    Int thread_num = memory/2

    command <<<
        set -euxo pipefail

        touch ~{bai}

        hiphase \
        --threads ~{thread_num} \
        --bam ~{bam} \
        --reference ~{ref_fasta} \
        --global-realignment-cputime 300 \
        --vcf ~{unphased_snp_vcf} \
        --output-vcf ~{samplename}_phased_snp.vcf.gz \
        --vcf ~{unphased_sv_vcf} \
        --output-vcf ~{samplename}_phased_sv.vcf.gz \
        --haplotag-file ~{samplename}_phased_sv_haplotag.tsv \
        --stats-file ~{samplename}.stats.csv \
        --blocks-file ~{samplename}.blocks.tsv \
        --summary-file ~{samplename}.summary.tsv \
        --verbose \
        ~{extra_args}
        
    >>>

    output {
        File phased_snp_vcf = "~{samplename}_phased_snp.vcf.gz"
        File phased_snp_vcf_idx = "~{samplename}_phased_snp.vcf.gz.tbi"
        File phased_sv_vcf   = "~{samplename}_phased_sv.vcf.gz"
        File phased_sv_vcf_idx = "~{samplename}_phased_sv.vcf.gz.tbi"
        File haplotag_file = "~{samplename}_phased_sv_haplotag.tsv"
        File hiphase_stats = "~{samplename}.stats.csv"
        File hiphase_blocks = "~{samplename}.blocks.tsv"
        File hiphase_summary = "~{samplename}.summary.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          thread_num,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/hiphase:v1.5.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task HiPhaseTRGT {

    meta {
        description: "Generates phased VCF. Note this runs fast so no need to parallize."
    }


    input {
        File bam
        File bai

        File unphased_snp_vcf
        File unphased_snp_idx
        File unphased_sv_vcf
        File unphased_sv_idx
        File unphased_trgt_vcf
        File unphased_trgt_idx

        File ref_fasta
        File ref_fasta_fai
        String samplename

        Int memory
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"

        String extra_args

        RuntimeAttr? runtime_attr_override
    }

    Int bam_sz = ceil(size(bam, "GB"))
	Int disk_size = bam_sz + 200 # if bam_sz > 200 then 2*bam_sz else bam_sz + 200
    Int thread_num = memory/2

    command <<<
        set -euxo pipefail

        touch ~{bai}

        hiphase \
        --threads ~{thread_num} \
        --bam ~{bam} \
        --reference ~{ref_fasta} \
        --global-realignment-cputime 300 \
        --vcf ~{unphased_snp_vcf} \
        --output-vcf ~{samplename}_phased_snp.vcf.gz \
        --vcf ~{unphased_sv_vcf} \
        --output-vcf ~{samplename}_phased_sv.vcf.gz \
        --vcf ~{unphased_trgt_vcf} \
        --output-vcf ~{samplename}_phased_trgt.vcf.gz \
        --haplotag-file ~{samplename}_phased_sv_haplotag.tsv \
        --stats-file ~{samplename}.stats.csv \
        --blocks-file ~{samplename}.blocks.tsv \
        --summary-file ~{samplename}.summary.tsv \
        --verbose \
        ~{extra_args}
        
    >>>

    output {
        File phased_snp_vcf = "~{samplename}_phased_snp.vcf.gz"
        File phased_snp_vcf_idx = "~{samplename}_phased_snp.vcf.gz.tbi"
        File phased_sv_vcf   = "~{samplename}_phased_sv.vcf.gz"
        File phased_sv_vcf_idx = "~{samplename}_phased_sv.vcf.gz.tbi"
        File phased_trgt_vcf   = "~{samplename}_phased_trgt.vcf.gz"
        File phased_trgt_vcf_idx = "~{samplename}_phased_trgt.vcf.gz.tbi"
        File haplotag_file = "~{samplename}_phased_sv_haplotag.tsv"
        File hiphase_stats = "~{samplename}.stats.csv"
        File hiphase_blocks = "~{samplename}.blocks.tsv"
        File hiphase_summary = "~{samplename}.summary.tsv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          thread_num,
        mem_gb:             memory,
        disk_gb:            disk_size,
        boot_disk_gb:       100,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/hangsuunc/hiphase:v1.5.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        zones: zones
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}