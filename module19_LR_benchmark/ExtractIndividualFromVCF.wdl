#wdl script to extract individual vcfs from joint vcfs that include variants of all sizes, ranging from 1bp SNVs to large SVs;
#this script splits variants in each individual to SNVs, Indels 1-30bp, Indels 30-50bp, Indels 1-50bp, SVs >30bp, and SVs>50bp

version 1.0

workflow ExtractIndividualFromVCF {
  input {
    File vcf_file
    Array[String] sample_ids
    String sv_pipeline_base_docker
  }

  scatter (sample_id in sample_ids) {
    call ExtractVariantIndividualGenome {
      input:
        vcf_file = vcf_file,
        sample_id = sample_id,
        docker_image = sv_pipeline_base_docker
    }

    call SplitVariantsBySize{
      input:
        input_vcf = ExtractVariantIndividualGenome.non_ref_vcf,
        docker_image = sv_pipeline_base_docker
    }

    call ConcatVcfs as concat_vcf_1to50bp{
      input:
        vcf1 = SplitVariantsBySize.indel_1_30_vcf,
        vcf2 = SplitVariantsBySize.indel_31_50_vcf,
        idx1 = SplitVariantsBySize.indel_1_30_vcf_idx,
        idx2 = SplitVariantsBySize.indel_31_50_vcf_idx,
        input_vcf = ExtractVariantIndividualGenome.non_ref_vcf,
        appendix = 'indel_1_50',
        docker_image = sv_pipeline_base_docker
    }

    call ConcatVcfs as concat_vcf_over30bp{
      input:
        vcf1 = SplitVariantsBySize.indel_31_50_vcf,
        vcf2 = SplitVariantsBySize.sv_vcf,
        idx1 = SplitVariantsBySize.indel_31_50_vcf_idx,
        idx2 = SplitVariantsBySize.sv_vcf_idx,
        input_vcf = ExtractVariantIndividualGenome.non_ref_vcf,
        appendix = 'sv_over30',
        docker_image = sv_pipeline_base_docker
    }
  }


  output {
    Array[File] all_snv_vcfs = SplitVariantsBySize.snv_vcf  
    Array[File] all_indel_1_30 = SplitVariantsBySize.indel_1_30_vcf
    Array[File] all_indel_31_50 = SplitVariantsBySize.indel_31_50_vcf
    Array[File] all_sv_over50 = SplitVariantsBySize.sv_vcf
    Array[File] all_sv_over30 = concat_vcf_over30bp.combined_vcf
    Array[File] all_indel_1_50 = concat_vcf_1to50bp.combined_vcf
  }
}

task ConcatVcfs {
  input {
    File vcf1
    File vcf2
    File idx1
    File idx2
    File input_vcf
    String appendix
    String docker_image
  }

  String prefix = basename(input_vcf, '.vcf.gz')

  command <<<
    set -euo pipefail

    bcftools concat -a -Oz -o "~{prefix}.~{appendix}.vcf.gz" ~{vcf1} ~{vcf2}
    tabix -p vcf "~{prefix}.~{appendix}.vcf.gz"

  >>>

  output {
    File combined_vcf = "~{prefix}.~{appendix}.vcf.gz"
    File combined_vcf_index = "~{prefix}.~{appendix}.vcf.gz.tbi"
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "2G"
  }
}

task ExtractVariantIndividualGenome {
  input {
    File vcf_file
    String sample_id
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euxo pipefail

    # Extract variants for the sample
    bcftools view -s ~{sample_id} ~{vcf_file} -Oz -o ~{sample_id}.vcf.gz
    tabix -p vcf ~{sample_id}.vcf.gz

    # Filter informative, alternative genotypes
    bcftools view -e 'GT="0|0" || GT=".|." || GT=".|0" || GT="0|."' ~{sample_id}.vcf.gz -Oz -o ~{sample_id}.non_ref.vcf.gz
    tabix -p vcf ~{sample_id}.non_ref.vcf.gz

  >>>

  output {
    File non_ref_vcf = "~{sample_id}.non_ref.vcf.gz"
    File non_ref_vcf_idx = "~{sample_id}.non_ref.vcf.gz.tbi"
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

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task SplitVariantsBySize {
  input {
    File input_vcf
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(input_vcf, ".vcf.gz")
  command <<<
    set -euxo pipefail

    # Index VCF if needed
    if [[ ! -f "~{input_vcf}.tbi" && ! -f "~{input_vcf}.csi" ]]; then
      bcftools index ~{input_vcf}
    fi

    python3 <<CODE
  
    import pysam
    import os

    def get_variant_type_size(ref, alt):
        len_ref = len(ref)
        len_alt = len(alt)
        svlen = abs(len_ref - len_alt)

        if len_ref == 1 and len_alt == 1:
            return "SNV", 0
        elif svlen <= 30:
            return "INDEL_1_30", svlen
        elif 30 < svlen <= 50:
            return "INDEL_30_50", svlen
        else:
            return "SV_GT_50", svlen

    def main(input_vcf):
        # Open input VCF
        vcf_in = pysam.VariantFile(input_vcf)

        # Prepare output files
        outputs = {
            "SNV": pysam.VariantFile(f"snvs.vcf", 'w', header=vcf_in.header),
            "INDEL_1_30": pysam.VariantFile(f"indels_1_30.vcf", 'w', header=vcf_in.header),
            "INDEL_30_50": pysam.VariantFile(f"indels_30_50.vcf", 'w', header=vcf_in.header),
            "SV_GT_50": pysam.VariantFile(f"svs_gt_50.vcf", 'w', header=vcf_in.header)
        }

        for rec in vcf_in.fetch():
            if len(rec.alts) != 1:
                continue  # skip multiallelics for now
            ref = rec.ref
            alt = rec.alts[0]
            vtype, svlen = get_variant_type_size(ref, alt)

            outputs[vtype].write(rec)

        # Close all files
        for f in outputs.values():
            f.close()
        vcf_in.close()

    main("~{input_vcf}")
    
    CODE



    mv snvs.vcf "~{prefix}.snv.vcf"
    mv indels_1_30.vcf "~{prefix}.indel_1_30.vcf"
    mv indels_30_50.vcf "~{prefix}.indel_31_50.vcf"
    mv svs_gt_50.vcf "~{prefix}.sv_over50.vcf"

    bgzip "~{prefix}.snv.vcf"
    bgzip "~{prefix}.indel_1_30.vcf"
    bgzip "~{prefix}.indel_31_50.vcf"
    bgzip "~{prefix}.sv_over50.vcf"

    tabix -p vcf "~{prefix}.snv.vcf.gz"
    tabix -p vcf "~{prefix}.indel_1_30.vcf.gz"
    tabix -p vcf "~{prefix}.indel_31_50.vcf.gz"
    tabix -p vcf "~{prefix}.sv_over50.vcf.gz"

  >>>

  output {
    File snv_vcf = "~{prefix}.snv.vcf.gz"
    File indel_1_30_vcf = "~{prefix}.indel_1_30.vcf.gz"
    File indel_31_50_vcf = "~{prefix}.indel_31_50.vcf.gz"
    File sv_vcf = "~{prefix}.sv_over50.vcf.gz"

    File snv_vcf_idx = "~{prefix}.snv.vcf.gz.tbi"
    File indel_1_30_vcf_idx = "~{prefix}.indel_1_30.vcf.gz.tbi"
    File indel_31_50_vcf_idx = "~{prefix}.indel_31_50.vcf.gz.tbi"
    File sv_vcf_idx = "~{prefix}.sv_over50.vcf.gz.tbi"
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

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task ExtractVariantSites {
  input {
    File input_vcf
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(input_vcf, '.vcf.gz')
  command <<<
    set -euxo pipefail

    # Output header
    echo -e "CHROM\tPOS\tEND\tID\tSVTYPE\tSVLEN" > ~{prefix}.variant_sites.tsv

    # Use awk to calculate END = POS + length(REF) - 1
    bcftools query -f '%CHROM\t%POS\t%REF\t%ID\t%INFO/SVTYPE\t%INFO/SVLEN\n' ~{input_vcf} | \
      awk '{
        ref_len = length($3);
        end_pos = $2 + ref_len - 1;
        print $1"\t"$2"\t"end_pos"\t"$4"\t"$5"\t"$6
      }' >> ~{prefix}.variant_sites.bed
  >>>

  output {
    File variant_sites = "~{prefix}.variant_sites.bed"
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

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task AnnotateGenomicContext {
  input {
    File variant_sites
    File anno_script_Rscript
    File anno_script_bash
    File repeat_mask
    File segmental_duplicates
    File simple_repeats
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(variant_sites, '.bed')
  command <<<
    set -euxo pipefail

    bash ~{anno_script_bash} ~{variant_sites} ~{prefix}.SVID_GC \
    --rm ~{repeat_mask} \
    --sd ~{segmental_duplicates} \
    --sr ~{simple_repeats}

  >>>

  output {
    File variant_anno = "~{prefix}.SVID_GC"
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

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: docker_image
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}