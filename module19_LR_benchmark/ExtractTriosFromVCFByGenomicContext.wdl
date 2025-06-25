version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/LR_genotype:ExtractTriosFromVCF/versions/14/plain-WDL/descriptor" as ExtractTriosFromVCF

workflow ExtractTriosFromVCFByGenomicContext {
  input {
    File input_vcf
    Array[Array[String]] families  # List of trios: [[fa1, mo1, ch1], [fa2, mo2, ch2], ...]
    String sv_pipeline_base_docker 
    File inheri_table
    String prefix

    File anno_script_bash
    File anno_script_Rscript
    File repeat_mask
    File simple_repeats
    File segmental_duplicates

    RuntimeAttr? runtime_attr_override
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_snv
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_sv
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_indel_lg
    RuntimeAttr? runtime_attr_ovr_calcu_inheri_table_indel_sm
    RuntimeAttr? runtime_attr_annotate_genomic_context
  }

  call ExtractVariantSites{
    input:
      input_vcf = input_vcf,
      docker_image = sv_pipeline_base_docker
  }

  call AnnotateGenomicContext{
    input:
      variant_sites = ExtractVariantSites.variant_sites,
      anno_script_bash = anno_script_bash,
      anno_script_Rscript = anno_script_Rscript,
      repeat_mask = repeat_mask,
      segmental_duplicates = segmental_duplicates,
      simple_repeats = simple_repeats,
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_annotate_genomic_context
  }

  call SplitVcfByAnnotation as split_vcf_by_annotation_sr {
    input:
      vcf_file = input_vcf,
      svid_annotation = AnnotateGenomicContext.variant_anno, 
      anno_list = "SR",
      midfix = "SR",
      docker_image = sv_pipeline_base_docker
  }

  call SplitVcfByAnnotation as split_vcf_by_annotation_sd {
    input:
      vcf_file = input_vcf,
      svid_annotation = AnnotateGenomicContext.variant_anno, 
      anno_list = "SD",
      midfix = "SD",
      docker_image = sv_pipeline_base_docker
  }

  call SplitVcfByAnnotation as split_vcf_by_annotation_clean {
    input:
      vcf_file = input_vcf,
      svid_annotation = AnnotateGenomicContext.variant_anno, 
      anno_list = "US,RM",
      midfix = "US_RM",
      docker_image = sv_pipeline_base_docker
  }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_SRs{
    input:
      input_vcf = split_vcf_by_annotation_sr.split_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.SR",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_ovr_calcu_inheri_table_sv = runtime_attr_ovr_calcu_inheri_table_sv,
      runtime_attr_ovr_calcu_inheri_table_snv = runtime_attr_ovr_calcu_inheri_table_snv,
      runtime_attr_ovr_calcu_inheri_table_indel_sm = runtime_attr_ovr_calcu_inheri_table_indel_sm,
      runtime_attr_ovr_calcu_inheri_table_indel_lg = runtime_attr_ovr_calcu_inheri_table_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_SDs{
    input:
      input_vcf = split_vcf_by_annotation_sd.split_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.SD",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_ovr_calcu_inheri_table_sv = runtime_attr_ovr_calcu_inheri_table_sv,
      runtime_attr_ovr_calcu_inheri_table_snv = runtime_attr_ovr_calcu_inheri_table_snv,
      runtime_attr_ovr_calcu_inheri_table_indel_sm = runtime_attr_ovr_calcu_inheri_table_indel_sm,
      runtime_attr_ovr_calcu_inheri_table_indel_lg = runtime_attr_ovr_calcu_inheri_table_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_USRMs{
    input:
      input_vcf = split_vcf_by_annotation_clean.split_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.US_RM",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_ovr_calcu_inheri_table_sv = runtime_attr_ovr_calcu_inheri_table_sv,
      runtime_attr_ovr_calcu_inheri_table_snv = runtime_attr_ovr_calcu_inheri_table_snv,
      runtime_attr_ovr_calcu_inheri_table_indel_sm = runtime_attr_ovr_calcu_inheri_table_indel_sm,
      runtime_attr_ovr_calcu_inheri_table_indel_lg = runtime_attr_ovr_calcu_inheri_table_indel_lg
    }

  call ExtractTriosFromVCF.ExtractTriosFromVCF as extract_trios_from_all{
    input:
      input_vcf = input_vcf,
      families = families,
      inheri_table = inheri_table,
      prefix = "~{prefix}.all",
      sv_pipeline_base_docker = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_override,
      runtime_attr_ovr_calcu_inheri_table_sv = runtime_attr_ovr_calcu_inheri_table_sv,
      runtime_attr_ovr_calcu_inheri_table_snv = runtime_attr_ovr_calcu_inheri_table_snv,
      runtime_attr_ovr_calcu_inheri_table_indel_sm = runtime_attr_ovr_calcu_inheri_table_indel_sm,
      runtime_attr_ovr_calcu_inheri_table_indel_lg = runtime_attr_ovr_calcu_inheri_table_indel_lg
    }
  
  output{
    Array[File] inheritance_table_inte_SR = extract_trios_from_SRs.inheritance_table_inte
    Array[File] inheritance_table_inte_SD = extract_trios_from_SDs.inheritance_table_inte
    Array[File] inheritance_table_inte_US_RM = extract_trios_from_USRMs.inheritance_table_inte
    Array[File] inheritance_table_inte_all = extract_trios_from_all.inheritance_table_inte
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

    python3 <<CODE
    import pysam

    def determine_svtype(ref, alt):
        if len(ref) == 1 and len(alt) == 1:
            return "SNV"
        elif len(ref) < len(alt):
            return "INS"
        elif len(ref) > len(alt):
            return "DEL"
        else:
            return "OTH"

    def determine_svlen(ref, alt):
        return str(abs(len(alt) - len(ref)))

    vcf_in = pysam.VariantFile("~{input_vcf}")
    with open("~{prefix}.variant_sites.bed", "w") as out:
        for rec in vcf_in.fetch():
            chrom = rec.contig
            pos = str(rec.pos)
            ref = rec.ref
            end = str(rec.pos + len(ref) -1)  # .stop is preferred over parsing INFO['END']
            alt = rec.alts[0] if rec.alts else "."
            svtype = determine_svtype(ref, alt)
            svlen = determine_svlen(ref, alt)
            ID = f"{chrom}_{pos}_{end}_{svtype}_{svlen}"
            out.write('\t'.join([chrom, pos, end, ID, svtype, svlen]) + "\n")
    CODE

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
    mem_gb: ceil(size(variant_sites, "GiB"))*10,
    disk_gb: 10 + ceil(size(variant_sites, "GiB"))*10,
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

task SplitVcfByAnnotation {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    String anno_list
    String midfix
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<

    python3 <<CODE
    import pysam
    import os
    from collections import defaultdict

    def determine_svtype(ref, alt):
        if len(ref) == 1 and len(alt) == 1:
            return "SNV"
        elif len(ref) < len(alt):
            return "INS"
        elif len(ref) > len(alt):
            return "DEL"
        else:
            return "OTH"

    def determine_svlen(ref, alt):
        return str(abs(len(alt) - len(ref)))

    # Load SVID -> annotation mapping
    annotation_list = "~{anno_list}".split(',')
    svid_list = []
    with open("~{svid_annotation}", "r") as f:
        for line in f:
            if line.strip():
                svid, annotation = line.strip().split('\t')
                if annotation in annotation_list:
                  svid_list.append(svid)

    # Prepare output VCF writers per annotation
    vcf_in = pysam.VariantFile("~{vcf_file}", "r")
    vcf_out = pysam.VariantFile("~{prefix}.~{midfix}.vcf.gz", 'w', header = vcf_in.header)

    for rec in vcf_in.fetch():
        chrom = rec.chrom
        pos = str(rec.pos)
        ref = rec.ref
        end = str(rec.pos + len(ref) -1)
        alt = rec.alts[0] if rec.alts else "N"

        svtype = determine_svtype(ref, alt)
        svlen = determine_svlen(ref, alt)
        svid = f"{chrom}_{pos}_{end}_{svtype}_{svlen}"

        if svid in svid_list:
            vcf_out.write(rec)

    vcf_out.close()

    CODE
  >>>

  output {
    File split_vcf = "~{prefix}.~{midfix}.vcf.gz"
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