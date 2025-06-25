version 1.0

workflow ExtractTriosFromVCF {
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
      docker_image = sv_pipeline_base_docker
  }


  call split_vcf_by_annotation{
    input:
      vcf_file = input_vcf,
      svid_annotation = AnnotateGenomicContext.variant_anno, 
      docker_image = sv_pipeline_base_docker
  }


  scatter (family in families) {
    call WriteTrioSampleFile {
      input:
        family = family,
        docker_image = sv_pipeline_base_docker
    }

    call ExtractTrioVCF {
      input:
        input_vcf = input_vcf,
        sample_file = WriteTrioSampleFile.sample_file,
        family_id = WriteTrioSampleFile.family_id,
        docker_image = sv_pipeline_base_docker
    }

    call SplitVariantsBySize {
      input:
        input_vcf = ExtractTrioVCF.output_vcf,
        docker_image = sv_pipeline_base_docker
    }

    call CalculateInheritanceTable as calcu_inheri_table_snv{
      input:
        input_vcf = SplitVariantsBySize.snv_vcf,
        input_vcf_idx = SplitVariantsBySize.snv_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_snv
    }

    call CalculateInheritanceTable as calcu_inheri_table_indel_sm{
      input:
        input_vcf = SplitVariantsBySize.indel_1_30_vcf,
        input_vcf_idx = SplitVariantsBySize.indel_1_30_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_indel_sm
    }

    call CalculateInheritanceTable as calcu_inheri_table_indel_lg{
      input:
        input_vcf = SplitVariantsBySize.indel_31_50_vcf,
        input_vcf_idx = SplitVariantsBySize.indel_31_50_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_indel_lg
    }
    
    call CalculateInheritanceTable as calcu_inheri_table_sv{
      input:
        input_vcf = SplitVariantsBySize.sv_vcf,
        input_vcf_idx = SplitVariantsBySize.sv_vcf_idx,
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_ovr_calcu_inheri_table_sv
    }

    call IntegrateInheritanceTable{
      input:
        inheri_table_snv = calcu_inheri_table_snv.inheri_stat,
        inheri_table_indel_sm = calcu_inheri_table_indel_sm.inheri_stat,
        inheri_table_indel_lg = calcu_inheri_table_indel_lg.inheri_stat,
        inheri_table_sv = calcu_inheri_table_sv.inheri_stat,
        inheri_table_ref = inheri_table,
        family_id = WriteTrioSampleFile.family_id,
        prefix   = prefix,
        docker_image = sv_pipeline_base_docker
    }    

  }

  output{
    Array[File] inheritance_table_inte = IntegrateInheritanceTable.integrated_inheri_stat
  }
}

task WriteTrioSampleFile {
  input {
    Array[String] family  # [father, mother, child]
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command {
    echo -e "${family[0]}\n${family[1]}\n${family[2]}" > trio_${family[0]}_${family[2]}.samples.txt
  }

  output {
    File sample_file = "trio_${family[0]}_${family[2]}.samples.txt"
    String family_id = "${family[0]}_${family[2]}"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 2,
    disk_gb: 5,
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

task ExtractTrioVCF {
  input {
    File input_vcf
    File sample_file
    String family_id
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command {
    bcftools view -S ~{sample_file} -Oz -o ~{family_id}.trio.vcf.gz ~{input_vcf}
  }

  output {
    File output_vcf = "~{family_id}.trio.vcf.gz"
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

task CalculateInheritanceTable {
  input {
    File input_vcf
    File input_vcf_idx
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 + ceil(size(input_vcf,"GB") * 2)
  Int mem_size =  ceil(size(input_vcf,"GB") * 2)

  String prefix = basename(input_vcf, ".vcf.gz")
  command <<<
    set -euxo pipefail

    bcftools view -H ~{input_vcf} | cut -f10- | sort | uniq -c > ~{prefix}.inheri.stat

  >>>

  output {
    File inheri_stat = "~{prefix}.inheri.stat"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_size,
    disk_gb: disk_size,
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

task IntegrateInheritanceTable {
  input {
    File inheri_table_snv
    File inheri_table_indel_sm
    File inheri_table_indel_lg
    File inheri_table_sv
    File inheri_table_ref
    String family_id
    String prefix

    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  Int disk_size = 10 
  Int mem_size =  15

  command <<<
    set -euxo pipefail

    Rscript -e '

    process_and_collapse <- function(df) {
      # Step 1: Replace "." with "0" in columns 2 to 4
      df[, 2:4] <- lapply(df[, 2:4], function(col) gsub("\\.", "0", col))
      
      # Step 2: Collapse rows by columns 2–4, summing column 1
      # Convert column 1 to numeric (if not already)
      df[, 1] <- as.numeric(df[, 1])
      
      # Use aggregate to group and sum
      collapsed_df <- aggregate(df[, 1], by = df[, 2:4], FUN = sum)
      
      # Rename the columns appropriately
      colnames(collapsed_df) <- c("fa", "mo", "pb", "count_variants")
      
      return(collapsed_df)
    }

    read_and_merge_five_tables <- function(file1, file2, file3, file4, file5) {
      # Read tables 1-4 without headers
      t1 <- read.table(file1, header = FALSE, stringsAsFactors = FALSE)
      t2 <- read.table(file2, header = FALSE, stringsAsFactors = FALSE)
      t3 <- read.table(file3, header = FALSE, stringsAsFactors = FALSE)
      t4 <- read.table(file4, header = FALSE, stringsAsFactors = FALSE)

      # Process tables 1–4
      p1 <- process_and_collapse(t1)
      p2 <- process_and_collapse(t2)
      p3 <- process_and_collapse(t3)
      p4 <- process_and_collapse(t4)

      # Rename count columns for clarity
      colnames(p1)[4] <- "count_snv"
      colnames(p2)[4] <- "count_indel_sm"
      colnames(p3)[4] <- "count_indel_lg"
      colnames(p4)[4] <- "count_sv"

      # Merge processed tables by col2, col3, col4
      merged <- Reduce(function(x, y) merge(x, y, by = c("fa", "mo", "pb"), all = TRUE),
                       list(p1, p2, p3, p4))

      # Read table 5 with header
      t5 <- read.table(file5, header = TRUE, stringsAsFactors = FALSE)

      # Final merge with table 5
      final <- merge(merged, t5, by.x = c("fa", "mo", "pb"), by.y = c("fa", "mo", "pb"), all = TRUE)

      return(final)
    }

    file1 = "~{inheri_table_snv}"
    file2 = "~{inheri_table_indel_sm}"
    file3 = "~{inheri_table_indel_lg}"
    file4 = "~{inheri_table_sv}"
    file5 = "~{inheri_table_ref}"

    integrated_inheri_table = read_and_merge_five_tables(file1, file2, file3, file4, file5)
    write.table(integrated_inheri_table, "~{family_id}.~{prefix}.integrated_inheritance_table.tsv", quote=F, sep="\t", col.names=T, row.names=F)
    
    '

  >>>

  output {
    File integrated_inheri_stat = "~{family_id}.~{prefix}.integrated_inheritance_table.tsv"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: mem_size,
    disk_gb: disk_size,
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

task split_vcf_by_annotation {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    RuntimeAttr? runtime_attr_override
  }


  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    pip install pysam

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
    annotation_map = {}
    with open("~{svid_annotation}", "r") as f:
        for line in f:
            if line.strip():
                svid, annotation = line.strip().split('\t')
                annotation_map[svid] = annotation

    # Prepare output VCF writers per annotation
    vcf_in = pysam.VariantFile("~{vcf_file}", "r")
    vcf_outs = defaultdict(lambda: None)

    for rec in vcf_in.fetch():
        chrom = rec.chrom
        pos = str(rec.pos)
        ref = rec.ref
        end = str(rec.pos + len(ref) -1)
        alt = rec.alts[0] if rec.alts else "N"

        svtype = determine_svtype(ref, alt)
        svlen = determine_svlen(ref, alt)
        svid = f"{chrom}_{pos}_{end}_{svtype}_{svlen}"

        if svid in annotation_map:
            annotation = annotation_map[svid]
            if vcf_outs[annotation] is None:
                out_path = f"~{prefix}.{annotation}.vcf"
                vcf_outs[annotation] = pysam.VariantFile(out_path, "w", header=vcf_in.header)
            vcf_outs[annotation].write(rec)

    # Close all output files
    for writer in vcf_outs.values():
        writer.close()

    # Save output file list
    with open("split_outputs.txt", "w") as out_list:
        for annotation in vcf_outs.keys():
            out_list.write(f"{annotation}.vcf\n")
    CODE
  >>>

  output {
    Array[File] split_vcfs = read_lines("split_outputs.txt")
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
