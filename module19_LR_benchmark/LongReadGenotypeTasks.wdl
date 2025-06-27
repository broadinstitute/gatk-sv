version 1.0


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

task AddGenomicContextToVcfPython {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
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
    svid_anno = {}

    print('reading in SVID and their genomic context ... ')
    with open("~{svid_annotation}", "r") as f:
        for line in f:
            if line.strip():
                svid, annotation = line.strip().split('\t')
                svid_anno[svid] = annotation

    print('complete')

    # Prepare output VCF writers per annotation
    vcf_in = pysam.VariantFile("~{vcf_file}", "r")

    header = vcf_in.header
    header.info.add("GC", number=1, type="String", description="Genomic context of the variant")

    vcf_out = pysam.VariantFile("~{prefix}.GC_anno.vcf.gz", 'w', header = header)

    for rec in vcf_in.fetch():
        svid = rec.id
        rec_new = rec
        rec_new.info["GC"] = svid_anno[svid]
        vcf_out.write(rec_new)

    vcf_out.close()

    CODE
  >>>

  output {
    File annotated_vcf = "~{prefix}.GC_anno.vcf.gz"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 15 + ceil(size(vcf_file)*2),
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

task AddGenomicContextToVcfR {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    set -e

    # use R script to add GC to the vcf
    R --vanilla <<EOF

    svid_gc <- read.table("~{svid_annotation}", header = TRUE)
    vcf_in <- read.table("~{vcf_file}", header = FALSE)
    colnames(vcf_in)[3] = 'SVID'
    vcf_out <- merge(vcf_in, svid_gc, by='SVID')
    vcf_out[,8] = paste(vcf_out[,8], paste('GC',vcf_out$GC, sep='='),sep=';')
    vcf_out_v2 = vcf_out[,c(2,3,1,4:(ncol(vcf_out)-1))]

    write.table(vcf_out, file = "~{prefix}.GC_anno.vcf", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)

    EOF

    #use python script to add GC to vcf header
    python3 <<CODE

    import pysam
    import os

    vcf_in = pysam.VariantFile("~{vcf_file}", "r")
    header = vcf_in.header
    header.info.add("GC", number=1, type="String", description="Genomic context of the variant")
    vcf_out = pysam.VariantFile("~{prefix}.header.vcf.gz", 'w', header = header)
    vcf_in.close()
    vcf_out.close()
    
    CODE

    cat <(zcat ~{prefix}.header.vcf.gz) ~{prefix}.GC_anno.vcf | bgzip > ~{prefix}.GC_anno.vcf.gz
    tabix -p vcf ~{prefix}.GC_anno.vcf.gz

  >>>


  output {
    File annotated_vcf = "~{prefix}.GC_anno.vcf.gz"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 15,
    disk_gb: 15 + ceil(size(vcf_file)*2),
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

task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File]? vcfs_idx
    Boolean merge_sort = false
    String? outfile_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = outfile_prefix + ".vcf.gz"
  String merge_flag = if merge_sort then "--allow-overlaps" else ""

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size(vcfs, "GB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    VCFS="~{write_lines(vcfs)}"
    if ~{!defined(vcfs_idx)}; then
      cat ${VCFS} | xargs -n1 tabix
    fi
    bcftools concat -a ~{merge_flag} --output-type z --file-list ${VCFS} --output "~{outfile_name}"
    tabix -p vcf -f "~{outfile_name}"
  >>>

  output {
    File concat_vcf = outfile_name
    File concat_vcf_idx = outfile_name + ".tbi"
  }
}

task MergeVcfs {
  input {
    Array[File] input_vcfs
    String output_name
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e

    mkdir -p indexed_vcfs

    # Ensure all input VCFs are indexed
    for vcf in ~{sep=' ' input_vcfs}; do
      cp "$vcf" indexed_vcfs/
      vcf_basename=$(basename "$vcf")
      if [ ! -f "${vcf}.tbi" ]; then
        echo "Index not found for $vcf, creating it..."
        tabix -p vcf "$vcf"
        cp "${vcf}.tbi" "indexed_vcfs/${vcf_basename}.tbi"
      else
        cp "${vcf}.tbi" "indexed_vcfs/${vcf_basename}.tbi"
      fi
    done

    # Merge the input VCFs
    bcftools concat -a \
      ~{sep=' ' input_vcfs} \
      -Oz -o merged.tmp.vcf.gz

    # Index the temp merged VCF
    tabix -p vcf merged.tmp.vcf.gz

    # Remove duplicates
    bcftools norm -d all \
      -Oz -o ~{output_name} merged.tmp.vcf.gz

    # Index the final deduplicated VCF
    tabix -p vcf ~{output_name}
  >>>

  output {
    File merged_vcf = "~{output_name}"
    File merged_vcf_index = "~{output_name}.tbi"
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
    #rewrite the vcfs in the a new file with the SVID updated to be consistent with the bed file
    vcf_out= pysam.VariantFile("~{prefix}.SVID_updated.vcf.gz", 'w', header = vcf_in.header)
    
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
            rec.id = ID
            out.write('\t'.join([chrom, pos, end, ID, svtype, svlen]) + "\n")
            vcf_out.write(rec)
    CODE

    tabix -p vcf "~{prefix}.SVID_updated.vcf.gz"

  >>>

  output {
    File variant_sites = "~{prefix}.variant_sites.bed"
    File updated_vcf = "~{prefix}.SVID_updated.vcf.gz"
    File updated_vcf_idx = "~{prefix}.SVID_updated.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 20,
    disk_gb: 20 + ceil(size(input_vcf,"GiB")*2),
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

task ExtractChromosomeVariants {
  input {
    File input_vcf       # .vcf.gz
    File input_vcf_index # .vcf.gz.tbi
    String chromosome    # e.g. "chr1" or "1"
    String output_name   # e.g. "chr1.vcf.gz"
    String docker_image = "biocontainers/bcftools:v1.17-1-deb-py3"
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e

    # Extract chromosome-specific variants
    bcftools view -r ~{chromosome} ~{input_vcf} -Oz -o ~{output_name}

    # Index the output VCF
    tabix -p vcf ~{output_name}
  >>>

  output {
    File chr_vcf = output_name
    File chr_vcf_index = output_name + ".tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb:  10,
    disk_gb: 10 + ceil(size(input_vcf,"GiB")*1.5),
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

task SplitVcfIntoShards {
  input {
    File input_vcf       # .vcf.gz
    File input_vcf_index # .vcf.gz.tbi
    Int variants_per_shard    # number of varians in each splitted variants
    String output_prefix   # e.g. "chr1.vcf.gz"
    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -e

    mkdir chunks

    # Extract header
    bcftools view -h ~{input_vcf} > chunks/header.vcf

    # Output body lines (no header), then split
    bcftools view -H ~{input_vcf} | split -l ~{variants_per_shard} - chunks/body_

    # Reconstruct chunked VCFs with header
    for body in chunks/body_*; do
      chunk_name=chunks/~{output_prefix}_$(basename "$body")
      cat chunks/header.vcf "$body" | bgzip -c > "${chunk_name}.vcf.gz"
      tabix -p vcf "${chunk_name}.vcf.gz"
    done
 
   >>>

  output {
    Array[File] split_vcfs = glob("chunks/*.vcf.gz")
    Array[File] split_vcf_indexes = glob("chunks/*.vcf.gz.tbi")
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb:  20,
    disk_gb: 20 + ceil(size(input_vcf,"GiB")*5),
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

task SplitVcfByAnnotationPython {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File vcf_idx
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

    print('reading in SVIDs ... ')
    with open("~{svid_annotation}", "r") as f:
        for line in f:
            if line.strip():
                svid, annotation = line.strip().split('\t')
                if annotation in annotation_list:
                  svid_list.append(svid)
    print('complete !')

    # Prepare output VCF writers per annotation
    vcf_in = pysam.VariantFile("~{vcf_file}", "r")
    vcf_out = pysam.VariantFile("~{prefix}.~{midfix}.vcf.gz", 'w', header = vcf_in.header)

    for rec in vcf_in.fetch():
        if rec.id in svid_list:
            print(svid)
            vcf_out.write(rec)

    vcf_out.close()

    CODE

    tabix -p vcf  "~{prefix}.~{midfix}.vcf.gz"
  >>>

  output {
    File split_vcf = "~{prefix}.~{midfix}.vcf.gz"
    File split_vcf_idx = "~{prefix}.~{midfix}.vcf.gz.tbi"
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

task SplitVcfByAnnotationR {
  input {
    File vcf_file         # Input VCF file (bgzipped or plain)
    File vcf_idx
    File svid_annotation  # 2-column TSV: SVID <tab> annotation
    String docker_image
    String anno_list
    String midfix
    RuntimeAttr? runtime_attr_override
  }

  String prefix = basename(vcf_file,'.vcf.gz')

  command <<<
    set -e

    R --vanilla <<EOF

    annotation_list <- strsplit("~{anno_list}", ',')[[1]]
    svid_gc <- read.table("~{svid_annotation}", header = TRUE)
    svid_list <- svid_gc[svid_gc[,2] %in% annotation_list, 1]
    vcf_in <- read.table("~{vcf_file}", header = FALSE)
    vcf_out <- vcf_in[vcf_in[,3] %in% svid_list, ]
    write.table(vcf_out, file = "~{prefix}.~{midfix}.vcf", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
    
    EOF

    # Re-attach VCF header and compress + index
    
    bcftools view -h ~{vcf_file} | cat - ~{prefix}.~{midfix}.vcf | bgzip > ~{prefix}.~{midfix}.vcf.gz
    tabix -p vcf ~{prefix}.~{midfix}.vcf.gz

  >>>

  output {
    File split_vcf = "~{prefix}.~{midfix}.vcf.gz"
    File split_vcf_idx = "~{prefix}.~{midfix}.vcf.gz.tbi"
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10 + ceil(size(vcf_file, "GiB")*5),
    disk_gb: 15 + ceil(size(vcf_file, "GiB")*5),
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
    disk_gb: 20 + ceil(size(input_vcf, "GiB")*5),
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
      for(i in c(2:4)){
        df[,i] = sapply(df[,i], function(x){strsplit(as.character(x),":")[[1]][1]})
      }
      df = df[df[,2]!="./." | df[,3]!="./." | df[,4]!="./.",]

      df[, 2:4] <- lapply(df[, 2:4], function(col) gsub("/", "|", col))
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

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
}


