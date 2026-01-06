version 1.0

import "Structs.wdl"

task SplitVariants {
  input {
    File vcf
    Int n_per_split
    Boolean generate_bca
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75, 
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    Array[File] lt5kb_beds = glob("lt5kb.*")
    Array[File] gt5kb_beds = glob("gt5kb.*")
    Array[File] bca_beds = glob("bca.*")
    Array[File] ins_beds = glob("ins.*")
  }
  command <<<
    set -euo pipefail
    svtk vcf2bed ~{vcf} bed_file.bed
    python /opt/sv-pipeline/04_variant_resolution/scripts/split_variants.py \
      --bed bed_file.bed \
      ~{"--n " + n_per_split} \
      ~{if generate_bca then "--bca" else ""}

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

task SplitVcf {
  input {
    File vcf
    Int n_per_split
    String evidence_type    # pe or sr
    Boolean bgzip           # bgzip output (for SR)
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String output_ext = if bgzip then "vcf.gz" else "vcf"

  output {
    Array[File] vcfs = glob("~{evidence_type}.*.~{output_ext}")
  }
  command <<<

    set -euxo pipefail
    if [[ ~{vcf} == *.gz ]] ; then
      zcat ~{vcf} | sed -n -e '/^#/p' > header.vcf;
      zcat ~{vcf} | sed -e '/^#/d' | split -l ~{n_per_split} - ~{evidence_type}.;
    else
      sed -n -e '/^#/p' ~{vcf} > header.vcf;
      sed -e '/^#/d' ~{vcf} | split -l ~{n_per_split} - ~{evidence_type}.;
    fi
    for f in ~{evidence_type}.*; do cat header.vcf $f ~{if bgzip then "| bgzip -c > $f.vcf.gz" else "> $f.vcf"}; done

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

task AddGenotypes {
  input {
    File vcf
    File genotypes
    File varGQ
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File genotyped_vcf = "~{prefix}.genotyped.vcf.gz"
    File genotyped_vcf_index = "~{prefix}.genotyped.vcf.gz.tbi"
  }
  command <<<

    set -euxo pipefail

    # in some cases a vargq cannot be computed and is returned as '.'. Remove these from the final vcf.
    gzip -cd ~{varGQ} | awk '$5 == "." {print $1}' > bad.vargq.list
    gzip -cd ~{vcf} | { grep -wvf bad.vargq.list || [[ $? == 1 ]]; } | bgzip > clean.vcf.gz
    gzip -cd ~{genotypes} | { grep -wvf bad.vargq.list || [[ $? == 1 ]]; } | bgzip > clean.genotypes.txt.gz
    gzip -cd ~{varGQ} | { grep -wvf bad.vargq.list || [[ $? == 1 ]]; } | bgzip > clean.vargq.txt.gz

    /opt/sv-pipeline/04_variant_resolution/scripts/add_genotypes.py \
      clean.vcf.gz \
      clean.genotypes.txt.gz \
      clean.vargq.txt.gz \
      ~{prefix}.genotyped.vcf

    mkdir tmp
    bcftools sort -T tmp/ ~{prefix}.genotyped.vcf -Oz -o ~{prefix}.genotyped.vcf.gz
    tabix ~{prefix}.genotyped.vcf.gz

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

task MakeSubsetVcf {
  input {
    File vcf
    File bed
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String prefix = basename(bed, ".bed")

  output {
    File subset_vcf = "${prefix}.vcf.gz"
  }
  command <<<

    set -euxo pipefail
    zcat ~{vcf} | fgrep -e "#" > ~{prefix}.vcf
    zcat ~{vcf} | { fgrep -w -f <(cut -f4 ~{bed}) || true; } >> ~{prefix}.vcf
    bgzip ~{prefix}.vcf

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

task ConcatGenotypedVcfs {
  input {
    Array[File] lt5kb_vcfs
    Array[File] gt5kb_vcfs
    Array[File] bca_vcfs
    Array[File] ins_vcfs
    String batch
    String evidence_type    # depth or pesr
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File genotyped_vcf = "${batch}.~{evidence_type}.vcf.gz"
    File genotyped_vcf_index = "${batch}.~{evidence_type}.vcf.gz.tbi"
  }
  command <<<

    set -euo pipefail
    vcf-concat ~{sep=" " lt5kb_vcfs} ~{sep=" " gt5kb_vcfs} ~{sep=" " bca_vcfs} \
      | vcf-sort -c \
      | bgzip -c > ~{batch}.~{evidence_type}.vcf.gz
    tabix -p vcf ~{batch}.~{evidence_type}.vcf.gz

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

task MergePESRCounts {
  input {
    Array[File] count_list
    Array[File] sum_list
    String evidence_type    # pe or sr
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File counts = "~{evidence_type}_counts.txt.gz"
    File sum = "~{evidence_type}_sum.txt.gz"
  }
  command <<<

    set -euo pipefail
    zcat ~{sep=" " count_list} | fgrep -v -e "name" | gzip -c > ~{evidence_type}_counts.txt.gz
    echo "" | gzip -c > empty_file.gz
    cat ~{sep=" " sum_list} empty_file.gz > ~{evidence_type}_sum.txt.gz

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

task RDTestGenotype {
  input {
    File bed
    File bin_exclude
    File bin_exclude_idx
    File coveragefile
    File? coveragefile_index
    File medianfile
    File? famfile
    File ref_dict
    Array[String] samples
    File gt_cutoffs
    Int n_bins
    String prefix
    Boolean generate_melted_genotypes
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    coveragefile: {
      localization_optional: true
    }
    coveragefile_index: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.7)

  output {
    File genotypes = "${prefix}.geno"
    File copy_states = "${prefix}.median_geno"
    File metrics = "${prefix}.metrics"
    File gq = "${prefix}.gq"
    File varGQ = "${prefix}.vargq"
    File melted_genotypes = "rd.geno.cnv.bed.gz"
  }
  command <<<

    set -euo pipefail

    if [ -s ~{bed} ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{coveragefile} \
        -L ~{bed} \
        -O local.RD.txt.gz

      tabix -f -0 -s1 -b2 -e3 local.RD.txt.gz
    else
      touch local.RD.txt
      bgzip local.RD.txt
      tabix -p bed local.RD.txt.gz
    fi

    Rscript /opt/RdTest/RdTest.R \
      -b ~{bed} \
      -c local.RD.txt.gz \
      -m ~{medianfile} \
      ~{"-f " + famfile} \
      -n ~{prefix} \
      -w ~{write_lines(samples)} \
      -i ~{n_bins} \
      -r ~{gt_cutoffs} \
      -y ~{bin_exclude} \
      -g TRUE

    # In case of empty output, these files are not created
    touch ~{prefix}.geno
    touch ~{prefix}.median_geno
    touch ~{prefix}.metrics
    touch ~{prefix}.gq
    touch ~{prefix}.vargq

    if [ ~{generate_melted_genotypes} == "true" ]; then
      /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py ~{prefix}.geno ~{prefix}.gq rd.geno.cnv.bed
      sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz
    else
      echo "" | bgzip -c > rd.geno.cnv.bed.gz
    fi

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task CountPE {
  input {
    File vcf
    File discfile
    File? discfile_index
    File medianfile
    File ref_dict
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    discfile: {
      localization_optional: true
    }
    discfile_index: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.7)

  String prefix = basename(vcf, ".vcf")

  output {
    File pe_counts = "${prefix}.pe_counts.txt.gz"
  }
  command <<<

    set -euo pipefail
    svtk vcf2bed --split-bnd --no-header ~{vcf} test.bed
    awk -v OFS="\t" -v window=500 '{if ($2-window>0){print $1,$2-window,$2+window}else{print $1,0,$2+window}}' test.bed  > region.bed
    awk -v OFS="\t" -v window=500 '{if ($3-window>0){print $1,$3-window,$3+window}else{print $1,0,$3+window}}' test.bed  >> region.bed
    sort -k1,1 -k2,2n region.bed > region.sorted.bed
    bedtools merge -i region.sorted.bed > region.merged.bed

    if [ -s region.merged.bed ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{discfile} \
        -L region.merged.bed \
        -O local.PE.txt.gz

      tabix -f -0 -s1 -b2 -e2 local.PE.txt.gz
    else
      touch local.PE.txt
      bgzip local.PE.txt
      tabix -0 -s1 -b2 -e2 local.PE.txt.gz
    fi

    svtk count-pe -s ~{write_lines(samples)} --medianfile ~{medianfile} ~{vcf} local.PE.txt.gz ~{prefix}.pe_counts.txt
    gzip ~{prefix}.pe_counts.txt

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}

task CountSR {
  input {
    File vcf
    File splitfile
    File? splitfile_index
    File medianfile
    File ref_dict
    Array[String] samples
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    splitfile: {
      localization_optional: true
    }
    splitfile_index: {
      localization_optional: true
    }
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  Float mem_gb = select_first([runtime_attr.mem_gb, default_attr.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.7)

  String prefix = basename(vcf, ".vcf")

  output {
    File sr_counts = "${prefix}.sr_counts.txt.gz"
    File sr_sum = "${prefix}.sr_sum.txt.gz"
  }
  command <<<

    set -euo pipefail
    svtk vcf2bed --split-bnd --no-header ~{vcf} test.bed
    awk -v OFS="\t" '{if ($2-250>0){print $1,$2-250,$2+250}else{print $1,0,$2+250}}' test.bed  >> region.bed
    awk -v OFS="\t" '{if ($3-250>0){print $1,$3-250,$3+250}else{print $1,0,$3+250}}' test.bed  >> region.bed
    sort -k1,1 -k2,2n region.bed > region.sorted.bed
    bedtools merge -i region.sorted.bed > region.merged.bed

    if [ -s region.merged.bed ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{splitfile} \
        -L region.merged.bed \
        -O local.SR.txt.gz

      tabix -f -0 -s1 -b2 -e2 local.SR.txt.gz
    else
      touch local.SR.txt
      bgzip local.SR.txt
      tabix -0 -s1 -b2 -e2 local.SR.txt.gz
    fi

    svtk count-sr -s ~{write_lines(samples)} --medianfile ~{medianfile} ~{vcf} local.SR.txt.gz ~{prefix}.sr_counts.txt
    /opt/sv-pipeline/04_variant_resolution/scripts/sum_SR.sh ~{prefix}.sr_counts.txt ~{prefix}.sr_sum.txt.gz
    gzip ~{prefix}.sr_counts.txt

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}
task AddBatchSamples {
  input {
    File batch_vcf
    File cohort_vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File updated_vcf = "${prefix}.vcf.gz"
  }
  command <<<

    set -euo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/add_batch_samples.py ~{batch_vcf} ~{cohort_vcf} ~{prefix}.vcf
    bgzip ~{prefix}.vcf
  
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
task IntegrateDepthGq {
  input {
    File vcf
    File RD_melted_genotypes
    File RD_vargq
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1, 
    mem_gb: 3.75,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File genotypes = "genotype.indiv.depth.txt.gz"
    File varGQ = "genotype.variant.depth.txt.gz"
  }
  command <<<

    /opt/sv-pipeline/04_variant_resolution/scripts/IntegrateGQ_depthonly.sh \
      ~{vcf} \
      ~{RD_melted_genotypes} \
      ~{RD_vargq}
  
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

task ReformatGenotypedVcf {
  input {
    File vcf
    File? script
    String output_prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + size(vcf, "GB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File out = "~{output_prefix}.vcf.gz"
    File out_index = "~{output_prefix}.vcf.gz.tbi"
  }
  command <<<
    set -euo pipefail
    python ~{default="/opt/sv-pipeline/04_variant_resolution/scripts/reformat_genotyped_vcf.py" script} --vcf ~{vcf} \
      | bgzip \
      > ~{output_prefix}.vcf.gz
    tabix ~{output_prefix}.vcf.gz

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