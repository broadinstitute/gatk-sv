version 1.0

import "Structs.wdl"

task FilterVcfBySampleGenotypeAndAddEvidenceAnnotation {
  input {
    File vcf_gz
    String sample_id
    String sv_base_mini_docker
    String evidence
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

  String filebase = basename(vcf_gz, ".vcf.gz")
  String outfile = "~{filebase}.~{sample_id}.vcf.gz"
  Array[String] sample_array = [sample_id]

  output {
    File out = "~{outfile}"
  }
  command <<<
    set -euo pipefail
    sampleIndex=`gzip -cd ~{vcf_gz} | grep '^#CHROM' | cut -f10- | tr "\t" "\n" | awk '$1 == "~{sample_id}" {found=1; print NR - 1} END { if (found != 1) { print "sample not found"; exit 1; }}'`

    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t~{evidence}\n" ~{vcf_gz} | bgzip -c > evidence_annotations.tab.gz
    tabix -s1 -b2 -e2 evidence_annotations.tab.gz

    echo '##INFO=<ID=EVIDENCE,Number=.,Type=String,Description="Classes of random forest support.">' > header_line.txt

    bcftools annotate \
        -i "GT[${sampleIndex}]=\"alt\"" \
        -a evidence_annotations.tab.gz \
        -c CHROM,POS,REF,ALT,EVIDENCE \
        -h header_line.txt \
        ~{vcf_gz} | bgzip -c > ~{outfile}
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


task FilterVcfForShortDepthCalls {
  input {
    File vcf_gz
    Int min_length
    String filter_name

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

  String filebase = basename(vcf_gz, ".vcf.gz")
  String outfile = "~{filebase}.filter_depth_lt_~{min_length}.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<
    set -euo pipefail

    bcftools filter ~{vcf_gz} \
      -i 'INFO/SVLEN >= ~{min_length} || INFO/ALGORITHMS[*] != "depth"' \
      -s ~{filter_name} |
         bgzip -c > ~{outfile}

    tabix ~{outfile}

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


task GetUniqueNonGenotypedDepthCalls {
  input {
    File vcf_gz
    String sample_id
    File ref_panel_dels
    File ref_panel_dups
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

  String filebase = basename(vcf_gz, ".vcf.gz")
  String outfile = "~{filebase}.unique_~{sample_id}_depth_non_alt_gt.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"

  }
  command <<<
    set -euo pipefail
    sampleIndex=`gzip -cd ~{vcf_gz} | grep '^#CHROM' | cut -f10- | tr "\t" "\n" | awk '$1 == "~{sample_id}" {found=1; print NR - 1} END { if (found != 1) { print "sample not found"; exit 1; }}'`

    bcftools filter \
        -i "FILTER == \"PASS\" && GT[${sampleIndex}]!=\"alt\" && INFO/ALGORITHMS[*] == \"depth\" && INFO/SVLEN >= 10000 && (INFO/SVTYPE == \"DEL\" || INFO/SVTYPE == \"DUP\")" \
        ~{vcf_gz} | bgzip -c > pass_depth_not_alt.vcf.gz

    bgzip -cd pass_depth_not_alt.vcf.gz | grep '#' > header.txt

    zcat ~{ref_panel_dels} ~{ref_panel_dups} | sort -k1,1V -k2,2n > ref_panel_depth_calls.bed
    intersectBed -a pass_depth_not_alt.vcf.gz -b ref_panel_depth_calls.bed -r -f .5 -v > unique_depth_records.txt

    cat header.txt unique_depth_records.txt | bgzip -c > ~{outfile}

    tabix ~{outfile}
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

task FilterVcfForCaseSampleGenotype {
  input {
    File vcf_gz
    String sample_id
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

  String filebase = basename(vcf_gz, ".vcf.gz")
  String outfile = "~{filebase}.filter_by_~{sample_id}_gt.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<
    set -euo pipefail
    sampleIndex=`gzip -cd ~{vcf_gz} | grep '^#CHROM' | cut -f10- | tr "\t" "\n" | awk '$1 == "~{sample_id}" {found=1; print NR - 1} END { if (found != 1) { print "sample not found"; exit 1; }}'`

    bcftools filter \
        -i "FILTER ~ \"MULTIALLELIC\" || GT[${sampleIndex}]=\"alt\"" \
        ~{vcf_gz} | bgzip -c > ~{outfile}

    tabix ~{outfile}
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

task FilterLargePESRCallsWithoutRawDepthSupport {
  input {
    File pesr_vcf
    File raw_dels
    File raw_dups

    Int? min_large_pesr_call_size_for_filtering
    Float? min_large_pesr_depth_overlap_fraction

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

  String filebase = basename(pesr_vcf, ".vcf.gz")
  String outfile = "~{filebase}.filter_large_pesr_by_depth.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  
  command <<<
    set -euo pipefail

    svtk vcf2bed ~{pesr_vcf} stdout | cut -f1-5 | awk '$3 - $2 > ~{if defined(min_large_pesr_call_size_for_filtering) then min_large_pesr_call_size_for_filtering else 1000000} && ($5 == "DEL")' \
        | coverageBed -a stdin -b ~{raw_dels} | awk '$NF < ~{if defined(min_large_pesr_depth_overlap_fraction) then min_large_pesr_depth_overlap_fraction else "0.3"} {print $4}' > large_dels_without_raw_depth_support.list

    svtk vcf2bed ~{pesr_vcf} stdout | cut -f1-5 | awk '$3 - $2 > ~{if defined(min_large_pesr_call_size_for_filtering) then min_large_pesr_call_size_for_filtering else 1000000} && ($5 == "DUP")' \
        | coverageBed -a stdin -b ~{raw_dups} | awk '$NF < ~{if defined(min_large_pesr_depth_overlap_fraction) then min_large_pesr_depth_overlap_fraction else "0.3"} {print $4}' > large_dups_without_raw_depth_support.list

    cat large_dels_without_raw_depth_support.list large_dups_without_raw_depth_support.list > large_pesr_without_raw_depth_support.list

    cat \
        <(gzip -cd ~{pesr_vcf} | grep -v '^#' | grep -w -f large_pesr_without_raw_depth_support.list | sed -e 's/SVTYPE=DEL/SVTYPE=BND/' -e 's/SVTYPE=DUP/SVTYPE=BND/' -e 's/<DEL>/<BND>/' -e 's/<DUP>/<BND>/') \
        <(gzip -cd ~{pesr_vcf} | grep -v '^#' | grep -v -w -f large_pesr_without_raw_depth_support.list) \
       | cat <(sed -n -e '/^#/p' <(zcat ~{pesr_vcf})) - \
       | vcf-sort -c \
       | bgzip -c \
       > ~{outfile}

    tabix ~{outfile}
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

task FilterVcfWithReferencePanelCalls {
  input {
    File single_sample_vcf
    File cohort_vcf
    String case_sample_id
    Float? max_ref_panel_carrier_freq
    Float? required_reciprocal_overlap_match_pct
    Float? required_cnv_coverage_pct

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

  String filebase = basename(single_sample_vcf, ".vcf.gz")
  String outfile = "~{filebase}.filter_by_ref_panel.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<
  set -euo pipefail
  BCFTOOLS=/usr/local/bin/bcftools
  $BCFTOOLS query -l ~{cohort_vcf} > ref_samples.list
  maxRefSamples=$(wc -l ref_samples.list | awk '{print sprintf("%.0f", $1 * ~{default="0.01" max_ref_panel_carrier_freq})}')

  /opt/sv-pipeline/scripts/single_sample/apply_ref_panel_genotypes_filter.py \
    ~{single_sample_vcf} \
    ~{cohort_vcf} \
    ~{case_sample_id} ~{default="0.5" required_cnv_coverage_pct} ${maxRefSamples} \
    del_dup_cpx_inv_filtered.vcf.gz

  $BCFTOOLS query \
      -i 'FILTER ~ "MULTIALLELIC"' \
      -f '%CHROM\t%POS\t%END\t%ID\t%SVLEN\n' \
      ~{cohort_vcf} | \
      awk '{OFS="\t"; $3 = $2 + $5; print}' \
      > cohort.gts.multiallelic.bed

  $BCFTOOLS query \
      -i 'FILTER ~ "MULTIALLELIC"' \
      -S ref_samples.list \
      -f '%CHROM\t%POS\t%END\t%ID\t%SVLEN\n' \
      ~{single_sample_vcf} | \
      awk '{OFS="\t"; $3 = $2 + $5; print}' \
      > case.ref_panel_variant.multiallelic.bed

  intersectBed \
      -a case.ref_panel_variant.multiallelic.bed \
      -b cohort.gts.multiallelic.bed \
      -f ~{default="0.5" required_reciprocal_overlap_match_pct} -r -v -wa | cut -f4 > multiallelics.list

  $BCFTOOLS query -i 'SVTYPE="INS"' \
      -f '%CHROM\t%POS\t%END\t%ID\t%SVLEN\n' \
      ~{cohort_vcf} | \
      awk '{OFS="\t"; $3 = $2 + $5; print}' > cohort.gts.ins.bed

  $BCFTOOLS query -i "SVTYPE == 'INS' && AC >= ${maxRefSamples}" \
      -S ref_samples.list \
      -f '%CHROM\t%POS\t%END\t%ID\t%SVLEN\n' \
      ~{single_sample_vcf} | \
      awk '{OFS="\t"; $3 = $2 + $5; print}'  > case.ref_panel_variant.ins.bed

  intersectBed \
      -a case.ref_panel_variant.ins.bed \
      -b cohort.gts.ins.bed \
      -f ~{default="0.5" required_reciprocal_overlap_match_pct} -r -v -wa | cut -f4 > case_variants_not_in_ref_panel.ins.list

  $BCFTOOLS query -i 'SVTYPE="BND"' \
      -f '%CHROM\t%POS\t%INFO/CHR2\t%INFO/END\t%ID\t\n' \
      ~{cohort_vcf} | \
      awk '{OFS="\t"; print $1,$2-50,$2+50,$3,$4-50,$4+50,$5}' > cohort.gts.bnd.bedpe

  $BCFTOOLS query -i "SVTYPE='BND' && AC >= ${maxRefSamples}" \
      -S ref_samples.list \
      -f '%CHROM\t%POS\t%INFO/CHR2\t%INFO/END\t%ID\t\n' \
      ~{single_sample_vcf} | \
      awk '{OFS="\t"; print $1,$2-50,$2+50,$3,$4-50,$4+50,$5}' > case.gts.bnd.bedpe

  pairToPair -a case.gts.bnd.bedpe \
      -b cohort.gts.bnd.bedpe \
      -type both |\
      cut -f7 | sort -u > case_bnds_to_keep.list

  cp case_variants_not_in_ref_panel.ins.list case_variants_not_in_ref_panel.list

  $BCFTOOLS filter \
      -e 'ID=@case_variants_not_in_ref_panel.list || ( SVTYPE="BND" && ID!=@case_bnds_to_keep.list ) || (FILTER ~ "MULTIALLELIC" && ID!=@multiallelics.list )' \
      -s REF_PANEL_GENOTYPES \
      -m + \
      del_dup_cpx_inv_filtered.vcf.gz | bgzip -c > ~{outfile}
  tabix ~{outfile}
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

task ResetFilter {
  input {
    File single_sample_vcf
    File single_sample_vcf_idx
    String filter_to_reset
    String info_header_line

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

  String filebase = basename(single_sample_vcf, ".vcf.gz")
  String outfile = "~{filebase}.reset_~{filter_to_reset}_filter.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<

  set -euo pipefail

  echo '~{info_header_line}' > header.txt

  bcftools filter -i 'FILTER ~ "~{filter_to_reset}"' ~{single_sample_vcf} | bgzip -c > sites.vcf.gz
  tabix sites.vcf.gz

  bcftools annotate \
    -k \
    -a sites.vcf.gz \
    -m ~{filter_to_reset} \
    -h header.txt \
    -x FILTER/~{filter_to_reset} \
    ~{single_sample_vcf} | bgzip -c > ~{outfile}

  tabix ~{outfile}
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

task FinalVCFCleanup {
  input {
    File single_sample_vcf
    File single_sample_vcf_idx
    File ref_fasta
    File ref_fasta_idx

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

  String filebase = basename(single_sample_vcf, ".vcf.gz")
  String outfile = "~{filebase}.final_cleanup.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }
  command <<<

     set -euo pipefail

     # clean up bad END tags
     bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/END\n' \
        ~{single_sample_vcf} | awk '$5 < $2 {OFS="\t"; $5 = $2 + 1; print}' | bgzip -c > bad_ends.txt.gz
     tabix -s1 -b2 -e2 bad_ends.txt.gz
     bcftools annotate \
        -a bad_ends.txt.gz \
        -c CHROM,POS,REF,ALT,END \
        ~{single_sample_vcf} | bgzip -c > fix_bad_ends.vcf.gz

     /opt/sv-pipeline/scripts/single_sample/update_variant_representations.py fix_bad_ends.vcf.gz ~{ref_fasta} \
        | vcf-sort -c \
        | bgzip -c > ~{outfile}

     tabix ~{outfile}
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

task RewriteSRCoords {
  input {
    File vcf
    File metrics
    File cutoffs
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 10,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File annotated_vcf = "${prefix}.corrected_coords.vcf.gz"
  }
  command <<<

    set -euo pipefail

    /opt/sv-pipeline/03_variant_filtering/scripts/rewrite_SR_coords.py ~{vcf} ~{metrics} ~{cutoffs} stdout \
      | vcf-sort -c \
      | bgzip -c \
      > ~{prefix}.corrected_coords.vcf.gz

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

task ConvertCNVsWithoutDepthSupportToBNDs {
  input {
    File genotyped_pesr_vcf
    File allosome_file
    String case_sample
    File merged_famfile
    Int? min_length
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 3.5,
    disk_gb: 10,
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  String filebase = basename(genotyped_pesr_vcf, ".vcf.gz")
  String outfile = "~{filebase}.convert_cnvs_to_bnd.vcf.gz"

  output {
    File out_vcf = "~{outfile}"
    File out_vcf_idx = "~{outfile}.tbi"
  }
  command <<<

    set -euo pipefail

    /opt/sv-pipeline/scripts/single_sample/convert_cnvs_without_depth_support_to_bnds.py \
        ~{genotyped_pesr_vcf} \
        ~{allosome_file} \
        ~{merged_famfile} \
        ~{case_sample} \
        ~{default="1000" min_length} \
        -o ~{outfile}

    tabix ~{outfile}

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

task VcfToBed {
  input {
    File vcf
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
    File bed = "${prefix}.bed"
  }
  command <<<

    svtk vcf2bed ~{vcf} ~{prefix}.bed -i ALL --include-filters

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

task SampleQC {
  input {
    File vcf
    File sample_filtering_qc_file
    String sv_pipeline_base_docker
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

  String filebase = basename(vcf, ".vcf.gz")
  String outfile = "~{filebase}.sample_qc.vcf.gz"

  output {
    File out = "~{outfile}"
    File out_idx = "~{outfile}.tbi"
  }

  command <<<
    set -euo pipefail

    wgdPF=`cat ~{sample_filtering_qc_file} | awk '$1 == "wgd_score_sample" {print $5}'`
    if [ $wgdPF = "FAIL" ]
    then
      bcftools filter -e 'SVTYPE!="BND"' -m + -s SAMPLE_WGD_OUTLIER ~{vcf} \
        | sed 's/ID=SAMPLE_WGD_OUTLIER,Description=".*"/ID=SAMPLE_WGD_OUTLIER,Description="Case sample is an outlier for WGD dosage score compared to reference panel"/' \
        | bgzip -c \
        > wgd_filtered.vcf.gz
    else
      cp ~{vcf} wgd_filtered.vcf.gz
    fi

    coveragePF=`cat ~{sample_filtering_qc_file} | awk '$1 == "rd_mean_sample" {print $5}'`
    if [ $coveragePF = "FAIL" ]
    then
      bcftools filter -e 'SVTYPE!="BND"' -m + -s SAMPLE_COVERAGE_OUTLIER wgd_filtered.vcf.gz \
        | sed 's/ID=SAMPLE_COVERAGE_OUTLIER,Description=".*"/ID=SAMPLE_COVERAGE_OUTLIER,Description="Case sample is an outlier for coverage compared to reference panel"/' \
        | bgzip -c \
      > ~{outfile}
    else
      cp wgd_filtered.vcf.gz ~{outfile}
    fi

    tabix ~{outfile}

  >>>
  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_base_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }

}


