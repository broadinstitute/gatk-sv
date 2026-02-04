version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "FormatVcfForGatk.wdl" as fvcf

workflow CleanVcfChromosome {
  input {
    File vcf
    String contig
    File background_list
    File bothsides_pass_list
    File? outlier_samples_list
    File ped_file
    File ploidy_table

    String chr_x
    String chr_y
    String prefix
    Int format_vcf_records_per_shard = 5000
    Int preprocess_records_per_shard = 5000
    Int postprocess_records_per_shard = 5000

    File contig_list
    File allosome_fai
    
    File HERVK_reference
    File LINE1_reference
    File intron_reference

    String gatk_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_format_to_clean_create_ploidy
    RuntimeAttr? runtime_attr_format_to_clean_scatter
    RuntimeAttr? runtime_attr_format_to_clean_format
    RuntimeAttr? runtime_attr_format_to_clean_concat
    RuntimeAttr? runtime_attr_scatter_preprocess
    RuntimeAttr? runtime_attr_preprocess
    RuntimeAttr? runtime_attr_concat_preprocess
    RuntimeAttr? runtime_attr_revise_overlapping_cnvs
    RuntimeAttr? runtime_attr_revise_large_cnvs
    RuntimeAttr? runtime_attr_revise_multiallelics
    RuntimeAttr? runtime_attr_scatter_postprocess
    RuntimeAttr? runtime_attr_postprocess
    RuntimeAttr? runtime_attr_concat_postprocess
    RuntimeAttr? runtime_override_drop_redundant_cnvs
    RuntimeAttr? runtime_override_sort_drop_redundant_cnvs
    RuntimeAttr? runtime_override_stitch_fragmented_cnvs
    RuntimeAttr? runtime_override_rescue_me_dels
    RuntimeAttr? runtime_attr_add_high_fp_rate_filters
    RuntimeAttr? runtime_attr_add_retro_del_filters
    RuntimeAttr? runtime_override_final_cleanup
    RuntimeAttr? runtime_attr_format_to_output_create_ploidy
    RuntimeAttr? runtime_attr_format_to_output_scatter
    RuntimeAttr? runtime_attr_format_to_output_format
    RuntimeAttr? runtime_attr_format_to_output_concat
  }

  call fvcf.FormatVcfForGatk as FormatVcfToClean {
    input:
      vcf=vcf,
      prefix="~{prefix}.formatted",
      ped_file=ped_file,
      records_per_shard=format_vcf_records_per_shard,
      contig_list=contig_list,
      bothside_pass_list=bothsides_pass_list,
      background_fail_list=background_list,
      chr_x=chr_x,
      chr_y=chr_y,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_create_ploidy=runtime_attr_format_to_clean_create_ploidy,
      runtime_attr_scatter=runtime_attr_format_to_clean_scatter,
      runtime_attr_format=runtime_attr_format_to_clean_format,
      runtime_override_concat=runtime_attr_format_to_clean_concat
  }

  call MiniTasks.ScatterVcf as ScatterPreprocess {
    input:
      vcf=FormatVcfToClean.gatk_formatted_vcf,
      vcf_index=FormatVcfToClean.gatk_formatted_vcf_index,
      prefix="~{prefix}.preprocess.scatter",
      records_per_shard=preprocess_records_per_shard,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_scatter_preprocess
  }

  scatter (shard in ScatterPreprocess.shards) {
    call CleanVcfPreprocess {
      input:
        vcf=shard,
        chr_x=chr_x,
        chr_y=chr_y,
        background_list=background_list,
        bothsides_pass_list=bothsides_pass_list,
        ped_file=ped_file,
        prefix="~{prefix}.preprocess",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_preprocess
    }
  }

  call MiniTasks.ConcatVcfs as ConcatPreprocess {
    input:
      vcfs=CleanVcfPreprocess.out,
      vcfs_idx=CleanVcfPreprocess.out_idx,
      allow_overlaps=true,
      outfile_prefix="~{prefix}.preprocess.concat",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_preprocess
  }

  call CleanVcfReviseOverlappingCnvs {
    input:
      vcf=ConcatPreprocess.concat_vcf,
      vcf_idx=ConcatPreprocess.concat_vcf_idx,
      prefix="~{prefix}.revise_overlapping_cnvs",
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_revise_overlapping_cnvs
  }

  call CleanVcfReviseMultiallelicCnvs {
    input:
      vcf=CleanVcfReviseOverlappingCnvs.out,
      vcf_idx=CleanVcfReviseOverlappingCnvs.out_idx,
      outlier_samples_list=outlier_samples_list,
      prefix="~{prefix}.revise_multiallelic_cnvs",
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_revise_large_cnvs
  }

  call CleanVcfReviseOverlappingMultiallelics {
    input:
      vcf=CleanVcfReviseMultiallelicCnvs.out,
      vcf_idx=CleanVcfReviseMultiallelicCnvs.out_idx,
      prefix="~{prefix}.revise_overlapping_multiallelics",
      gatk_docker=gatk_docker,
      runtime_attr_override=runtime_attr_revise_multiallelics
  }

  call MiniTasks.ScatterVcf as ScatterPostprocess {
    input:
      vcf=CleanVcfReviseOverlappingMultiallelics.out,
      vcf_index=CleanVcfReviseOverlappingMultiallelics.out_idx,
      prefix="~{prefix}.postprocess.scatter",
      records_per_shard=postprocess_records_per_shard,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_scatter_postprocess
  }

  scatter (shard in ScatterPostprocess.shards) {
    call CleanVcfPostprocess {
      input:
        vcf=shard,
        ped_file=ped_file,
        prefix="~{prefix}.postprocess",
        sv_pipeline_docker=sv_pipeline_docker,
        runtime_attr_override=runtime_attr_postprocess
    }
  }

  call MiniTasks.ConcatVcfs as ConcatPostprocess {
    input:
      vcfs=CleanVcfPostprocess.out,
      vcfs_idx=CleanVcfPostprocess.out_idx,
      allow_overlaps=true,
      outfile_prefix="~{prefix}.postprocess.concat",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_attr_concat_postprocess
  }

  call DropRedundantCnvs {
    input:
      vcf=ConcatPostprocess.concat_vcf,
      prefix="~{prefix}.drop_redundant_cnvs",
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_drop_redundant_cnvs
  }

  call MiniTasks.SortVcf as SortDropRedundantCnvs {
    input:
      vcf=DropRedundantCnvs.out,
      outfile_prefix="~{prefix}.drop_redundant_cnvs.sorted",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_sort_drop_redundant_cnvs
  }

  call StitchFragmentedCnvs {
    input:
      vcf=SortDropRedundantCnvs.out,
      prefix="~{prefix}.stitch_fragmented_cnvs",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_stitch_fragmented_cnvs
  }

  call RescueMobileElementDeletions {
    input:
      vcf = StitchFragmentedCnvs.stitched_vcf_shard,
      prefix = "~{prefix}.rescue_me_dels",
      LINE1 = LINE1_reference,
      HERVK = HERVK_reference,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_override_rescue_me_dels
  }

  call AddHighFDRFilters {
    input:
      vcf=RescueMobileElementDeletions.out,
      prefix="~{prefix}.high_fdr_filtered",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_add_high_fp_rate_filters
  }

  call AddRetroDelFilters {
    input:
      vcf=AddHighFDRFilters.out,
      intron_reference=intron_reference,
      contig=contig,
      prefix="~{prefix}.retro_del_filtered",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_attr_add_retro_del_filters
  }

  call FinalCleanup {
    input:
      vcf=AddRetroDelFilters.out,
      contig=contig,
      prefix="~{prefix}.final_cleanup",
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_final_cleanup
  }

  call fvcf.FormatVcfForGatk as FormatVcfToOutput {
    input:
      vcf=FinalCleanup.final_cleaned_shard,
      prefix="~{prefix}.final_format",
      ped_file=ped_file,
      records_per_shard=format_vcf_records_per_shard,
      contig_list=contig_list,
      bothside_pass_list=bothsides_pass_list,
      background_fail_list=background_list,
      chr_x=chr_x,
      chr_y=chr_y,
      sv_base_mini_docker=sv_base_mini_docker,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_create_ploidy=runtime_attr_format_to_output_create_ploidy,
      runtime_attr_scatter=runtime_attr_format_to_output_scatter,
      runtime_attr_format=runtime_attr_format_to_output_format,
      runtime_override_concat=runtime_attr_format_to_output_concat
  }
    
  output {
    File out = FormatVcfToOutput.gatk_formatted_vcf
    File out_idx = FormatVcfToOutput.gatk_formatted_vcf_index
  }
}

task CleanVcfPreprocess {
  input {
    File vcf
    String chr_x
    String chr_y
    File background_list
    File bothsides_pass_list
    File ped_file
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 4,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
  String output_vcf = "~{prefix}.vcf.gz"

  command <<<
    set -euxo pipefail

    bcftools view --header-only ~{vcf} | grep '^##' > header.txt

    cat <<EOF >> header.txt
    ##FILTER=<ID=UNRESOLVED,Description="Variant is unresolved">
    ##INFO=<ID=HIGH_SR_BACKGROUND,Number=0,Type=Flag,Description="Variant has high number of SR splits in background samples">
    ##INFO=<ID=BOTHSIDES_SUPPORT,Number=0,Type=Flag,Description="Variant has read-level support for both sides of breakpoint">
    ##INFO=<ID=REVISED_EVENT,Number=0,Type=Flag,Description="Variant has been revised due to a copy number mismatch">
    EOF

    bcftools view --header-only ~{vcf} | grep '^#CHROM' >> header.txt

    bcftools view ~{vcf} | bcftools reheader -h header.txt | bgzip -c > processed.reheader.vcf.gz

    rm header.txt
    
    python /opt/sv-pipeline/04_variant_resolution/scripts/cleanvcf_preprocess.py \
      -V processed.reheader.vcf.gz \
      -O ~{output_vcf} \
      --chrX ~{chr_x} \
      --chrY ~{chr_y} \
      --fail-list ~{background_list} \
      --pass-list ~{bothsides_pass_list} \
      --ped-file ~{ped_file}

    tabix -p vcf ~{output_vcf}
  >>>

  output {
    File out="~{output_vcf}"
    File out_idx="~{output_vcf}.tbi"
  }
}

task CleanVcfReviseOverlappingCnvs {
  input {
    File vcf
    File vcf_idx
    String prefix
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 4,
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
    docker: gatk_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
  String output_vcf = "~{prefix}.vcf.gz"

  command <<<
    set -euo pipefail
    
    gatk --java-options "-Xmx~{java_mem_mb}m" SVReviseOverlappingCnvs \
      -V ~{vcf} \
      -O ~{output_vcf}
  >>>

  output {
    File out="~{output_vcf}"
    File out_idx="~{output_vcf}.tbi"
  }
}

task CleanVcfReviseMultiallelicCnvs {
  input {
    File vcf
    File vcf_idx
    File? outlier_samples_list
    String prefix
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 4,
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
    docker: gatk_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
  String output_vcf = "~{prefix}.vcf.gz"

  command <<<
    set -euo pipefail
    
    gatk --java-options "-Xmx~{java_mem_mb}m" SVReviseMultiallelicCnvs \
      -V ~{vcf} \
      -O ~{output_vcf} \
      ~{if defined(outlier_samples_list) then "--outlier-samples ~{outlier_samples_list}" else "" }
  >>>

  output {
    File out="~{output_vcf}"
    File out_idx="~{output_vcf}.tbi"
  }
}

task CleanVcfReviseOverlappingMultiallelics {
  input {
    File vcf
    File vcf_idx
    String prefix
    String gatk_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 4,
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
    docker: gatk_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
  String output_vcf = "~{prefix}.vcf.gz"

  command <<<
    set -euo pipefail
    
    gatk --java-options "-Xmx~{java_mem_mb}m" SVReviseOverlappingMultiallelics \
      -V ~{vcf} \
      -O ~{output_vcf}
  >>>

  output {
    File out="~{output_vcf}"
    File out_idx="~{output_vcf}.tbi"
  }
}

task CleanVcfPostprocess {
  input {
    File vcf
    File? vcf_idx
    File ped_file
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 4,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)
  String output_vcf = "~{prefix}.vcf.gz"

  command <<<
    set -euo pipefail

    if [ ! -f "~{vcf}.tbi" ]; then
      tabix -p vcf ~{vcf}
    fi

    python /opt/sv-pipeline/04_variant_resolution/scripts/cleanvcf_postprocess.py \
      -V ~{vcf} \
      -O processed.vcf.gz \
      --ped-file ~{ped_file}

    bcftools annotate -x INFO/MULTIALLELIC,INFO/UNRESOLVED,INFO/EVENT,INFO/REVISED_EVENT,INFO/MULTI_CNV,INFO/varGQ processed.vcf.gz -o processed.annotated.vcf.gz -O z

    bcftools view -h processed.annotated.vcf.gz | grep "^##" | \
      grep -v -E "CIPOS|CIEND|RMSSTD|source|bcftools|GATKCommandLine|##FORMAT=<ID=EV>|##ALT=<ID=UNR>|##INFO=<ID=(MULTIALLELIC|UNRESOLVED|EVENT|REVISED_EVENT|MULTI_CNV|varGQ)" > temp_header.txt
    echo '##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description="Class of unresolved variant.">' >> temp_header.txt
    echo '##ALT=<ID=CNV,Description="Copy Number Polymorphism">' >> temp_header.txt

    bcftools view -h processed.annotated.vcf.gz | grep "^#CHROM" > chrom_header.txt

    cat temp_header.txt chrom_header.txt > header.txt
    
    bcftools reheader -h header.txt processed.annotated.vcf.gz -o ~{output_vcf}
    
    tabix -p vcf ~{output_vcf}
  >>>

  output {
    File out="~{output_vcf}"
    File out_idx="~{output_vcf}.tbi"
  }
}

task RescueMobileElementDeletions {
  input {
    File vcf
    String prefix
    File LINE1
    File HERVK
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 2,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python <<CODE
import os
import pysam
fin=pysam.VariantFile("~{vcf}")
fo=pysam.VariantFile("~{prefix}.bnd_del.vcf.gz", 'w', header = fin.header)
for record in fin:
    if record.info['SVTYPE'] in ['BND'] and record.info['STRANDS']=="+-" and record.chrom == record.info['CHR2'] and record.info['END2'] - record.start < 10000:
        record.info['SVLEN'] = record.info['END2'] - record.start
        fo.write(record)
fin.close()
fo.close()
CODE

    tabix -p vcf ~{prefix}.bnd_del.vcf.gz

    svtk vcf2bed ~{prefix}.bnd_del.vcf.gz -i ALL --include-filters ~{prefix}.bnd_del.bed
    bgzip ~{prefix}.bnd_del.bed

    bedtools coverage -wo -a ~{prefix}.bnd_del.bed.gz -b ~{LINE1} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_LINE1/' > manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv
    bedtools coverage -wo -a ~{prefix}.bnd_del.bed.gz -b ~{HERVK} | awk '{if ($NF>.5) print}' | cut -f4 | sed -e 's/$/\tDEL\tPASS\toverlap_HERVK/' >> manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv

    python <<CODE
import pysam
def SVID_MEI_DEL_readin(MEI_DEL_reset):
    out={}
    fin=open(MEI_DEL_reset)
    for line in fin:
        pin=line.strip().split()
        if not pin[0] in out.keys():
            out[pin[0]] = pin[3]
    fin.close()
    return out

hash_MEI_DEL_reset = SVID_MEI_DEL_readin("manual_revise.MEI_DEL_from_BND.SVID_SVTYPE_FILTER_INFO.tsv")
fin=pysam.VariantFile("~{vcf}")
fo=pysam.VariantFile("~{prefix}.vcf.gz", 'w', header = fin.header)
for record in fin:
    if record.id in hash_MEI_DEL_reset.keys():
        del record.filter['UNRESOLVED']
        record.info['SVTYPE'] = 'DEL'
        record.info['SVLEN'] = record.info['END2'] - record.start
        record.stop = record.info['END2']
        record.info.pop("CHR2")
        record.info.pop("END2")
        record.info.pop("UNRESOLVED_TYPE")
        if hash_MEI_DEL_reset[record.id] == 'overlap_LINE1':
            record.alts = ('<DEL:ME:LINE1>',)
        if hash_MEI_DEL_reset[record.id] == 'overlap_HERVK':
            record.alts = ('<DEL:ME:HERVK>',)
    fo.write(record)
fin.close()
fo.close()
CODE
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}


# Remove CNVs that are redundant with CPX events or other CNVs
task DropRedundantCnvs {
  input {
    File vcf
    String prefix
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 2,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    /opt/sv-pipeline/04_variant_resolution/scripts/resolve_cpx_cnv_redundancies.py \
      ~{vcf} ~{prefix}.vcf.gz --temp-dir ./tmp
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}


# Stitch fragmented RD-only calls found in 100% of the same samples
task StitchFragmentedCnvs {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 2,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  Int java_mem_mb = ceil(select_first([runtime_override.mem_gb, runtime_default.mem_gb]) * 1000 * 0.7)

  command <<<
    set -euo pipefail
    echo "First pass..."
    java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 ~{vcf} \
      | bgzip \
      > tmp.vcf.gz
    rm ~{vcf}
    echo "Second pass..."
    java -Xmx~{java_mem_mb}M -jar ${STITCH_JAR} 0.2 200000 0.2 tmp.vcf.gz \
      | bgzip \
      > ~{prefix}.vcf.gz
  >>>

  output {
    File stitched_vcf_shard = "~{prefix}.vcf.gz"
  }
}

# Add FILTER status for pockets of variants with high FP rate: wham-only DELs and Scramble-only SVAs with HIGH_SR_BACKGROUND
task AddHighFDRFilters {
  input {
    File vcf
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 2,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python <<CODE
import pysam
with pysam.VariantFile("~{vcf}", 'r') as fin:
  header = fin.header
  header.add_line("##FILTER=<ID=HIGH_ALGORITHM_FDR,Description=\"Categories of variants with low precision including Wham-only deletions and certain Scramble SVAs\">")
  with pysam.VariantFile("~{prefix}.vcf.gz", 'w', header=header) as fo:
    for record in fin:
        if (record.info['ALGORITHMS'] == ('wham',) and record.info['SVTYPE'] == 'DEL') or \
          (record.info['ALGORITHMS'] == ('scramble',) and record.info['HIGH_SR_BACKGROUND'] and record.alts == ('<INS:ME:SVA>',)):
            record.filter.add('HIGH_ALGORITHM_FDR')
        fo.write(record)
CODE
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}

# Add FILTER status for variants that are close to an intron
task AddRetroDelFilters {
  input {
    File vcf
    File intron_reference
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  
  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 2,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    
    python /opt/sv-pipeline/04_variant_resolution/scripts/add_retro_del_filters.py \
      ~{vcf} \
      ~{intron_reference} \
      ~{contig} \
      ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
  }
}


# Final VCF cleanup
task FinalCleanup {
  input {
    File vcf
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: ceil(5.0 + size(vcf, "GB") * 1.5),
    disk_gb: ceil(10.0 + size(vcf, "GB") * 3.0),
    cpu_cores: 2,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail
    
    /opt/sv-pipeline/04_variant_resolution/scripts/rename_after_vcfcluster.py \
      --chrom ~{contig} \
      --prefix ~{prefix} \
      ~{vcf} stdout \
      | bcftools annotate --no-version -e 'SVTYPE=="CNV" && SVLEN<5000' -x INFO/MEMBERS -Oz -o ~{prefix}.vcf.gz
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File final_cleaned_shard = "~{prefix}.vcf.gz"
    File final_cleaned_shard_idx = "~{prefix}.vcf.gz.tbi"
  }
}