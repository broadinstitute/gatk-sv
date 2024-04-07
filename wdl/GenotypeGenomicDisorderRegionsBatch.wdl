version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as tasks_cohort
import "Utils.wdl" as util

workflow GenotypeGenomicDisorderRegionsBatch {
  input {
    String output_prefix
    String batch_name
    File rd_file
    File median_file
    File depth_sepcutoff_file

    Array[File] cohort_vcfs

    File ploidy_table
    File preprocessed_genomic_disorder_regions_bed
    File genomic_disorder_regions_bed
    File par_bed

    Float min_gdr_overlap_frac
    Boolean plot_subdivisions

    String? revise_args
    File? revise_script
    File? report_script

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_attr_subset_by_samples
    RuntimeAttr? runtime_override_concat_batch
    RuntimeAttr? runtime_gdr_overlapping_variants
    RuntimeAttr? runtime_rdtest_full
    RuntimeAttr? runtime_rdtest_subdiv
    RuntimeAttr? runtime_revise_vcf_batch
    RuntimeAttr? runtime_vcf2bed_before_revise
    RuntimeAttr? runtime_vcf2bed_after_revise
    RuntimeAttr? runtime_rdtest_before_revise
    RuntimeAttr? runtime_rdtest_after_revise
    RuntimeAttr? runtime_create_report
  }

  call util.GetSampleIdsFromMedianCoverageFile {
    input:
      median_file = median_file,
      name = batch_name,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_ids_from_median
  }

  scatter (i in range(length(cohort_vcfs))) {
    call util.SubsetVcfBySamplesList {
      input:
        vcf = cohort_vcfs[i],
        list_of_samples = GetSampleIdsFromMedianCoverageFile.out_file,
        outfile_name = "~{output_prefix}.shard_~{i}",
        remove_samples = false,
        remove_private_sites = true,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_subset_by_samples
    }
  }
  call tasks_cohort.ConcatVcfs {
    input:
      vcfs = SubsetVcfBySamplesList.vcf_subset,
      vcfs_idx = SubsetVcfBySamplesList.vcf_subset_index,
      outfile_prefix = "~{output_prefix}.concat",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_override_concat_batch
  }
  call GetGDROverlappingVariants {
    input:
    vcf = ConcatVcfs.concat_vcf,
    genomic_disorder_regions_bed = genomic_disorder_regions_bed,
    prefix = "~{output_prefix}",
    min_gdr_overlap_frac_plotting = min_gdr_overlap_frac,
    sv_pipeline_docker = sv_pipeline_docker,
    runtime_attr_override = runtime_gdr_overlapping_variants
  }
  # Run RdTest and generate plots on variants overlapping one or more GDRs (plotted carriers highlighted)
  call RunRdTest as RunRdTestVariantsOverlappingGDR {
    input:
      output_prefix = "rdtest_var2gdr_~{batch_name}",
      rdtest_bed = GetGDROverlappingVariants.variants_bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_full
  }
  # Run RdTest and generate plots on GDRs overlapping one or more variants (plotted carriers highlighted)
  call RunRdTest as RunRdTestGDROverlappingVariants {
    input:
      output_prefix = "rdtest_gdr2var_~{batch_name}",
      rdtest_bed = GetGDROverlappingVariants.gdr_bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_full
  }

  # Run RdTest and generate plots over the all GD regions (plotted carriers NOT highlighted)
  call RunRdTest as RunRdTestFullRegions {
    input:
      output_prefix = "rdtest_full_~{batch_name}",
      rdtest_bed = genomic_disorder_regions_bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      inject_sample = true,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_full
  }
  # Run RdTest and generate plots over subdivided GD regions (plotted carriers NOT highlighted)
  call RunRdTest as RunRdTestSubdivision {
    input:
      output_prefix = "rdtest_subdiv_~{batch_name}",
      rdtest_bed = preprocessed_genomic_disorder_regions_bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      depth_sepcutoff = depth_sepcutoff_file,
      inject_sample = true,
      do_plot = plot_subdivisions,
      do_genotyping = true,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_subdiv
  }
  call ReviseGenomicDisorderRegions {
    input:
      prefix = "~{output_prefix}.revise_gdr",
      rdtest_tars = [RunRdTestSubdivision.out],
      batch_name = batch_name,
      vcf = ConcatVcfs.concat_vcf,
      new_record_prefix = "~{output_prefix}_revise_gdr_",
      ploidy_table = ploidy_table,
      genomic_disorder_regions_bed = preprocessed_genomic_disorder_regions_bed,
      par_bed = par_bed,
      min_gdr_overlap_frac = min_gdr_overlap_frac,
      args = revise_args,
      script = revise_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_revise_vcf_batch
  }
  call VcfToBed as VcfToBedRevisedBeforeUpdateRecords {
    input:
      vcf = ReviseGenomicDisorderRegions.revised_before_update_vcf,
      prefix = "~{output_prefix}.before_revise",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_vcf2bed_before_revise
  }
  call VcfToBed as VcfToBedRevisedAfterUpdateRecords {
    input:
      vcf = ReviseGenomicDisorderRegions.revised_after_update_vcf,
      prefix = "~{output_prefix}.after_revise",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_vcf2bed_after_revise
  }
  call VcfToBed as VcfToBedNewRecords {
    input:
      vcf = ReviseGenomicDisorderRegions.revised_new_records_vcf,
      prefix = "~{output_prefix}.new_records",
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_vcf2bed_after_revise
  }
  call RunRdTest as RunRdTestBeforeUpdate {
    input:
      output_prefix = "rdtest_before_revise_~{batch_name}",
      rdtest_bed = VcfToBedRevisedBeforeUpdateRecords.bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_before_revise
  }
  call RunRdTest as RunRdTestAfterUpdate {
    input:
      output_prefix = "rdtest_after_revise_~{batch_name}",
      rdtest_bed = VcfToBedRevisedAfterUpdateRecords.bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_after_revise
  }
  call RunRdTest as RunRdTestNewRecords {
    input:
      output_prefix = "rdtest_new_~{batch_name}",
      rdtest_bed = VcfToBedNewRecords.bed,
      rd_file = rd_file,
      rd_index = rd_file + ".tbi",
      median_file = median_file,
      do_plot = true,
      do_genotyping = false,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_rdtest_after_revise
  }

  call CreateRevisionReport {
    input:
      prefix = "~{output_prefix}.gdr_report",
      rdtest_tars = [RunRdTestVariantsOverlappingGDR.out, RunRdTestGDROverlappingVariants.out,
                    RunRdTestFullRegions.out, RunRdTestSubdivision.out, RunRdTestBeforeUpdate.out,
                    RunRdTestAfterUpdate.out, RunRdTestNewRecords.out],
      batch_name = batch_name,
      revised_genotypes_tsv = ReviseGenomicDisorderRegions.revised_genotypes_tsv,
      revision_manifest_tsv = ReviseGenomicDisorderRegions.revision_manifest_tsv,
      genomic_disorder_regions_bed = genomic_disorder_regions_bed,
      script = report_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_create_report
  }
  output{
    File batch_gdr_report = CreateRevisionReport.out

    File batch_gdr_revised_before_update_vcf = ReviseGenomicDisorderRegions.revised_before_update_vcf
    File batch_gdr_revised_before_update_index = ReviseGenomicDisorderRegions.revised_before_update_index

    File batch_gdr_revised_after_update_vcf = ReviseGenomicDisorderRegions.revised_after_update_vcf
    File batch_gdr_revised_after_update_index = ReviseGenomicDisorderRegions.revised_after_update_index

    File batch_gdr_new_records_vcf = ReviseGenomicDisorderRegions.revised_new_records_vcf
    File batch_gdr_new_records_index = ReviseGenomicDisorderRegions.revised_new_records_index

    File batch_gdr_revised_genotypes_tsv = ReviseGenomicDisorderRegions.revised_genotypes_tsv
    File batch_gdr_revision_manifest_tsv = ReviseGenomicDisorderRegions.revision_manifest_tsv

    File batch_subsetted_vcf = ConcatVcfs.concat_vcf
    File batch_subsetted_index = ConcatVcfs.concat_vcf_idx
  }
}

task RunRdTest {
  input{
    String output_prefix
    File rdtest_bed
    File rd_file
    File rd_index
    File median_file
    File? depth_sepcutoff  # Required if do_genotyping = true
    Boolean inject_sample = false  # If missing a samples column, fill one using the first sample from the median_file
    Boolean do_plot
    Boolean do_genotyping
    Int large_size_cutoff = 1000000
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50.0 + size(rd_file, "GiB") + size(rdtest_bed, "GiB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    if ~{inject_sample}; then
      # Inject one sample from the batch into the 5th column
      SAMPLE=$(awk -F'\t' '{ if (NR==1) {print $1} }' ~{median_file})
      awk -F'\t' -v OFS='\t' -v s="$SAMPLE" '{print $1,$2,$3,$4,s,$5}' ~{rdtest_bed} > intervals.bed
    else
      # Remove sites with no carriers, which isn't currently supported by RdTest
      awk -F'\t' -v OFS='\t' '$5!=""' ~{rdtest_bed} > intervals.bed
    fi
    mkdir ~{output_prefix}/
    Rscript /opt/RdTest/RdTest.R \
      ~{if do_genotyping then "-g TRUE -v TRUE" else ""} \
      ~{if do_plot then "-p TRUE" else ""} \
      ~{"-r " + depth_sepcutoff} \
      -b intervals.bed \
      -c ~{rd_file} \
      -m ~{median_file} \
      -n ~{output_prefix} \
      -s ~{large_size_cutoff} \
      -o ~{output_prefix}
    tar czvf ~{output_prefix}.tar.gz ~{output_prefix}/
  >>>
  output{
    File out = "~{output_prefix}.tar.gz"
  }
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

task ReviseGenomicDisorderRegions {
  input{
    String prefix
    Array[File] rdtest_tars
    String batch_name
    File vcf
    File ploidy_table
    File genomic_disorder_regions_bed
    File par_bed
    String new_record_prefix
    Float min_gdr_overlap_frac
    String? args
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50.0 + size(vcf, "GiB") * 3 + size(rdtest_tars, "GiB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    mkdir rdtest/
    while read -r FILE; do
      tar xzf $FILE -C rdtest/
    done < ~{write_lines(rdtest_tars)}
    ls rdtest/*/*.geno > geno_files.list
    python ~{default="/opt/src/sv-pipeline/scripts/revise_genomic_disorder_regions.py" script} \
      ~{args} \
      --min-region-overlap ~{min_gdr_overlap_frac} \
      --vcf ~{vcf} \
      --batch ~{batch_name} \
      --geno-file-list geno_files.list \
      --ploidy-table ~{ploidy_table} \
      --region-bed ~{genomic_disorder_regions_bed} \
      --par-bed ~{par_bed} \
      --out ~{prefix}
  >>>
  output {
    File revised_before_update_vcf = "~{prefix}.revised_before_update.vcf.gz"
    File revised_before_update_index = "~{prefix}.revised_before_update.vcf.gz.tbi"
    File revised_after_update_vcf = "~{prefix}.revised_after_update.vcf.gz"
    File revised_after_update_index = "~{prefix}.revised_after_update.vcf.gz.tbi"
    File revised_new_records_vcf = "~{prefix}.new_records.vcf.gz"
    File revised_new_records_index = "~{prefix}.new_records.vcf.gz.tbi"
    File revised_genotypes_tsv = "~{prefix}.revised_genotypes.tsv.gz"
    File revision_manifest_tsv = "~{prefix}.revision_manifest.tsv.gz"
  }
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

task CreateRevisionReport {
  input{
    String prefix
    String batch_name
    Array[File] rdtest_tars
    File revised_genotypes_tsv
    File revision_manifest_tsv
    File genomic_disorder_regions_bed
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(100.0 + size(rdtest_tars, "GiB") * 5),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
  command <<<
    set -euxo pipefail
    mkdir rdtest/
    while read -r FILE; do
      tar xzf $FILE -C rdtest/
    done < ~{write_lines(rdtest_tars)}
    find rdtest/ -name '*.jpg' > images.list
    echo "~{batch_name}" > batches.list
    python ~{default="/opt/src/sv-pipeline/scripts/create_gd_report.py" script} \
      --image-list images.list \
      --batch-list batches.list \
      --region-bed ~{genomic_disorder_regions_bed} \
      --genotypes ~{revised_genotypes_tsv} \
      --manifest ~{revision_manifest_tsv} \
      --out ~{prefix}
    # Create output tarball
    mkdir ~{prefix}/
    mv rdtest ~{prefix}/
    mv ~{prefix}.~{batch_name}.html ~{prefix}/
    tar czf ~{prefix}.tar.gz ~{prefix}/
  >>>
  output {
    File out = "~{prefix}.tar.gz"
  }
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

task GetGDROverlappingVariants {
  input {
    File vcf
    File genomic_disorder_regions_bed
    Float min_gdr_overlap_frac_plotting
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(50 + size(vcf, "GiB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File variants_bed = "~{prefix}.variants_with_gdr_overlaps.bed"
    File gdr_bed = "~{prefix}.gdr_with_variant_overlaps.bed"
  }
  command <<<
    set -euxo pipefail

    # Min size should be half the smallest interval length
    MIN_SIZE=$(awk -F'\t' '{print $3-$2}' ~{genomic_disorder_regions_bed} | sort -n | awk -v f=~{min_gdr_overlap_frac_plotting} '{if (NR==1) {print int($1 * f)}}')

    # Get non-ref DEL/DUP bed records
    # Note we need to use double-quotes around the filtering expression since we are referencing a shell variable
    bcftools view -i "(SVTYPE==\"DEL\" || SVTYPE==\"DUP\") && SVLEN>=$MIN_SIZE && COUNT(GT=\"alt\")>0" ~{vcf} \
      | svtk vcf2bed --no-header - intervals.bed

    # Separate DEL and DUP records
    awk -F'\t' -v OFS='\t' '$5=="DEL"' intervals.bed > intervals.DEL.bed
    awk -F'\t' -v OFS='\t' '$5=="DUP"' intervals.bed > intervals.DUP.bed
    awk -F'\t' -v OFS='\t' '$5=="DEL"' ~{genomic_disorder_regions_bed} > gdr.DEL.bed
    awk -F'\t' -v OFS='\t' '$5=="DUP"' ~{genomic_disorder_regions_bed} > gdr.DUP.bed

    # Get variants overlapping at least X% of a GDR
    # Records are named <VARIANT_ID>__<GDR_ID>
    # Note we swap columns 5 and 6 (SVTYPE and SAMPLES) for RdTest
    bedtools intersect -wo -F ~{min_gdr_overlap_frac_plotting} -a intervals.DEL.bed -b gdr.DEL.bed \
      | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4"__"$10,$6,$5}' \
      > intervals.DEL.gdr_overlaps.bed
    bedtools intersect -wo -F ~{min_gdr_overlap_frac_plotting} -a intervals.DUP.bed -b gdr.DUP.bed \
    | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$4"__"$10,$6,$5}' \
    > intervals.DUP.gdr_overlaps.bed
    cat intervals.DEL.gdr_overlaps.bed intervals.DUP.gdr_overlaps.bed \
      | sort -k1,1V -k2,2n -k3,3n \
      > ~{prefix}.variants_with_gdr_overlaps.bed

    # Get GDRs overlapped at least X% by a variant
    # Records are named <GDR_ID>__<VARIANT_ID>
    # Note we swap columns 5 and 6 (SVTYPE and SAMPLES) for RdTest
    bedtools intersect -wo -F ~{min_gdr_overlap_frac_plotting} -a intervals.DEL.bed -b gdr.DEL.bed \
      | awk -F'\t' -v OFS='\t' '{print $7,$8,$9,$10"__"$4,$6,$5}' \
      > gdr.DEL.variant_overlaps.bed
    bedtools intersect -wo -F ~{min_gdr_overlap_frac_plotting} -a intervals.DUP.bed -b gdr.DUP.bed \
      | awk -F'\t' -v OFS='\t' '{print $7,$8,$9,$10"__"$4,$6,$5}' \
      > gdr.DUP.variant_overlaps.bed
    cat gdr.DEL.variant_overlaps.bed gdr.DUP.variant_overlaps.bed \
      | sort -k1,1V -k2,2n -k3,3n \
      > ~{prefix}.gdr_with_variant_overlaps.bed

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
                               disk_gb: ceil(50 + size(vcf, "GiB") * 2),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  output {
    File bed = "~{prefix}.bed"
  }
  command <<<
    set -euxo pipefail
    svtk vcf2bed --no-header ~{vcf} intervals.bed
    # Swap last two columns for RdTest
    awk -F"\t" -v OFS="\t" '{print $1, $2, $3, $4, $6, $5}' intervals.bed > ~{prefix}.bed
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