version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as Utils

# Workflow to perform depth-based genotyping per batch
# on predicted CPX CNVs

workflow GenotypeCpxCnvsPerBatch {
  input {
    File bin_exclude
    File cpx_bed
    File rd_depth_sep_cutoff
    Int n_per_split_small
    Int n_per_split_large
    Int n_rd_test_bins
    String batch
    File median_file
    File ped_file
    File coverage_file
    File ref_dict

    String linux_docker
    String sv_base_mini_docker
    String sv_pipeline_rdtest_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_ids_from_median
    RuntimeAttr? runtime_override_split_bed_by_size
    RuntimeAttr? runtime_override_rd_genotype

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_concat_melted_genotypes
  }

  File coverage_file_idx = coverage_file + ".tbi"

  call Utils.GetSampleIdsFromMedianCoverageFile {
    input:
      median_file = median_file,
      name = batch,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_override_ids_from_median
  }
  
  call SplitBedBySize {
    input:
      bed=cpx_bed,
      n_per_split_small=n_per_split_small,
      n_per_split_large=n_per_split_large,
      samples_list=GetSampleIdsFromMedianCoverageFile.out_file,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_split_bed_by_size
  }

  scatter (split_bed_file in flatten([SplitBedBySize.small_beds, SplitBedBySize.large_beds])) {
    call RdTestGenotype as RdGenotype {
      input:
        bin_exclude=bin_exclude,
        bed=split_bed_file,
        coverage_file=coverage_file,
        coverage_file_idx=coverage_file_idx,
        median_file=median_file,
        ped_file=ped_file,
        samples_list=GetSampleIdsFromMedianCoverageFile.out_file,
        gt_cutoffs=rd_depth_sep_cutoff,
        n_bins=n_rd_test_bins,
        prefix=basename(split_bed_file, ".bed"),
        ref_dict = ref_dict,
        sv_pipeline_rdtest_docker=sv_pipeline_rdtest_docker,
        runtime_attr_override=runtime_override_rd_genotype
    }
  }

  call MiniTasks.ZcatCompressedFiles as ConcatMeltedGenotypes {
    input:
      shards=RdGenotype.melted_genotypes,
      outfile_name=batch + ".rd_genos.bed.gz",
      filter_command="sort -Vk1,1 -k2,2n -k3,3n",
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_concat_melted_genotypes
  }

  output {
    File genotypes = ConcatMeltedGenotypes.outfile
  }
}


task SplitBedBySize {
  input {
    File bed
    Int n_per_split_small
    Int n_per_split_large
    File samples_list
    Int? size_cutoff_kb
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Int size_cutoff=select_first([size_cutoff_kb, 5])

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  Float input_size = size([bed, samples_list], "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + 2.0 * compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GiB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu -o pipefail

    #First, replace samples in input bed with full list of all samples in batch
    MERGED_BED="newBed_wSamples.bed"
    zcat ~{bed} \
      | awk -F "\t" -v OFS="\t" -v samples=$(paste -s -d, ~{samples_list}) \
        '{ print $1, $2, $3, $4, samples, $6 }' \
      | sort -Vk1,1 -k2,2n -k3,3n \
      > $MERGED_BED
    rm -f ~{bed}

    #Second, split by small vs large CNVs
    SMALL_BED="lt~{size_cutoff}kb.bed"
    LARGE_BED="gt~{size_cutoff}kb.bed"
    touch $SMALL_BED
    touch $LARGE_BED
    awk -F "\t" -v OFS="\t" -v SMALL_BED=$SMALL_BED -v LARGE_BED=$LARGE_BED\
      '{if($3 - $2 < ~{size_cutoff}) {print > SMALL_BED} else {print > LARGE_BED}}' \
      $MERGED_BED
    rm -f $MERGED_BED

    echo -e "$SMALL_BED ~{n_per_split_small}\n$LARGE_BED ~{n_per_split_large}" \
      | while read SIZED_BED LINES_PER_SPLIT; do
          PREFIX="${SIZED_BED%.*}." # remove "bed" from end of file

          N_LINES=$(wc -l < $SIZED_BED)
          N_CHUNKS=$((N_LINES / LINES_PER_SPLIT))
          if [ $N_CHUNKS -gt 1 ]; then
            # split file into multiple chunks
            N_DIGITS=${#N_CHUNKS}

            split -d -a $N_DIGITS -n l/$N_CHUNKS \
              --numeric-suffixes=$(printf "%0${N_DIGITS}d" 1) \
              $SIZED_BED \
              $PREFIX
            rm -f $SIZED_BED
          else
            # only one chunk, just rename file to look like a chunk
            mv $SIZED_BED ${PREFIX}1
          fi
        done
  >>>

  output {
    Array[File] small_beds = glob("lt~{size_cutoff}kb.*")
    Array[File] large_beds = glob("gt~{size_cutoff}kb.*")
  }
}

# Run depth-based genotyping
task RdTestGenotype {
  input {
    File bed
    File coverage_file
    File coverage_file_idx
    File median_file
    File ped_file
    File samples_list
    File bin_exclude
    File gt_cutoffs
    File ref_dict
    Int n_bins
    String prefix
    String sv_pipeline_rdtest_docker
    RuntimeAttr? runtime_attr_override
  }

  parameter_meta {
    coverage_file: {
      localization_optional: true
    }
  }

  # when filtering/sorting/etc, memory usage will likely go up (much of the data will have to
  # be held in memory or disk while working, potentially in a form that takes up more space)
  # NOTE: doubled representation of "bed" below is NOT A TYPO
  #  bincov information is remote-tabixed in, preventing accurate measurement of size, but the bincov
  #  info is at lower resolution than the input bed file, so an upper estimate can be obtained by just
  #  doubling the representation of "bed" in the input size
  Float input_size = size([bed, bed, median_file, ped_file, samples_list, gt_cutoffs], "GiB")
  Float compression_factor = 5.0
  Float base_disk_gb = 5.0
  Float base_mem_gb = 2.0
  RuntimeAttr runtime_default = object {
    mem_gb: base_mem_gb + compression_factor * input_size,
    disk_gb: ceil(base_disk_gb + input_size * (2.0 + 2.0 * compression_factor)),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  Float mem_gb = select_first([runtime_override.mem_gb, runtime_default.mem_gb])
  Int java_mem_mb = ceil(mem_gb * 1000 * 0.8)

  runtime {
    memory: mem_gb + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_rdtest_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -eu

    grep -v "^#" ~{bed} | sort -k1,1V -k2,2n | bedtools merge -i stdin -d 1000000 > merged.bed

    set -o pipefail

    if [ -s merged.bed ]; then
      java -Xmx~{java_mem_mb}M -jar ${GATK_JAR} PrintSVEvidence \
        --sequence-dictionary ~{ref_dict} \
        --evidence-file ~{coverage_file} \
        -L merged.bed \
        -O local.RD.txt.gz
    else
      touch local.RD.txt
      bgzip local.RD.txt
    fi

    tabix -p bed local.RD.txt.gz
    tabix -p bed ~{bin_exclude}

    Rscript /opt/RdTest/RdTest.R \
      -b ~{bed} \
      -c local.RD.txt.gz \
      -m ~{median_file} \
      -f ~{ped_file} \
      -n ~{prefix} \
      -w ~{samples_list} \
      -i ~{n_bins} \
      -r ~{gt_cutoffs} \
      -y ~{bin_exclude} \
      -g TRUE

    if [ -f "~{prefix}.geno" ] && [ -f "~{prefix}.gq" ] ; then
      /opt/sv-pipeline/04_variant_resolution/scripts/merge_RdTest_genotypes.py \
        ~{prefix}.geno ~{prefix}.gq rd.geno.cnv.bed
      sort -k1,1V -k2,2n rd.geno.cnv.bed | uniq | bgzip -c > rd.geno.cnv.bed.gz
    else
      # In case RdTest does not produce output because the input is empty
      echo "" | bgzip -c > rd.geno.cnv.bed.gz
    fi
  >>>

  output {
    File melted_genotypes = "rd.geno.cnv.bed.gz"
  }
}

