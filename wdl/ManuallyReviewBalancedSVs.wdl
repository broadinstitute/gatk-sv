version 1.0

import "Structs.wdl"
import "ManuallyReviewBalancedSVsPerBatch.wdl" as batch_rev
import "TasksMakeCohortVcf.wdl" as tasks

workflow ManuallyReviewBalancedSVs {
  input {
    String prefix

    Array[File] cohort_vcfs  # cohort vcf or vcfs sharded by contig
    Array[File] batch_manta_tloc_vcfs
    Array[File] batch_pe_files

    Array[String] batches
    Array[File] samples_in_batches

    Int min_size

    File cytobands
    File cytobands_index

    Boolean process_ctx
    Boolean process_cpx
    Boolean process_inv

    File generate_pe_tabix_py_script # for development
    File calculate_pe_stats_script # for development

    String sv_base_mini_docker
    String sv_pipeline_docker

    RuntimeAttr? runtime_attr_select_ctx
    RuntimeAttr? runtime_attr_select_cpx
    RuntimeAttr? runtime_attr_select_inv
    RuntimeAttr? runtime_attr_subset_samples
    RuntimeAttr? runtime_attr_combine_tlocs
    RuntimeAttr? runtime_attr_vcf2bed
    RuntimeAttr? runtime_attr_generate_script
    RuntimeAttr? runtime_attr_collect_pe
    RuntimeAttr? runtime_attr_collect_pe_background
    RuntimeAttr? runtime_attr_concat_ctx
    RuntimeAttr? runtime_attr_concat_ctx_background
    RuntimeAttr? runtime_attr_concat_cpx
    RuntimeAttr? runtime_attr_concat_inv
    RuntimeAttr? runtime_attr_calculate_ctx_background
    RuntimeAttr? runtime_attr_calculate_ctx_stats
    RuntimeAttr? runtime_attr_calculate_cpx_stats
    RuntimeAttr? runtime_attr_calculate_inv_stats

  }

  # select svs cohort-wide and merge across contigs
  # then select samples in batch in separate step

  if (process_ctx) {
    scatter (i in range(length(cohort_vcfs))) {
      call SelectSVType as SelectCTX {
        input:
          vcf = cohort_vcfs[i],
          svtype = "CTX",
          prefix="~{prefix}.shard_~{i}",
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_attr_select_ctx
      }
    }

    call tasks.ConcatVcfs as ConcatCTX {
      input:
        vcfs = SelectCTX.svtype_vcf,
        vcfs_idx = SelectCTX.svtype_vcf_index,
        naive = true,
        outfile_prefix = "~{prefix}.CTX",
        sv_base_mini_docker = sv_base_mini_docker
    }

    call LabelArms {
      input:
        ctx_vcf = ConcatCTX.concat_vcf,
        ctx_vcf_index = ConcatCTX.concat_vcf_idx,
        prefix = prefix,
        cytobands = cytobands,
        cytobands_index = cytobands_index,
        sv_pipeline_docker = sv_pipeline_docker
    }

    scatter (i in range(length(batches))) {
      call batch_rev.ManuallyReviewBalancedSVsPerBatch as ManuallyReviewCTXPerBatch {
        input:
          batch = batches[i],
          svtype = "CTX",
          cohort_vcf = LabelArms.witharms_vcf,
          cohort_vcf_index = LabelArms.witharms_vcf_index,
          batch_pe_file = batch_pe_files[i],
          batch_manta_tloc_vcf = batch_manta_tloc_vcfs[i],
          collect_background_pe = true,
          batch_samples = samples_in_batches[i],
          generate_pe_tabix_py_script=generate_pe_tabix_py_script,
          sv_pipeline_docker = sv_pipeline_docker,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_subset_samples=runtime_attr_subset_samples,
          runtime_attr_combine_tlocs=runtime_attr_combine_tlocs,
          runtime_attr_vcf2bed=runtime_attr_vcf2bed,
          runtime_attr_generate_script=runtime_attr_generate_script,
          runtime_attr_collect_pe=runtime_attr_collect_pe,
          runtime_attr_collect_pe_background=runtime_attr_collect_pe_background
      }
    }

    # concatenate per-batch evidence files
    call ConcatEvidences as ConcatCTXEvidences {
      input:
        prefix = "~{prefix}.CTX",
        evidences = ManuallyReviewCTXPerBatch.batch_pe_evidence,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_concat_ctx
    }

    # no overall file header so use different task
    call tasks.ZcatCompressedFiles as ConcatCTXBackground {
      input:
        outfile_name = "~{prefix}.CTX.background_PE.tsv.gz",
        shards = select_all(ManuallyReviewCTXPerBatch.batch_pe_background),
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_concat_ctx_background
    }

    # compute stats about PE evidence
    call CalculatePEBackground as CalculateCTXBackground {
      input:
        prefix = "~{prefix}.CTX",
        background_pe = ConcatCTXBackground.outfile,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_calculate_ctx_background
    }

    call CalculatePEStats as CalculateCTXStats {
      input:
        prefix = "~{prefix}.CTX",
        evidence = ConcatCTXEvidences.concat_evidence,
        background = CalculateCTXBackground.stats,
        calculate_pe_stats_script = calculate_pe_stats_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_calculate_ctx_stats
    }

  }

  if (process_cpx) {
    scatter (i in range(length(cohort_vcfs))) {
      call SelectSVType as SelectCPX {
        input:
          vcf = cohort_vcfs[i],
          svtype = "CPX",
          min_size=min_size,
          prefix="~{prefix}.shard_~{i}",
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_attr_select_cpx
      }
    }

    call tasks.ConcatVcfs as ConcatCPX {
      input:
        vcfs = SelectCPX.svtype_vcf,
        vcfs_idx = SelectCPX.svtype_vcf_index,
        naive = false,
        outfile_prefix = "~{prefix}.CPX",
        sv_base_mini_docker = sv_base_mini_docker
    }

    scatter (i in range(length(batches))) {
      call batch_rev.ManuallyReviewBalancedSVsPerBatch as ManuallyReviewCPXPerBatch {
        input:
          batch = batches[i],
          svtype = "CPX",
          cohort_vcf = ConcatCPX.concat_vcf,
          cohort_vcf_index = ConcatCPX.concat_vcf_idx,
          batch_pe_file = batch_pe_files[i],
          batch_samples = samples_in_batches[i],
          generate_pe_tabix_py_script=generate_pe_tabix_py_script,
          sv_pipeline_docker = sv_pipeline_docker,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_subset_samples=runtime_attr_subset_samples,
          runtime_attr_combine_tlocs=runtime_attr_combine_tlocs,
          runtime_attr_vcf2bed=runtime_attr_vcf2bed,
          runtime_attr_generate_script=runtime_attr_generate_script,
          runtime_attr_collect_pe=runtime_attr_collect_pe
      }
    }

    # concatenate per-batch evidence files
    call ConcatEvidences as ConcatCPXEvidences {
      input:
        prefix = "~{prefix}.CPX",
        evidences = ManuallyReviewCPXPerBatch.batch_pe_evidence,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_concat_cpx
    }

    # compute stats about PE evidence
    call CalculatePEStats as CalculateCPXStats {
      input:
        prefix = "~{prefix}.CPX",
        evidence = ConcatCPXEvidences.concat_evidence,
        calculate_pe_stats_script = calculate_pe_stats_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_calculate_cpx_stats
    }
  }

  if (process_inv) {
    scatter (i in range(length(cohort_vcfs))) {
      call SelectSVType as SelectINV {
        input:
          vcf = cohort_vcfs[i],
          svtype = "INV",
          min_size=min_size,
          prefix="~{prefix}.shard_~{i}",
          sv_pipeline_docker=sv_pipeline_docker,
          runtime_attr_override=runtime_attr_select_inv
      }
    }

    call tasks.ConcatVcfs as ConcatINV {
      input:
        vcfs = SelectINV.svtype_vcf,
        vcfs_idx = SelectINV.svtype_vcf_index,
        naive = true,
        outfile_prefix = "~{prefix}.INV",
        sv_base_mini_docker = sv_base_mini_docker
    }

    scatter (i in range(length(batches))) {
      call batch_rev.ManuallyReviewBalancedSVsPerBatch as ManuallyReviewINVPerBatch {
        input:
          batch = batches[i],
          svtype = "INV",
          cohort_vcf = ConcatINV.concat_vcf,
          cohort_vcf_index = ConcatINV.concat_vcf_idx,
          batch_pe_file = batch_pe_files[i],
          batch_samples = samples_in_batches[i],
          generate_pe_tabix_py_script=generate_pe_tabix_py_script,
          sv_pipeline_docker = sv_pipeline_docker,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_subset_samples=runtime_attr_subset_samples,
          runtime_attr_combine_tlocs=runtime_attr_combine_tlocs,
          runtime_attr_vcf2bed=runtime_attr_vcf2bed,
          runtime_attr_generate_script=runtime_attr_generate_script,
          runtime_attr_collect_pe=runtime_attr_collect_pe
      }
    }

    # concatenate per-batch evidence files
    call ConcatEvidences as ConcatINVEvidences {
      input:
        prefix = "~{prefix}.INV",
        evidences = ManuallyReviewINVPerBatch.batch_pe_evidence,
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_attr_concat_inv
    }

    # compute stats about PE evidence
    call CalculatePEStats as CalculateINVStats {
      input:
        prefix = "~{prefix}.INV",
        evidence = ConcatINVEvidences.concat_evidence,
        calculate_pe_stats_script = calculate_pe_stats_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_calculate_inv_stats
    }
  }

  output {
    File? ctx_evidence = ConcatCTXEvidences.concat_evidence
    File? cpx_evidence = ConcatCPXEvidences.concat_evidence
    File? inv_evidence = ConcatINVEvidences.concat_evidence

    File? ctx_stats = CalculateCTXStats.stats
    File? cpx_stats = CalculateCPXStats.stats
    File? inv_stats = CalculateINVStats.stats
  }
}


task SelectSVType {
  input {
    String prefix
    File vcf
    String svtype
    Int? min_size
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String size_filter_cmd = if defined(min_size) then " && INFO/SVLEN>=~{min_size}" else ""

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 1.0,
                               disk_gb: ceil(10 + 2 * size(vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    bcftools view \
      -i "INFO/SVTYPE=='~{svtype}'~{size_filter_cmd}" \
      ~{vcf} \
      | bcftools +fill-tags \
      -O z \
      -o "~{prefix}.~{svtype}.vcf.gz" \
      -- -t AC,AN,AF

    tabix ~{prefix}.~{svtype}.vcf.gz
  >>>

  output {
    File svtype_vcf = "~{prefix}.~{svtype}.vcf.gz"
    File svtype_vcf_index = "~{prefix}.~{svtype}.vcf.gz.tbi"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: "4 GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    docker: sv_pipeline_docker
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
  }
}


task LabelArms {
  input {
    String prefix
    File ctx_vcf
    File ctx_vcf_index
    File cytobands
    File cytobands_index
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + 2 * size(ctx_vcf, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python <<CODE
import pysam
import gzip

def get_arms(record, cytobands):
  regionA = '{0}:{1}-{1}'.format(record.chrom, record.pos)
  regionB = '{0}:{1}-{1}'.format(record.info['CHR2'], record.info['END2'])

  def _get_arm(region):
    print(region)
    arm = next(cytobands.fetch(region))
    return arm.split()[3][0]

  return _get_arm(regionA), _get_arm(regionB)

cytobands = pysam.TabixFile("~{cytobands}")
vcf = pysam.VariantFile("~{ctx_vcf}")
vcf.header.add_line('##INFO=<ID=END2,Number=1,Type=Integer,Description="Second position">')
outvcf = pysam.VariantFile("~{prefix}.CTX.witharms.vcf.gz", 'w', header=vcf.header)
for record in vcf:
  armA, armB = get_arms(record, cytobands)
  # get CTX arm information
  if armA == armB:
    record.info['CPX_TYPE'] = "CTX_PP/QQ"
  else:
    record.info['CPX_TYPE'] = "CTX_PQ/QP"
  outvcf.write(record)

CODE

    tabix ~{prefix}.CTX.witharms.vcf.gz
  >>>

  output {
    File witharms_vcf = "~{prefix}.CTX.witharms.vcf.gz"
    File witharms_vcf_index = "~{prefix}.CTX.witharms.vcf.gz.tbi"
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


task ConcatEvidences{
  input {
    Array[File] evidences
    String prefix
    String sv_base_mini_docker
   RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5,
    disk_gb: ceil(10 + 2 * size(evidences, "GB")),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    # get only one header
    zcat ~{evidences[0]} | awk 'NR==1' | bgzip -c > ~{prefix}.PE_review.txt.gz
    while read FILE; do
      zcat $FILE | sed '1d'
    done < ~{write_lines(evidences)} \
      | bgzip -c \
      >> ~{prefix}.PE_review.txt.gz

  >>>

  output {
    File concat_evidence = "~{prefix}.PE_review.txt.gz"
  }

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


task CalculatePEStats {
  input {
    String prefix
    File evidence
    File? background
    File calculate_pe_stats_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + 2 * size(evidence, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{calculate_pe_stats_script} \
      -p ~{evidence} \
      ~{"-b " + background} \
      -o ~{prefix}.stats.tsv

    bgzip ~{prefix}.stats.tsv
  >>>

  output {
    File stats = "~{prefix}.stats.tsv.gz"
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

task CalculatePEBackground {
  input {
    String prefix
    File background_pe
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + 2 * size(background_pe, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python << CODE
import gzip
import argparse


def write_background(out_file, background):
    with open(out_file, 'w') as out:
        out.write("\t".join("#SVID background_min1 background_min4 background_min10".split()) + "\n")
        for svid in background:
            counts = background[svid].values()
            min1 = len(counts)
            min4 = len([x for x in counts if x > 3])
            min10 = len([x for x in counts if x > 9])
            out.write(f"{svid}\t{min1}\t{min4}\t{min10}\n")


def increment_count(background, line, curr_svid, curr_samples, pe_header):
    fields = line.strip().split("\t")
    # ignore carriers when counting background PE
    sample_id = fields[pe_header["sample"]]
    if sample_id not in curr_samples:
        if sample_id in background[curr_svid]:
            background[curr_svid][sample_id] += 1
        else:
            background[curr_svid][sample_id] = 1


def process(pe_background, out_file):
    with gzip.open(pe_background, 'rt') as pe:
        pe_header = {x:i for i,x in enumerate("chrom1 pos1 dir1 chrom2 pos2 dir2 sample".split())}
        curr_svid = None
        curr_samples = None
        curr_lines = None
        background = dict()  # {svid: {non_carrier_sample: pe_count}}
        for line in pe:
            if line.startswith("#"):
                curr_svid, samp_string = line.strip().lstrip("#").split("\t")
                curr_samples = set(samp_string.split(","))
                curr_lines = set()
                if curr_svid not in background:
                    background[curr_svid] = dict()
            else:
                if line not in curr_lines:
                    curr_lines.add(line)
                    increment_count(background, line, curr_svid, curr_samples, pe_header)
    write_background(out_file, background)

process("~{background_pe}", "~{prefix}.background.stats.tsv")
CODE

    bgzip ~{prefix}.background.stats.tsv
  >>>

  output {
    File stats = "~{prefix}.background.stats.tsv.gz"
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
