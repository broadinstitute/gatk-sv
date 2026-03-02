version 1.0

import "Structs.wdl"
import "ManuallyReviewBalancedSVsPerBatch.wdl" as batch_rev
import "TasksMakeCohortVcf.wdl" as tasks

workflow ManuallyReviewBalancedSVs {
  input {
    String prefix

    Array[File] cohort_vcfs  # cohort vcf or vcfs sharded by contig
    Array[File] batch_pe_files

    Array[String] batches
    Array[File] samples_in_batches
    File batch_membership

    String vid_include_cmd

    File generate_pe_tabix_py_script # for development
    File calculate_pe_stats_script # for development
    File plot_pe_script

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
    RuntimeAttr? runtime_attr_concat_cpx
    RuntimeAttr? runtime_attr_calculate_cpx_stats
    RuntimeAttr? runtime_attr_plot_pe
    RuntimeAttr? runtime_attr_get_batch_idxs

  }

  # select svs cohort-wide and merge across contigs
  # then select samples in batch in separate step


  scatter (i in range(length(cohort_vcfs))) {
    call SelectSVType as SelectCPX {
      input:
        vcf = cohort_vcfs[i],
        vid_include_cmd=vid_include_cmd,
        batch_membership=batch_membership,
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

  call GetBatchIndexes {
    input:
      batches = batches,
      carrier_batches = SelectCPX.carrier_batches,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override = runtime_attr_get_batch_idxs
  }

  scatter (idx in range(length(GetBatchIndexes.batch_idxs))) {
    Int i = GetBatchIndexes.batch_idxs[idx]
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

  call PlotPEEvidence {
    input:
      prefix = "~{prefix}.CPX",
      evidence = ConcatCPXEvidences.concat_evidence,
      plot_pe_script = plot_pe_script,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_plot_pe
  }


  output {
    File cpx_evidence = ConcatCPXEvidences.concat_evidence
    File cpx_stats = CalculateCPXStats.stats
    File pe_plots = PlotPEEvidence.plots
  }
}


task SelectSVType {
  input {
    String prefix
    File vcf
    File batch_membership
    String vid_include_cmd
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

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
      -i '~{vid_include_cmd}' \
      ~{vcf} \
      -O z \
      -o "~{prefix}.selected.vcf.gz"

    tabix ~{prefix}.selected.vcf.gz

    python <<CODE
import pysam
carriers = set()
with pysam.VariantFile("~{prefix}.selected.vcf.gz") as vcf:
    for record in vcf:
        counter = 0
        for s,gt in record.samples.items():
            if any([a is not None and a > 0 for a in gt['GT']]):
                carriers.add(s)
                counter += 1
                if counter >= 3:
                    break

samp_to_batch = dict()
with open("~{batch_membership}", 'r') as inp:
    for line in inp:
        samp, batch = line.strip("\n").split("\t")
        samp_to_batch[samp] = batch

carrier_batches = {samp_to_batch[x] for x in carriers}

with open("~{prefix}.batches.txt", 'w') as out:
    for i in carrier_batches:
        out.write(f"{i}\n")

CODE


  >>>

  output {
    File svtype_vcf = "~{prefix}.selected.vcf.gz"
    File svtype_vcf_index = "~{prefix}.selected.vcf.gz.tbi"
    File carrier_batches = "~{prefix}.batches.txt"
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


task GetBatchIndexes {
  input {
    Array[String] batches
    Array[File] carrier_batches
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

  command <<<
    set -euo pipefail

    python <<CODE
carrier_batches = set()
for f in ["~{sep='", "' carrier_batches}"]:
    with open(f, 'r') as inp:
        for line in inp:
            carrier_batches.add(line.strip("\n"))

batch_idxs = [i for i, x in enumerate(["~{sep='", "' batches}"]) if x in carrier_batches]
with open("batch_idxs.txt", 'w') as out:
    for i in batch_idxs:
        out.write(f"{i}\n")
CODE


  >>>

  output {
    Array[Int] batch_idxs = read_lines("batch_idxs.txt")
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


task PlotPEEvidence {
  input {
    String prefix
    File evidence
    File plot_pe_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
                               cpu_cores: 1,
                               mem_gb: 3.75,
                               disk_gb: ceil(10 + 10 * size(evidence, "GB")),
                               boot_disk_gb: 10,
                               preemptible_tries: 3,
                               max_retries: 1
                             }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python <<CODE
import gzip
with gzip.open(~{evidence}, 'rt') as pe, open("filenames.txt", 'w') as fn:
    first = True
    curr_lines = None
    curr_fname = None
    for line in pe:
        fields = line.strip().lstrip("#").split("\t")
        if first:
            first = False
        elif line.startswith("#"):
            if curr_lines is not None and len(curr_lines) > 0:
                with open(curr_fname, 'w') as out:
                    for l in curr_lines:
                        out.write(l)
                fn.write(curr_fname + '\n')
            vid = fields[3]
            samp = fields[13]
            curr_fname = f"{vid}.{samp}.pe"
            curr_lines = []
        else:
            if line not in curr_lines:
                curr_lines.append(line)
    # handle last variant
    if len(curr_lines) > 0:
        with open(curr_fname, 'w') as out:
            for l in curr_lines:
                out.write(l)
        fn.write(curr_fname + '\n')
CODE

    mkdir plots
    while read f; do
      Rscript ~{plot_pe_script} $f "plots/$f.pdf"
    done < filenames.txt

    tar -czvf ~{prefix}.pe_plots.tar.gz plots/
  >>>

  output {
    File plots = "~{prefix}.plots.tar.gz"
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
