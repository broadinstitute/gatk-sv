version 1.0

import "Structs.wdl"
import "TasksMakeCohortVcf.wdl" as MiniTasks
import "Utils.wdl" as util
import "AnnotateFunctionalConsequences.wdl" as func
import "AnnotateExternalAFPerShard.wdl" as eaf

# Perform annotation per contig

workflow ShardedAnnotateVcf {

  input {
    File vcf
    File vcf_idx
    String prefix
    String contig

    File ploidy_table
    File? apply_filters_script
    Float? no_call_rate_cutoff

    File? protein_coding_gtf
    File? noncoding_bed
    Int? promoter_window
    Int? max_breakend_as_cnv_length
    String? svannotate_additional_args

    File? sample_pop_assignments  # Two-column file with sample ID & pop assignment. "." for pop will ignore sample
    File? sample_keep_list
    File? sample_id_rename_map
    File? ped_file                # Used for M/F AF calculations
    File par_bed
    File? allosomes_list
    File? compute_afs_script
    Int   sv_per_shard

    File? ref_bed              # File with external allele frequencies
    String? ref_prefix         # prefix name for external AF call set (required if ref_bed set)
    Array[String]? population  # populations to annotate external AF for (required if ref_bed set)

    String sv_pipeline_docker
    String sv_base_mini_docker
    String gatk_docker

    RuntimeAttr? runtime_attr_svannotate
    RuntimeAttr? runtime_attr_compute_AFs
    RuntimeAttr? runtime_attr_subset_vcf_by_samples_list
    RuntimeAttr? runtime_attr_modify_vcf
    RuntimeAttr? runtime_attr_split_ref_bed
    RuntimeAttr? runtime_attr_split_query_vcf
    RuntimeAttr? runtime_attr_bedtools_closest
    RuntimeAttr? runtime_attr_select_matched_svs
    RuntimeAttr? runtime_attr_scatter_vcf
    RuntimeAttr? runtime_attr_estimate_rd_cn_af
    RuntimeAttr? runtime_attr_rename_samples
    RuntimeAttr? runtime_attr_apply_filters
  }

  if (defined(ref_bed)) {
    call eaf.SplitRefBed {
      input:
        bed = select_first([ref_bed]),
        contig = contig,
        sv_base_mini_docker = sv_base_mini_docker,
        runtime_attr_override = runtime_attr_split_ref_bed
    }
  }

  call MiniTasks.ScatterVcf {
    input:
      vcf = vcf,
      vcf_index = vcf_idx,
      prefix = prefix,
      records_per_shard = sv_per_shard,
      contig = contig,
      sv_pipeline_docker = sv_pipeline_docker,
      runtime_attr_override = runtime_attr_scatter_vcf
  }

  scatter (i in range(length(ScatterVcf.shards))) {
    String shard_prefix = "~{prefix}.~{contig}.~{i}"

    if (defined(sample_keep_list)) {
      call util.SubsetVcfBySamplesList {
        input:
          vcf = ScatterVcf.shards[i],
          list_of_samples = select_first([sample_keep_list]),
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_subset_vcf_by_samples_list
      }
    }

    if (defined(sample_id_rename_map)) {
      call RenameVcfSamplesTask {
        input:
          vcf = select_first([SubsetVcfBySamplesList.vcf_subset, ScatterVcf.shards[i]]),
          vcf_idx = SubsetVcfBySamplesList.vcf_subset_index,
          sample_id_rename_map = select_first([sample_id_rename_map]),
          prefix = shard_prefix,
          check_rename_all_samples = false,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_override = runtime_attr_rename_samples
      }
    }

    if (defined (protein_coding_gtf) || defined (noncoding_bed)) {
      call func.AnnotateFunctionalConsequences {
        input:
          vcf = select_first([RenameVcfSamplesTask.out, SubsetVcfBySamplesList.vcf_subset, ScatterVcf.shards[i]]),
          vcf_index = select_first([RenameVcfSamplesTask.out_index, SubsetVcfBySamplesList.vcf_subset_index]),
          prefix = shard_prefix,
          protein_coding_gtf = protein_coding_gtf,
          noncoding_bed = noncoding_bed,
          promoter_window = promoter_window,
          max_breakend_as_cnv_length = max_breakend_as_cnv_length,
          additional_args = svannotate_additional_args,
          gatk_docker = gatk_docker,
          runtime_attr_svannotate = runtime_attr_svannotate
      }
    }

    call ApplyFilters {
      input:
        vcf = select_first([AnnotateFunctionalConsequences.annotated_vcf, RenameVcfSamplesTask.out, SubsetVcfBySamplesList.vcf_subset, ScatterVcf.shards[i]]),
        prefix = shard_prefix,
        no_call_rate_cutoff = no_call_rate_cutoff,
        ploidy_table = ploidy_table,
        filter_reference_artifacts = true,
        remove_zero_carrier_sites = true,
        apply_filters_script = apply_filters_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_apply_filters
    }

    # Compute AC, AN, and AF per population & sex combination
    call ComputeAFs {
      input:
        vcf = ApplyFilters.filtered_vcf,
        prefix = shard_prefix,
        sample_pop_assignments = sample_pop_assignments,
        ped_file = ped_file,
        par_bed = par_bed,
        allosomes_list = allosomes_list,
        script = compute_afs_script,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_compute_AFs
    }

    call EstimateRdCnFrequency {
      input:
        vcf = ComputeAFs.af_vcf,
        vcf_idx = ComputeAFs.af_vcf_idx,
        ploidy_table = ploidy_table,
        par_bed = par_bed,
        sv_pipeline_docker = sv_pipeline_docker,
        runtime_attr_override = runtime_attr_estimate_rd_cn_af
    }

    if (defined(ref_bed)) {
      call eaf.AnnotateExternalAFPerShard {
        input:
          vcf = ComputeAFs.af_vcf,
          vcf_idx = ComputeAFs.af_vcf_idx,
          split_ref_bed_del = select_first([SplitRefBed.del]),
          split_ref_bed_dup = select_first([SplitRefBed.dup]),
          split_ref_bed_ins = select_first([SplitRefBed.ins]),
          split_ref_bed_inv = select_first([SplitRefBed.inv]),
          split_ref_bed_bnd = select_first([SplitRefBed.bnd]),
          population = select_first([population]),
          ref_prefix = select_first([ref_prefix]),
          prefix = "~{prefix}.~{contig}.~{i}",
          sv_base_mini_docker = sv_base_mini_docker,
          sv_pipeline_docker = sv_pipeline_docker,
          runtime_attr_modify_vcf = runtime_attr_modify_vcf,
          runtime_attr_split_query_vcf = runtime_attr_split_query_vcf,
          runtime_attr_bedtools_closest = runtime_attr_bedtools_closest,
          runtime_attr_select_matched_svs = runtime_attr_select_matched_svs
      }
    }
  }

  output {
    Array[File] sharded_annotated_vcf = if (defined (ref_bed)) then select_all(AnnotateExternalAFPerShard.annotated_vcf) else EstimateRdCnFrequency.rd_cn_freq_vcf
    Array[File] sharded_annotated_vcf_idx = if (defined (ref_bed)) then select_all(AnnotateExternalAFPerShard.annotated_vcf_tbi) else EstimateRdCnFrequency.rd_cn_freq_vcf
  }
}


task ApplyFilters {
  input {
    File vcf
    String prefix
    File ploidy_table
    String? cohort_id
    Int? shard_index
    Float? no_call_rate_cutoff
    Boolean filter_reference_artifacts
    Boolean remove_zero_carrier_sites
    File? apply_filters_script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 5,
    disk_gb: ceil(10.0 + 2 * size(vcf, "GiB")),
    boot_disk_gb: 30,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python ~{select_first([apply_filters_script, "/opt/sv-pipeline/scripts/apply_ncr_and_ref_artifact_filters.py"])} \
      --vcf ~{vcf} \
      --out ~{prefix}.vcf.gz \
      --ploidy-table ~{ploidy_table} \
      --ncr-threshold ~{no_call_rate_cutoff} \
      ~{"--cohort-id " + cohort_id} \
      ~{"--shard-index " + shard_index} \
      ~{if (filter_reference_artifacts) then "--filter-reference-artifacts" else ""} \
      ~{if (remove_zero_carrier_sites) then "--remove-zero-carrier-sites" else ""}

    touch ~{cohort_id}.vid_map.tsv

  >>>

  output {
    File filtered_vcf = "~{prefix}.vcf.gz"
    File id_rename_map = "~{cohort_id}.vid_map.tsv"
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


task RenameVcfSamplesTask {
  input {
    File vcf
    File? vcf_idx
    File sample_id_rename_map
    String prefix
    Boolean? check_rename_all_samples = true
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 1.0,
    disk_gb: ceil(10 + size(vcf, "GiB") * 2.0),
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail
    if ~{check_rename_all_samples}; then
      bcftools query -l ~{vcf} | sort > samples.list
      python <<CODE
with open("~{sample_id_rename_map}", 'r') as rename, open("samples.list", 'r') as header:
  all_to_rename = set()
  for line in rename:
    all_to_rename.add(line.strip("\n").split("\t")[0])  # create set of sample IDs to rename
  for line in header:
    sample = line.strip("\n")
    if sample not in all_to_rename:
      raise ValueError(f"Sample {sample} is in the VCF header but not in the renaming map")
CODE
    fi
    bcftools reheader --samples ~{sample_id_rename_map} -o ~{prefix}.vcf.gz ~{vcf}
    tabix ~{prefix}.vcf.gz
  >>>

  output {
    File out = "~{prefix}.vcf.gz"
    File out_index = "~{prefix}.vcf.gz.tbi"
  }
}


task EstimateRdCnFrequency {
  input {
    File vcf
    File? vcf_idx
    File par_bed
    File ploidy_table
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String output_file = basename(vcf, ".vcf.gz") + ".rd_cn_freq.vcf.gz"

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10.0 + size(vcf, "GB") * 2.0),
    cpu_cores: 1,
    preemptible_tries: 3,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_attr.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_attr.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_attr.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
import pysam
from collections import defaultdict


PLOIDY = dict()
PAR = defaultdict(list)


with open("~{par_bed}", 'r') as inp:
    for line in inp:
        chrom, start, end, label = line.strip("\n").split("\t")
        PAR[chrom].append((int(start), int(end)))  # PAR[chrom] = [(start1, end1), (start2, end2)]


with open("~{ploidy_table}", 'r') as inp:
    first = True
    header = []
    for line in inp:
        fields = line.strip("\n").split("\t")
        if first:
            header = fields
            first = False
            continue
        PLOIDY[fields[0]] = dict()  # first column = sample
        for i in range(1, len(fields)):
            PLOIDY[fields[0]][header[i]] = int(fields[i])  # PLOIDY[sample][contig] = ploidy


with pysam.VariantFile("~{vcf}") as vcf_in:
    out_header = vcf_in.header
    out_header.add_line('##INFO=<ID=RD_CN_ESTIMATED_AF,Number=1,Type=Float,Description="Estimated AF from RD_CN for CNVs">')
    with pysam.VariantFile("~{output_file}", 'w', header=out_header) as vcf_out:
        for record in vcf_in:
            svtype = record.info.get("SVTYPE", None)
            if svtype in ("DEL", "DUP", "CNV"):
                ac = 0
                an = 0
                par = False
                if record.chrom == "chrX" or record.chrom == "chrY":
                    for start, end in PAR[record.chrom]:
                        if float(min(record.stop, end) - max(record.pos, start))/(record.stop - record.pos) > 0.5:
                            par = True  # if variant is >50% in PAR
                for sample in record.samples:
                    rd_cn = record.samples[sample].get("RD_CN", None)
                    ecn = PLOIDY[sample][record.chrom]
                    if par and ecn == 1:
                        ecn = 2  # males expected CN in PAR on chrX = 2; no variants contained within par on chrY in gnomAD
                    if rd_cn is None:
                        continue
                    else:
                        an += ecn
                    if rd_cn == ecn - 1 or rd_cn == ecn + 1:
                        ac += 1
                    elif rd_cn < ecn - 1 or rd_cn > ecn + 1:
                        ac += ecn

                if an > 0:
                    record.info["RD_CN_ESTIMATED_AF"] = float(ac)/float(an)
            vcf_out.write(record)


CODE

    tabix -p vcf "~{output_file}"
  >>>

  output {
    File rd_cn_freq_vcf = output_file
    File rd_cn_freq_vcf_idx = output_file + ".tbi"
  }
}


task ComputeAFs {
  input {
    File vcf
    String prefix
    File? sample_pop_assignments
    File? ped_file
    File? par_bed
    File? allosomes_list
    File? script
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }
  RuntimeAttr default_attr = object {
    cpu_cores: 1,
    mem_gb: 1.5,
    disk_gb: ceil(20 + size(vcf, "GB") * 2),
    boot_disk_gb: 10,
    preemptible_tries: 3,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail
    ~{default="/opt/sv-pipeline/05_annotation/scripts/compute_AFs.py" script} "~{vcf}" stdout \
      ~{"-p " + sample_pop_assignments} \
      ~{"-f " + ped_file} \
      ~{"--par " + par_bed} \
      ~{"--allosomes-list " + allosomes_list} \
    | bgzip -c \
    > "~{prefix}.wAFs.vcf.gz"

    tabix -p vcf "~{prefix}.wAFs.vcf.gz"
  >>>

  output {
    File af_vcf = "~{prefix}.wAFs.vcf.gz"
    File af_vcf_idx = "~{prefix}.wAFs.vcf.gz.tbi"
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
