version 1.0

# Computes SV-SNV linkage disequilibrium (r-squared) between common SVs and
# nearby common SNVs, stratified by population, and annotates SVs with
# genomic context (RepeatMasker / SegDup / SimpleRepeat overlap).
#
# Adapted from the Talkowski_Pangenie/KAGE 1KG Terra notebook
# compute_ld_dv_intsv_trgt_v2.ipynb, generalized to take separate genome-wide
# SNV and SV VCFs (rather than one combined VCF) and to loop over populations.

import "Structs.wdl"

workflow ComputeSvSnvLd {
  input {
    File snv_vcf
    File sv_vcf
    File sample_pop_map           # TSV: sample \t population

    File repeatmasker_bed         # hg38 RepeatMasker regions, e.g. module18_annotate_genomic_context/references/hg38.RM.sorted.merged.bed.gz
    File segdup_bed               # hg38 segmental duplications
    File simplerepeat_bed         # hg38 simple repeats
    File annotate_genomic_context_script   # module18_annotate_genomic_context/annotate_genomic_context.sh
    File annotate_genomic_context_helper_r # module18_annotate_genomic_context/annotate_genomic_context_helper.R

    String output_prefix
    String staging_gcs_dir        # gs:// path for Hail MatrixTable checkpoints (like the notebook's temp_file_dir)

    Float sv_af_threshold = 0.005   # "common" SV cutoff
    Float snv_af_threshold = 0.005
    Int sv_mac_threshold = 2
    Int window_size = 1000000
    Float strong_ld_rsq = 0.8

    Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7",
                              "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14",
                              "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                              "chr21", "chr22"]

    String hail_docker
    String rdpesr_benchmark_docker

    RuntimeAttr? runtime_attr_prepare_mt
    RuntimeAttr? runtime_attr_export_sv_bed
    RuntimeAttr? runtime_attr_annotate_gc
    RuntimeAttr? runtime_attr_compute_ld
    RuntimeAttr? runtime_attr_combine_plot
  }

  call PrepareMatrixTables {
    input:
      snv_vcf = snv_vcf,
      sv_vcf = sv_vcf,
      sample_pop_map = sample_pop_map,
      staging_gcs_dir = staging_gcs_dir,
      hail_docker = hail_docker,
      runtime_attr_override = runtime_attr_prepare_mt
  }

  Array[String] populations = read_lines(PrepareMatrixTables.population_list)

  call ExportSvAnnotationBed {
    input:
      staging_gcs_dir = staging_gcs_dir,
      done_marker = PrepareMatrixTables.done_marker,
      hail_docker = hail_docker,
      runtime_attr_override = runtime_attr_export_sv_bed
  }

  call AnnotateGenomicContext {
    input:
      sv_sites_bed = ExportSvAnnotationBed.sv_sites_bed,
      extra_annotations_tsv = ExportSvAnnotationBed.extra_annotations_tsv,
      repeatmasker_bed = repeatmasker_bed,
      segdup_bed = segdup_bed,
      simplerepeat_bed = simplerepeat_bed,
      annotate_genomic_context_script = annotate_genomic_context_script,
      annotate_genomic_context_helper_r = annotate_genomic_context_helper_r,
      rdpesr_benchmark_docker = rdpesr_benchmark_docker,
      runtime_attr_override = runtime_attr_annotate_gc
  }

  scatter (contig in contigs) {
    call ComputeLdPerContig {
      input:
        contig = contig,
        staging_gcs_dir = staging_gcs_dir,
        populations = populations,
        sv_af_threshold = sv_af_threshold,
        snv_af_threshold = snv_af_threshold,
        sv_mac_threshold = sv_mac_threshold,
        window_size = window_size,
        strong_ld_rsq = strong_ld_rsq,
        done_marker = PrepareMatrixTables.done_marker,
        hail_docker = hail_docker,
        runtime_attr_override = runtime_attr_compute_ld
    }
  }

  call CombineAndPlot {
    input:
      max_ld_reports = flatten(ComputeLdPerContig.max_ld_reports),
      sv_annotations_with_gc = AnnotateGenomicContext.sv_annotations_with_gc,
      output_prefix = output_prefix,
      rdpesr_benchmark_docker = rdpesr_benchmark_docker,
      runtime_attr_override = runtime_attr_combine_plot
  }

  output {
    File max_ld_annotated_tsv = CombineAndPlot.max_ld_annotated_tsv
    Array[File] strong_ld_reports = flatten(ComputeLdPerContig.strong_ld_reports)
    File violin_plot_pdf = CombineAndPlot.violin_plot_pdf
    File sv_annotations_with_gc = AnnotateGenomicContext.sv_annotations_with_gc
  }
}

task PrepareMatrixTables {
  input {
    File snv_vcf
    File sv_vcf
    File sample_pop_map
    String staging_gcs_dir
    String hail_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 8,
    mem_gb: 52,
    disk_gb: ceil(200 + size(snv_vcf, "GB") * 3 + size(sv_vcf, "GB") * 3),
    boot_disk_gb: 20,
    preemptible_tries: 0,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python <<CODE
import hail as hl

hl.init()


def get_info_or_missing(mt, field, dtype):
    if field in mt.row.info.dtype.fields:
        return mt.info[field]
    return hl.missing(dtype)


# SV VCF: infer SVTYPE/END/allele_length for TRGT tandem-repeat records that
# lack standard SV INFO fields (allele_type == "trv"), same logic as the
# source notebook. Harmless no-op if the SV VCF has no TRGT records.
sv_mt = hl.import_vcf("~{sv_vcf}", reference_genome="GRCh38",
                       array_elements_required=False, force_bgz=True)

sv_mt = sv_mt.annotate_rows(
    SVTYPE=get_info_or_missing(sv_mt, "SVTYPE", hl.tstr),
    END=get_info_or_missing(sv_mt, "END", hl.tint32),
    SOURCE=get_info_or_missing(sv_mt, "SOURCE", hl.tstr),
    SVLEN=get_info_or_missing(sv_mt, "SVLEN", hl.tint32),
    allele_length=get_info_or_missing(sv_mt, "allele_length", hl.tint32),
    allele_type=get_info_or_missing(sv_mt, "allele_type", hl.tstr),
    MOTIFS=get_info_or_missing(sv_mt, "MOTIFS", hl.tstr),
)
sv_mt = sv_mt.drop("info")
sv_mt = sv_mt.select_entries("GT")

sv_mt = sv_mt.annotate_rows(
    allele_length=hl.if_else(
        sv_mt.allele_type == "trv",
        sv_mt.alleles[1].length() - sv_mt.alleles[0].length(),
        sv_mt.allele_length))
sv_mt = sv_mt.annotate_rows(
    SVTYPE=hl.if_else(
        hl.is_missing(sv_mt.SVTYPE),
        hl.if_else(
            sv_mt.allele_type == "trv",
            hl.if_else(sv_mt.allele_length < 0, "DEL",
                       hl.if_else(sv_mt.allele_length > 0, "INS", "BAL")),
            hl.if_else(hl.is_missing(sv_mt.allele_type), "SV", sv_mt.allele_type.upper())),
        sv_mt.SVTYPE))
sv_mt = sv_mt.annotate_rows(
    END=hl.if_else(
        hl.is_missing(sv_mt.END),
        hl.if_else(sv_mt.SVTYPE == "INS", sv_mt.locus.position + 1,
                   hl.if_else(sv_mt.SVTYPE == "DEL", sv_mt.locus.position - sv_mt.allele_length,
                              sv_mt.locus.position + sv_mt.alleles[0].length() - 1)),
        sv_mt.END))

# SNV VCF: every row is already an SNV by construction (dedicated SNV callset).
snv_mt = hl.import_vcf("~{snv_vcf}", reference_genome="GRCh38",
                        array_elements_required=False, force_bgz=True)
snv_mt = snv_mt.select_entries("GT")
snv_mt = snv_mt.annotate_rows(SVTYPE="SNV", END=snv_mt.locus.position)

# Harmonize samples: intersect and reorder columns so SV/SNV matrices align.
sv_samples = sv_mt.s.collect()
snv_samples = snv_mt.s.collect()
sv_sample_set = set(sv_samples)

common_order = [s for s in snv_samples if s in sv_sample_set]
new_sv_order = [sv_samples.index(s) for s in common_order]
new_snv_order = [snv_samples.index(s) for s in common_order]

sv_mt = sv_mt.choose_cols(new_sv_order)
snv_mt = snv_mt.choose_cols(new_snv_order)

pop_map = {}
with open("~{sample_pop_map}") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2 and parts[0]:
            pop_map[parts[0]] = parts[1]

pop_dict = hl.literal(pop_map, dtype=hl.tdict(hl.tstr, hl.tstr))
sv_mt = sv_mt.annotate_cols(POP=pop_dict.get(sv_mt.s, "UNK"))
snv_mt = snv_mt.annotate_cols(POP=pop_dict.get(snv_mt.s, "UNK"))

sv_mt = sv_mt.checkpoint("~{staging_gcs_dir}/svs.mt", overwrite=True)
snv_mt = snv_mt.checkpoint("~{staging_gcs_dir}/snps.mt", overwrite=True)

pops = sorted(set(pop_map.get(s, "UNK") for s in common_order))
with open("population_list.txt", "w") as f:
    f.write("ALL\n")
    for p in pops:
        f.write(p + "\n")

with open("done_marker.txt", "w") as f:
    f.write("done\n")
CODE
  >>>

  output {
    File done_marker = "done_marker.txt"
    File population_list = "population_list.txt"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: hail_docker
  }
}

task ExportSvAnnotationBed {
  input {
    String staging_gcs_dir
    File done_marker
    String hail_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 4,
    mem_gb: 26,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python <<CODE
import hail as hl

hl.init()

sv_mt = hl.read_matrix_table("~{staging_gcs_dir}/svs.mt")
sv_mt = hl.variant_qc(sv_mt)

chroms = sv_mt.locus.contig.collect()
poss = sv_mt.locus.position.collect()
ends = sv_mt.END.collect()
orig_ids = sv_mt.rsid.collect()
sources = sv_mt.SOURCE.collect()
svlens = sv_mt.SVLEN.collect()
alens = sv_mt.allele_length.collect()
svtypes = sv_mt.SVTYPE.collect()
atypes = sv_mt.allele_type.collect()
motifs = sv_mt.MOTIFS.collect()
acs = sv_mt.variant_qc.AC.collect()
ans = sv_mt.variant_qc.AN.collect()
afs = sv_mt.variant_qc.AF.collect()

with open("extra_annotations.tsv", "w") as f, open("sv_sites.bed", "w") as bed:
    f.write("\t".join(["sv_uid", "chrom", "pos", "end", "orig_id", "source", "svlen",
                        "allele_length", "svtype", "allele_type", "motifs",
                        "af", "ac", "an"]) + "\n")
    for j in range(len(chroms)):
        chrom = str(chroms[j])
        pos = poss[j]
        end = ends[j] if ends[j] is not None else pos
        svtype = str(svtypes[j])
        sv_uid = "{}:{}-{}_{}".format(chrom, pos, end, svtype)
        af = afs[j][1] if afs[j] is not None else ""
        ac = acs[j][1] if acs[j] is not None else ""
        an = ans[j] if ans[j] is not None else ""
        f.write("\t".join([sv_uid, chrom, str(pos), str(end), str(orig_ids[j]),
                            str(sources[j]), str(svlens[j]), str(alens[j]),
                            svtype, str(atypes[j]), str(motifs[j]),
                            str(af), str(ac), str(an)]) + "\n")

        start = max(pos - 1, 0)
        svlen = svlens[j] if svlens[j] is not None else (alens[j] if alens[j] is not None else 0)
        bed.write("\t".join([chrom, str(start), str(end), sv_uid, svtype, str(svlen)]) + "\n")
CODE
  >>>

  output {
    File sv_sites_bed = "sv_sites.bed"
    File extra_annotations_tsv = "extra_annotations.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: hail_docker
  }
}

task AnnotateGenomicContext {
  input {
    File sv_sites_bed
    File extra_annotations_tsv
    File repeatmasker_bed
    File segdup_bed
    File simplerepeat_bed
    File annotate_genomic_context_script
    File annotate_genomic_context_helper_r
    String rdpesr_benchmark_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 8,
    disk_gb: ceil(50 + size(repeatmasker_bed, "GB") * 6),
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    # annotate_genomic_context.sh has a bgzip-detection bug when optional
    # flags are passed (it checks $1 after `shift 2`, which by then is the
    # first flag, not INPUT). Always feed it uncompressed inputs so that
    # broken detection path never matters.
    gunzip -c ~{repeatmasker_bed} > RM.bed
    gunzip -c ~{segdup_bed} > SD.bed
    gunzip -c ~{simplerepeat_bed} > SR.bed

    bash ~{annotate_genomic_context_script} ~{sv_sites_bed} sv_gc.tsv \
      --rm RM.bed --sd SD.bed --sr SR.bed

    python3 <<CODE
gc = {}
with open("sv_gc.tsv") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if len(parts) >= 2:
            gc[parts[0]] = parts[1]

with open("~{extra_annotations_tsv}") as fin, open("sv_annotations_with_gc.tsv", "w") as fout:
    header = fin.readline().rstrip("\n")
    fout.write(header + "\tgenomic_context\n")
    for line in fin:
        parts = line.rstrip("\n").split("\t")
        sv_uid = parts[0]
        fout.write(line.rstrip("\n") + "\t" + gc.get(sv_uid, "US") + "\n")
CODE
  >>>

  output {
    File sv_annotations_with_gc = "sv_annotations_with_gc.tsv"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: rdpesr_benchmark_docker
  }
}

task ComputeLdPerContig {
  input {
    String contig
    String staging_gcs_dir
    Array[String] populations
    Float sv_af_threshold
    Float snv_af_threshold
    Int sv_mac_threshold
    Int window_size
    Float strong_ld_rsq
    File done_marker
    String hail_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 4,
    mem_gb: 26,
    disk_gb: 100,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python <<CODE
import hail as hl

hl.init()

contig = "~{contig}"
window_size = ~{window_size}
sv_af_threshold = ~{sv_af_threshold}
snv_af_threshold = ~{snv_af_threshold}
sv_mac_threshold = ~{sv_mac_threshold}
strong_ld_thresh = ~{strong_ld_rsq}
populations = "~{sep=',' populations}".split(",")
staging_dir = "~{staging_gcs_dir}"

grch38_lengths = hl.get_reference("GRCh38").lengths
contig_length = grch38_lengths[contig]

snv_mt_contig = hl.read_matrix_table(staging_dir + "/snps.mt")
snv_mt_contig = snv_mt_contig.filter_rows(snv_mt_contig.locus.contig == contig)
sv_mt_contig = hl.read_matrix_table(staging_dir + "/svs.mt")
sv_mt_contig = sv_mt_contig.filter_rows(sv_mt_contig.locus.contig == contig)


def find_ranges(sv_mins, sv_maxes, snp_positions):
    # Two-pointer sweep: for each SV window [sv_mins[k], sv_maxes[k]), find
    # the [range_mins[k], range_maxes[k]] index range into snp_positions
    # (which must be sorted ascending) that falls inside the window.
    n_sv = len(sv_mins)
    range_mins = [0] * n_sv
    range_maxes = [0] * n_sv
    open_intervals = set()
    added_idx = 0
    for i, pos in enumerate(snp_positions):
        while added_idx < n_sv and sv_maxes[added_idx] < pos:
            added_idx += 1
        while added_idx < n_sv and sv_mins[added_idx] <= pos and sv_maxes[added_idx] > pos:
            open_intervals.add(added_idx)
            range_mins[added_idx] = i
            added_idx += 1
        for sv_idx in list(open_intervals):
            if sv_maxes[sv_idx] < pos:
                open_intervals.discard(sv_idx)
            else:
                range_maxes[sv_idx] = i
    return range_mins, range_maxes


def process_ld(snv_mt_in, sv_mt_in, pop):
    if pop == "ALL":
        snv_mt = snv_mt_in
        sv_mt = sv_mt_in
    else:
        snv_mt = snv_mt_in.filter_cols(snv_mt_in.POP == pop)
        sv_mt = sv_mt_in.filter_cols(sv_mt_in.POP == pop)

    snv_mt = hl.variant_qc(snv_mt)
    snv_mt = snv_mt.filter_rows(hl.min(snv_mt.variant_qc.AF) >= snv_af_threshold)
    snv_mt = snv_mt.checkpoint("snv_{}_{}.mt".format(pop, contig), overwrite=True)

    sv_mt = hl.variant_qc(sv_mt)
    sv_mt = sv_mt.filter_rows(hl.min(sv_mt.variant_qc.AF) >= sv_af_threshold)
    sv_mt = sv_mt.filter_rows(hl.min(sv_mt.variant_qc.AC) >= sv_mac_threshold)
    sv_mt = sv_mt.checkpoint("sv_{}_{}.mt".format(pop, contig), overwrite=True)

    n_sv = sv_mt.count_rows()
    n_snv = snv_mt.count_rows()
    if n_sv == 0 or n_snv == 0:
        return

    sv_pos = sv_mt.locus.position.collect()
    sv_end = sv_mt.END.collect()
    sv_rsid = sv_mt.rsid.collect()
    sv_svtype = sv_mt.SVTYPE.collect()
    sv_filters = sv_mt.filters.collect()

    snp_pos = snv_mt.locus.position.collect()
    snp_alleles = snv_mt.alleles.collect()

    sv_mins = [max(0, p - window_size) for p in sv_pos]
    sv_maxes = [min(contig_length, (e if e is not None else p) + window_size)
                for p, e in zip(sv_pos, sv_end)]
    range_mins, range_maxes = find_ranges(sv_mins, sv_maxes, snp_pos)

    sv_bm = hl.linalg.BlockMatrix.from_entry_expr(
        sv_mt.GT.n_alt_alleles(), mean_impute=True, center=True,
        axis="rows", normalize=True, block_size=2048)
    snv_bm = hl.linalg.BlockMatrix.from_entry_expr(
        snv_mt.GT.n_alt_alleles(), mean_impute=True, center=True,
        axis="rows", normalize=True, block_size=2048)

    cor_mat = sv_bm @ snv_bm.T
    cor_mat = cor_mat.sparsify_row_intervals(range_mins, range_maxes)
    rsq_mat = cor_mat ** 2
    rsq_entries = rsq_mat.entries()

    maxes = rsq_entries.aggregate(hl.agg.group_by(rsq_entries.i, hl.agg.max(rsq_entries.entry)))
    strong = rsq_entries.filter(rsq_entries.entry > strong_ld_thresh)
    strong_grouped = strong.aggregate(
        hl.agg.group_by(strong.i, hl.agg.collect_as_set(
            hl.struct(j=strong.j, rsq=strong.entry))))

    with open("{}_{}_max_ld.tsv".format(pop, contig), "w") as f:
        f.write("\t".join(["sv_uid", "contig", "pos", "end", "orig_id", "svtype",
                            "filtered", "pop", "rsq"]) + "\n")
        for i in range(n_sv):
            if i not in maxes:
                continue
            end = sv_end[i] if sv_end[i] is not None else sv_pos[i]
            sv_uid = "{}:{}-{}_{}".format(contig, sv_pos[i], end, sv_svtype[i])
            f.write("\t".join([sv_uid, contig, str(sv_pos[i]), str(end),
                                str(sv_rsid[i]), str(sv_svtype[i]),
                                str(len(sv_filters[i]) > 0), pop, str(maxes[i])]) + "\n")

    with open("{}_{}_strong_ld.tsv".format(pop, contig), "w") as f:
        f.write("\t".join(["sv_uid", "contig", "pos", "end", "orig_id", "svtype",
                            "pop", "snv_chrom", "snv_pos", "snv_ref", "snv_alt", "rsq"]) + "\n")
        for i, links in strong_grouped.items():
            end = sv_end[i] if sv_end[i] is not None else sv_pos[i]
            sv_uid = "{}:{}-{}_{}".format(contig, sv_pos[i], end, sv_svtype[i])
            for link in links:
                alleles = snp_alleles[link.j]
                f.write("\t".join([sv_uid, contig, str(sv_pos[i]), str(end),
                                    str(sv_rsid[i]), str(sv_svtype[i]), pop,
                                    contig, str(snp_pos[link.j]), str(alleles[0]),
                                    str(alleles[1]), str(link.rsq)]) + "\n")


for pop in populations:
    process_ld(snv_mt_contig, sv_mt_contig, pop)
CODE
  >>>

  output {
    Array[File] max_ld_reports = glob("*_max_ld.tsv")
    Array[File] strong_ld_reports = glob("*_strong_ld.tsv")
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: hail_docker
  }
}

task CombineAndPlot {
  input {
    Array[File] max_ld_reports
    File sv_annotations_with_gc
    String output_prefix
    String rdpesr_benchmark_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    cpu_cores: 2,
    mem_gb: 8,
    disk_gb: 50,
    boot_disk_gb: 10,
    preemptible_tries: 1,
    max_retries: 1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 <<CODE
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["font.family"] = "Arial"

max_files = "~{sep=',' max_ld_reports}".split(",")
frames = [pd.read_csv(f, sep="\t") for f in max_files if f]
maxes = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
    columns=["sv_uid", "contig", "pos", "end", "orig_id", "svtype", "filtered", "pop", "rsq"])

annot = pd.read_csv("~{sv_annotations_with_gc}", sep="\t")[["sv_uid", "genomic_context", "af"]]

merged = maxes.merge(annot, on="sv_uid", how="left")
merged.to_csv("~{output_prefix}.max_ld_annotated.tsv", sep="\t", index=False)

merged["filtered"] = merged["filtered"].astype(str) == "True"
plot_df = merged[~merged["filtered"]]

pops = sorted(plot_df["pop"].unique()) if len(plot_df) else ["ALL"]

fig, axes = plt.subplots(1, len(pops), figsize=(8 * len(pops), 6), squeeze=False)
for ax, pop in zip(axes[0], pops):
    sub = plot_df[plot_df["pop"] == pop]
    svtypes = sorted(sub["svtype"].unique())
    data = [sub.loc[sub["svtype"] == st, "rsq"].dropna().values for st in svtypes]
    if any(len(d) > 0 for d in data):
        ax.violinplot(data, showmedians=True)
        ax.set_xticks(range(1, len(svtypes) + 1))
        ax.set_xticklabels(svtypes, rotation=45, ha="right")
    ax.set_title("Max SNP r-squared by SVTYPE ({})".format(pop))
    ax.set_ylabel("r-squared")

plt.tight_layout()
plt.savefig("~{output_prefix}.maxrsq_by_svtype.pdf")
CODE
  >>>

  output {
    File max_ld_annotated_tsv = "~{output_prefix}.max_ld_annotated.tsv"
    File violin_plot_pdf = "~{output_prefix}.maxrsq_by_svtype.pdf"
  }

  runtime {
    cpu: select_first([runtime_attr.cpu_cores, default_attr.cpu_cores])
    memory: select_first([runtime_attr.mem_gb, default_attr.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " SSD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, default_attr.boot_disk_gb])
    preemptible: select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries: select_first([runtime_attr.max_retries, default_attr.max_retries])
    docker: rdpesr_benchmark_docker
  }
}
