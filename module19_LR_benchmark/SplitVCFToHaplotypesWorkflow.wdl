version 1.0

import "Structs.wdl"

workflow SplitVCFToHaplotypesWorkflow {
  input {
    File input_vcf
    File input_vcf_idx
    File input_gtf  
    String output_prefix

    String sv_pipeline_base_docker

    RuntimeAttr? runtime_attr_annotate_vcf_with_genes
    RuntimeAttr? runtime_attr_split_vcf_to_haplotypes
    RuntimeAttr? runtime_attr_gene_gt_pattern
  }

  call AnnotateVCFWithGenes{
    input:
      input_vcf = input_vcf,
      input_vcf_idx = input_vcf_idx,
      input_gtf = input_gtf,
      output_prefix = "~{output_prefix}.with_gene_anno",
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_annotate_vcf_with_genes
  }

  call SplitVCFToHaplotypes {
    input:
      input_vcf = AnnotateVCFWithGenes.annotated_vcf,
      input_vcf_idx = AnnotateVCFWithGenes.annotated_vcf_idx,
      output_prefix = "~{output_prefix}.with_gene_anno.haplotypes",
      docker_image = sv_pipeline_base_docker,
      runtime_attr_override = runtime_attr_split_vcf_to_haplotypes
  }

  call PrepareAndIndexGTF {
    input:
        input_gtf_gz = input_gtf,
        docker_image = sv_pipeline_base_docker
  }

  call GeneGTPattern {
    input:
        input_vcf = SplitVCFToHaplotypes.output_vcf,
        input_gtf = PrepareAndIndexGTF.sorted_gtf_gz, 
        input_gtf_idx = PrepareAndIndexGTF.sorted_gtf_tbi, 
        docker_image = sv_pipeline_base_docker,
        runtime_attr_override = runtime_attr_gene_gt_pattern
  }


  output {
    File hap_vcf = SplitVCFToHaplotypes.output_vcf
    File gene_pattern_table = GeneGTPattern.pattern_table
    File gene_pattern_summary = GeneGTPattern.gene_summary
  }
}


task PrepareAndIndexGTF {
    input {
        File input_gtf_gz
        String docker_image
        RuntimeAttr? runtime_attr_override
    }

    String prefix = basename(input_gtf_gz, ".gtf.gz")
    command <<<
        set -euo pipefail

        echo "Decompressing GTF..."
        gunzip -c ~{input_gtf_gz} > ~{prefix}.gtf

        echo "Sorting GTF..."
        sort -k1,1 -k4,4n ~{prefix}.gtf | bgzip > ~{prefix}.sorted.gtf.gz

        echo "Indexing with tabix..."
        tabix -p gff ~{prefix}.sorted.gtf.gz

        echo "Done."
    >>>

    output {
        File sorted_gtf_gz = "~{prefix}.sorted.gtf.gz"
        File sorted_gtf_tbi = "~{prefix}.sorted.gtf.gz.tbi"
    }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(input_gtf_gz,"GiB"))*2,
        disk_gb: 20 + ceil(size(input_gtf_gz,"GiB"))*2,
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



task SplitVCFToHaplotypes {
  input {
    File input_vcf
    File input_vcf_idx
    String output_prefix
    String output_type = "z"  # v, z, or b

    String docker_image

    RuntimeAttr? runtime_attr_override

  }

  command <<<
    set -euo pipefail

    cat << 'EOF' > split_vcf_haplotypes.py
#!/usr/bin/env python3

import argparse
import pysam
import sys
import os

def parse_args():
    parser = argparse.ArgumentParser(
        description="Split phased diploid genotypes into haplotypes"
    )
    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--output-type", choices=["v","z","b"], default="z")
    return parser.parse_args()

def main():
    args = parse_args()

    mode_map = {"v": "w", "z": "wz", "b": "wb"}
    mode = mode_map[args.output_type]

    vcf_in = pysam.VariantFile(args.input)
    old_samples = list(vcf_in.header.samples)

    new_header = pysam.VariantHeader()

    # Copy header records
    for rec in vcf_in.header.records:
        if rec.key != "SAMPLE":
            new_header.add_record(rec)

    # Add haplotype samples
    for s in old_samples:
        new_header.add_sample(f"{s}_hap1")
        new_header.add_sample(f"{s}_hap2")

    vcf_out = pysam.VariantFile(args.output, mode, header=new_header)

    for rec in vcf_in:
        new_rec = vcf_out.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=rec.alleles,
            id=rec.id,
            qual=rec.qual,
            filter=rec.filter.keys()
        )

        for key in rec.info:
            new_rec.info[key] = rec.info[key]

        for s in old_samples:
            gt = rec.samples[s].get("GT")

            if gt is None or None in gt:
                h1, h2 = (None, None)
            else:
                if len(gt) == 2:
                    h1, h2 = gt
                else:
                    h1, h2 = (None, None)

            new_rec.samples[f"{s}_hap1"]["GT"] = (h1,)
            new_rec.samples[f"{s}_hap2"]["GT"] = (h2,)

        vcf_out.write(new_rec)

    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    main()
EOF

    python3 split_vcf_haplotypes.py \
        -i ~{input_vcf} \
        -o ~{output_prefix}.vcf.gz \
        --output-type ~{output_type}
    tabix -p vcf ~{output_prefix}.vcf.gz
  >>>

  output {
    File output_vcf = "~{output_prefix}.vcf.gz"
    File output_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
  }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(input_vcf,"GiB"))*2,
        disk_gb: 20 + ceil(size(input_vcf,"GiB"))*2,
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


task AnnotateVCFWithGenes {
  input {
    File input_vcf
    File input_vcf_idx
    File input_gtf
    String output_prefix

    String output_type = "z"   # v, z, b

    String docker_image
    RuntimeAttr? runtime_attr_override
  }

  command <<<
    set -euo pipefail

    cat << 'EOF' > annotate_vcf_with_genes.py
#!/usr/bin/env python3

import argparse
import pysam
import sys
import gzip
from collections import defaultdict
import os


def parse_args():
    parser = argparse.ArgumentParser(description="Annotate VCF with gene names from GTF")
    parser.add_argument("-v", "--vcf", required=True)
    parser.add_argument("-g", "--gtf", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--output-type", choices=["v","z","b"], default="z")
    return parser.parse_args()


def parse_gtf_attributes(attr_string):
    attrs = {}
    for item in attr_string.strip().split(";"):
        if item.strip() == "":
            continue
        key, value = item.strip().split(" ", 1)
        attrs[key] = value.replace('"', '').strip()
    return attrs


def load_genes(gtf_file):
    genes = defaultdict(list)

    if os.path.exists(gtf_file + ".tbi"):
        tbx = pysam.TabixFile(gtf_file)
        for chrom in tbx.contigs:
            for line in tbx.fetch(chrom):
                if line.startswith("#"):
                    continue
                fields = line.split("\t")
                if fields[2] != "gene":
                    continue
                start = int(fields[3])
                end = int(fields[4])
                attrs = parse_gtf_attributes(fields[8])
                gene_name = attrs.get("gene_name", attrs.get("gene_id", "UNKNOWN"))
                genes[chrom].append((start, end, gene_name))

    elif gtf_file.endswith(".gz"):
        with gzip.open(gtf_file, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] != "gene":
                    continue
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attrs = parse_gtf_attributes(fields[8])
                gene_name = attrs.get("gene_name", attrs.get("gene_id", "UNKNOWN"))
                genes[chrom].append((start, end, gene_name))

    else:
        with open(gtf_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] != "gene":
                    continue
                chrom = fields[0]
                start = int(fields[3])
                end = int(fields[4])
                attrs = parse_gtf_attributes(fields[8])
                gene_name = attrs.get("gene_name", attrs.get("gene_id", "UNKNOWN"))
                genes[chrom].append((start, end, gene_name))

    return genes


def find_hits(chrom, pos, genes):
    if chrom not in genes:
        return []
    return [g for start, end, g in genes[chrom] if start <= pos <= end]


def main():
    args = parse_args()

    mode_map = {"v": "w", "z": "wz", "b": "wb"}
    mode = mode_map[args.output_type]

    genes = load_genes(args.gtf)

    vcf_in = pysam.VariantFile(args.vcf)
    header = vcf_in.header.copy()

    if "GENE" not in header.info:
        header.info.add("GENE", number=".", type="String",
                        description="Gene(s) overlapping variant")

    vcf_out = pysam.VariantFile(args.output, mode, header=header)

    for rec in vcf_in:
        hits = find_hits(rec.contig, rec.pos, genes)

        new_rec = vcf_out.new_record(
            contig=rec.contig,
            start=rec.start,
            stop=rec.stop,
            alleles=rec.alleles,
            id=rec.id,
            qual=rec.qual,
            filter=rec.filter.keys()
        )

        # copy INFO
        for k in rec.info:
            new_rec.info[k] = rec.info[k]

        if hits:
            new_rec.info["GENE"] = ",".join(sorted(set(hits)))

        # copy GT safely only
        for s in rec.samples:
            gt = rec.samples[s].get("GT")
            if gt is None or None in gt:
                new_rec.samples[s]["GT"] = (None, None)
            else:
                new_rec.samples[s]["GT"] = tuple(gt)

        vcf_out.write(new_rec)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    main()
EOF

    python3 annotate_vcf_with_genes.py \
        -v ~{input_vcf} \
        -g ~{input_gtf} \
        -o ~{output_prefix}.vcf.gz \
        --output-type ~{output_type}

    tabix -p vcf ~{output_prefix}.vcf.gz
  >>>

  output {
    File annotated_vcf = "~{output_prefix}.vcf.gz"
    File annotated_vcf_idx = "~{output_prefix}.vcf.gz.tbi"
  }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(input_vcf,"GiB"))*2,
        disk_gb: 20 + ceil(size(input_vcf,"GiB"))*2,
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


task GeneGTPattern {
  input {
    # Inputs
    File input_vcf           # bgzipped VCF (.vcf.gz)
    File input_gtf           # bgzipped GTF (.gtf.gz, tabix-indexed)
    File input_gtf_idx

    String docker_image
    RuntimeAttr? runtime_attr_override
    # Optional runtime resources
    Int cpu = 2
    Int memory_gb = 8
  }

  String output_prefix = basename(input_vcf, ".vcf.gz")
  command <<<
    set -euo pipefail

    cat << 'EOF' > calcu_gene_GT_patterns.py
#!/usr/bin/env python3

import pysam
import argparse
from collections import defaultdict, Counter

parser = argparse.ArgumentParser()
parser.add_argument("--input_vcf", required=True)
parser.add_argument("--input_gtf", required=True)
parser.add_argument("--out_prefix", required=True)
args = parser.parse_args()

# ----------------------------
# GTF parsing
# ----------------------------
gene_size = {}
exon_count = defaultdict(int)

gtf_file = args.input_gtf
if gtf_file.endswith(".gz"):
    gtf = pysam.TabixFile(gtf_file)
    iterator = (line for chrom in gtf.contigs for line in gtf.fetch(chrom))
else:
    iterator = open(gtf_file)

for line in iterator:
    if line.startswith("#"): continue
    fields = line.strip().split("\t")
    feature, start, end, attr = fields[2], int(fields[3]), int(fields[4]), fields[8]
    gene = None
    for item in attr.split(";"):
        if "gene_name" in item:
            gene = item.split('"')[1]; break
    if gene is None: continue
    if feature=="gene": gene_size[gene] = end-start+1
    elif feature=="exon": exon_count[gene] += 1

# ----------------------------
# Stream VCF
# ----------------------------
vcf = pysam.VariantFile(args.input_vcf)
samples = list(vcf.header.samples)

gene_gt = defaultdict(lambda: defaultdict(list))
gene_variant_count = defaultdict(int)

for rec in vcf:
    gene = rec.info.get("GENE")
    if gene is None: continue
    genes = gene if isinstance(gene, tuple) else [gene]

    gt_cache = {}
    for s in samples:
        gt = rec.samples[s].get("GT")
        val = 0 if gt is None or gt[0] is None else int(gt[0])
        gt_cache[s] = val

    for g in genes:
        gene_variant_count[g] += 1
        for s in samples:
            gene_gt[g][s].append(str(gt_cache[s]))

# ----------------------------
# Compute patterns & summary
# ----------------------------
pattern_out = open(f"{args.out_prefix}.patterns.tsv","w")
pattern_out.write("gene\tgenotype_combination\tAC\n")

summary_out = open(f"{args.out_prefix}.summary.tsv","w")
summary_out.write(
    "gene\tgene_size\tn_exons\tn_variants\tmax_pattern_ac\t"
    "n_patterns\tsingletons\tdoubletons\tcommon_patterns\n"
)

for gene, sample_dict in gene_gt.items():
    combos = [";".join(sample_dict[s]) for s in samples]
    counts = Counter(combos)
    max_pattern_ac = max(counts.values()) if counts else 0

    # write pattern table
    for combo, ac in counts.items():
        pattern_out.write(f"{gene}\t{combo}\t{ac}\n")

    # summary stats
    n_patterns = len(counts)
    singletons = sum(1 for x in counts.values() if x==1)
    doubletons = sum(1 for x in counts.values() if x==2)
    common = sum(1 for x in counts.values() if x>10)

    summary_out.write(
        f"{gene}\t"
        f"{gene_size.get(gene,'NA')}\t"
        f"{exon_count.get(gene,'NA')}\t"
        f"{gene_variant_count.get(gene,0)}\t"
        f"{max_pattern_ac}\t"
        f"{n_patterns}\t"
        f"{singletons}\t"
        f"{doubletons}\t"
        f"{common}\n"
    )

pattern_out.close()
summary_out.close()
EOF


    python3 calcu_gene_GT_patterns.py \
      --input_vcf ~{input_vcf} \
      --input_gtf ~{input_gtf} \
      --out_prefix ~{output_prefix}

  >>>

  output {
    File pattern_table = "~{output_prefix}.patterns.tsv"
    File gene_summary = "~{output_prefix}.summary.tsv"
  }

    RuntimeAttr default_attr = object {
        cpu_cores: 4,
        mem_gb: 16 + ceil(size(input_vcf,"GiB"))*2,
        disk_gb: 20 + ceil(size(input_vcf,"GiB"))*2,
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


