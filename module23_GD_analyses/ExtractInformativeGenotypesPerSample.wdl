version 1.0

## Given a VCF (+ index) and a list of sample IDs (one per line):
##   1. Split out each sample's informative variants -- sites where that
##      sample's GT is neither missing nor hom-ref -- into a single-sample
##      VCF.
##   2. Pull POS/GT/LAA/LAD for those sites.
##   3. Plot ref-vs-alt allele depth across the region (all informative
##      sites), and allele balance (AD_smaller / AD_sum) restricted to
##      heterozygous sites, as a mosaicism/allelic-imbalance screen.
##
## Sample subsetting and GT filtering in step 1 are done as two separate
## bcftools passes: combining `-s` and `-i 'GT="alt"'` in a single
## invocation applies the GT expression before subsetting takes effect, so
## it matches sites where ANY of the original samples is non-ref rather
## than the target sample alone.
##
## LAD is local-allele-indexed: LAD[0] is always the REF depth, and LAD[i]
## (i>=1) is the depth for the allele whose global ALT index is LAA[i-1].
## For a heterozygous genotype a/b, the two allele depths are looked up
## through this LAA->LAD map rather than assumed to be REF vs. first ALT,
## since a/b need not include the REF allele (e.g. a multiallelic 13/14).

workflow ExtractInformativeGenotypesPerSample {
  input {
    File vcf
    File vcf_idx
    File sample_list
    String bcftools_docker = "staphb/bcftools:1.19"
    String python_docker = "python:3.11-slim"
  }

  Array[String] samples = read_lines(sample_list)

  scatter (sample in samples) {
    call ExtractSampleVariants {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        sample = sample,
        docker = bcftools_docker,
    }

    call ExtractGenotypeFields {
      input:
        vcf = ExtractSampleVariants.sample_vcf,
        vcf_idx = ExtractSampleVariants.sample_vcf_idx,
        sample = sample,
        docker = bcftools_docker,
    }

    call PlotAlleleDepth {
      input:
        genotype_tsv = ExtractGenotypeFields.genotype_tsv,
        sample = sample,
        docker = python_docker,
    }
  }

  output {
    Array[File] sample_vcfs = ExtractSampleVariants.sample_vcf
    Array[File] sample_vcf_indices = ExtractSampleVariants.sample_vcf_idx
    Array[File] ref_alt_depth_tables = PlotAlleleDepth.ref_alt_depth_tsv
    Array[File] ref_alt_depth_plots = PlotAlleleDepth.ref_alt_depth_png
    Array[File] het_balance_tables = PlotAlleleDepth.het_balance_tsv
    Array[File] het_balance_plots = PlotAlleleDepth.het_balance_png
  }
}

task ExtractSampleVariants {
  input {
    File vcf
    File vcf_idx
    String sample
    String docker
  }

  Int disk_gb = ceil(size(vcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    ln -s ~{vcf} input.vcf.gz
    ln -s ~{vcf_idx} input.vcf.gz.tbi

    bcftools view -s ~{sample} input.vcf.gz -Ou \
      | bcftools view -i 'GT="alt"' -Oz -o ~{sample}.informative.vcf.gz
    bcftools index -t ~{sample}.informative.vcf.gz
  >>>

  output {
    File sample_vcf = "~{sample}.informative.vcf.gz"
    File sample_vcf_idx = "~{sample}.informative.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task ExtractGenotypeFields {
  input {
    File vcf
    File vcf_idx
    String sample
    String docker
  }

  command <<<
    set -euo pipefail

    ln -s ~{vcf} input.vcf.gz
    ln -s ~{vcf_idx} input.vcf.gz.tbi

    bcftools query -f '%POS\t[%GT]\t[%LAA]\t[%LAD]\n' input.vcf.gz > ~{sample}.genotypes.tsv
  >>>

  output {
    File genotype_tsv = "~{sample}.genotypes.tsv"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "2 GiB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}

task PlotAlleleDepth {
  input {
    File genotype_tsv
    String sample
    String docker
  }

  command <<<
    set -euo pipefail
    pip install --quiet --no-cache-dir matplotlib

    python3 <<CODE
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE = "#2a78d6"
ORANGE = "#eb6834"

depth_rows = []
het_rows = []

with open("~{genotype_tsv}") as f:
    for line in f:
        pos, gt, laa, lad = line.rstrip("\n").split("\t")
        if gt in (".", "./.", ".|."):
            continue
        alleles = re.split(r"[/|]", gt)
        if len(alleles) != 2 or "." in alleles:
            continue
        a1, a2 = int(alleles[0]), int(alleles[1])
        if a1 == 0 and a2 == 0:
            continue
        if lad == ".":
            continue
        lad_vals = [int(x) for x in lad.split(",")]
        depth_map = {0: lad_vals[0]}
        if laa != ".":
            laa_vals = [int(x) for x in laa.split(",")]
            for i, global_idx in enumerate(laa_vals):
                depth_map[global_idx] = lad_vals[i + 1]

        ref_depth = depth_map.get(0, 0)
        alt_depth = sum(v for k, v in depth_map.items() if k != 0)
        depth_rows.append((int(pos), gt, ref_depth, alt_depth))

        if a1 != a2 and a1 in depth_map and a2 in depth_map:
            ad_a, ad_b = depth_map[a1], depth_map[a2]
            ad_sum = ad_a + ad_b
            if ad_sum > 0:
                ad_smaller = min(ad_a, ad_b)
                het_rows.append((int(pos), gt, ad_a, ad_b, ad_smaller, ad_sum, ad_smaller / ad_sum))

depth_rows.sort()
het_rows.sort()

with open("~{sample}.ref_alt_depth.tsv", "w") as f:
    f.write("POS\tGT\tREF_DEPTH\tALT_DEPTH\n")
    for r in depth_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

with open("~{sample}.het_balance.tsv", "w") as f:
    f.write("POS\tGT\tAD_A\tAD_B\tAD_SMALLER\tAD_SUM\tBALANCE\n")
    for r in het_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

if depth_rows:
    xs = [r[0] for r in depth_rows]
    fig, ax = plt.subplots(figsize=(12, 5), dpi=150)
    ax.scatter(xs, [r[2] for r in depth_rows], s=16, color=BLUE, label="Ref depth", zorder=3)
    ax.scatter(xs, [r[3] for r in depth_rows], s=16, color=ORANGE, label="Alt depth", zorder=3)
    ax.set_xlabel("Position")
    ax.set_ylabel("Read depth")
    ax.set_title("Ref vs. alt allele depth across informative sites (~{sample})")
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#e3e2dc", linewidth=1, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig("~{sample}.ref_alt_depth.png")
else:
    fig, ax = plt.subplots(figsize=(12, 5), dpi=150)
    ax.set_title(f"No informative sites for ~{sample}")
    fig.savefig("~{sample}.ref_alt_depth.png")

if het_rows:
    xs = [r[0] for r in het_rows]
    ys = [r[6] for r in het_rows]
    fig, ax = plt.subplots(figsize=(12, 5), dpi=150)
    ax.scatter(xs, ys, s=16, color=BLUE, zorder=3)
    ax.axhline(0.5, color="#b8b6ac", linewidth=1, zorder=2)
    ax.text(max(xs), 0.505, "expected balance (0.5)", ha="right", va="bottom", fontsize=9, color="#8a8980")
    ax.set_xlabel("Position")
    ax.set_ylabel("AD_smaller / AD_sum")
    ax.set_ylim(0, 0.55)
    ax.set_title("Allele balance at heterozygous sites (~{sample})")
    ax.grid(axis="y", color="#e3e2dc", linewidth=1, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig("~{sample}.het_balance.png")
else:
    fig, ax = plt.subplots(figsize=(12, 5), dpi=150)
    ax.set_title(f"No heterozygous sites for ~{sample}")
    fig.savefig("~{sample}.het_balance.png")
CODE
  >>>

  output {
    File ref_alt_depth_tsv = "~{sample}.ref_alt_depth.tsv"
    File ref_alt_depth_png = "~{sample}.ref_alt_depth.png"
    File het_balance_tsv = "~{sample}.het_balance.tsv"
    File het_balance_png = "~{sample}.het_balance.png"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}
