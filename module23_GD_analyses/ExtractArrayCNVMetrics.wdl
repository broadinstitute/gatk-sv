version 1.0

## Given one single-sample SNP-array VCF per sample and a region of
## interest, scatter across samples and run extract_array_cnv_metrics.sh
## per sample to pull out copy-number-relevant metrics (GT, BAF, LRR by
## default). Pivot the per-sample long-format tables into one wide matrix
## per metric: rows are sites, columns are samples, "." for a sample
## missing at a site.
##
## extract_array_cnv_metrics.sh handles compression and indexing itself, so
## each input VCF can be plain, gzipped, or bgzipped, indexed or not.
##
## Given a list of germline CNV carriers and (optionally empty) mosaic CNV
## carriers, also select a random subset of n_ref_samples reference
## samples (anyone not in either carrier list) and, for that combined set:
##   - plot a 2-panel BAF (top) / estimated-CN-from-LRR (bottom) figure per
##     sample, with a smoothed line on the CN panel
##   - plot everyone's smoothed CN together in one figure, colored by
##     carrier status (ref / mosaic / germline)

workflow ExtractArrayCNVMetrics {
  input {
    Array[File] vcfs
    String region
    File script
    Boolean extra_fields = false
    String output_basename = "array_cnv_metrics"
    String docker = "staphb/bcftools:1.19"
    String python_docker = "python:3.11-slim"

    Array[String] germline_carriers
    Array[String] mosaic_carriers = []
    Int n_ref_samples = 10
    Int random_seed = 42
    Array[Int] breakpoints = []
    Int smoothing_window = 10
  }

  scatter (vcf in vcfs) {
    call ExtractPerSample {
      input:
        vcf = vcf,
        region = region,
        script = script,
        extra_fields = extra_fields,
        docker = docker,
    }
  }

  call PivotToWideMatrices {
    input:
      tables = ExtractPerSample.metrics_tsv,
      output_basename = output_basename,
      docker = python_docker,
  }

  call PrepareSampleSelection {
    input:
      all_samples_file = write_lines(ExtractPerSample.sample_id),
      all_tables = ExtractPerSample.metrics_tsv,
      germline_file = write_lines(germline_carriers),
      mosaic_file = write_lines(mosaic_carriers),
      n_ref_samples = n_ref_samples,
      random_seed = random_seed,
      docker = python_docker,
  }

  scatter (i in range(length(PrepareSampleSelection.selected_samples))) {
    call PlotSampleBAFCN {
      input:
        metrics_tsv = PrepareSampleSelection.selected_tables[i],
        sample = PrepareSampleSelection.selected_samples[i],
        region = region,
        breakpoints = breakpoints,
        smoothing_window = smoothing_window,
        docker = python_docker,
    }
  }

  call PlotGroupCN {
    input:
      tables = PrepareSampleSelection.selected_tables,
      samples = PrepareSampleSelection.selected_samples,
      groups = PrepareSampleSelection.selected_groups,
      region = region,
      breakpoints = breakpoints,
      smoothing_window = smoothing_window,
      output_basename = output_basename,
      docker = python_docker,
  }

  output {
    Array[File] per_sample_tables = ExtractPerSample.metrics_tsv
    Array[File] wide_matrices = PivotToWideMatrices.wide_matrices
    Array[String] selected_samples = PrepareSampleSelection.selected_samples
    Array[String] selected_groups = PrepareSampleSelection.selected_groups
    Array[File] baf_cn_plots = PlotSampleBAFCN.plot_png
    File group_cn_plot = PlotGroupCN.plot_png
  }
}

task ExtractPerSample {
  input {
    File vcf
    String region
    File script
    Boolean extra_fields
    String docker
  }

  Int disk_gb = ceil(size(vcf, "GB") * 2) + 20

  command <<<
    set -euo pipefail

    sample=$(bcftools query -l ~{vcf} | head -n1)
    echo "$sample" > sample_id.txt
    bash ~{script} ~{vcf} ~{region} "${sample}.array_cnv_metrics.tsv" ~{true="--extra" false="" extra_fields}
  >>>

  output {
    String sample_id = read_string("sample_id.txt")
    File metrics_tsv = glob("*.array_cnv_metrics.tsv")[0]
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 2
  }
}

task PivotToWideMatrices {
  input {
    Array[File] tables
    String output_basename
    String docker
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
import csv

table_files = "~{sep=',' tables}".split(",")
site_key_cols = ["CHROM", "POS", "ID", "REF", "ALT"]

sites = set()
samples = set()
metric_cols = None
data = {}

for path in table_files:
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        if metric_cols is None:
            metric_cols = [c for c in reader.fieldnames if c not in site_key_cols and c != "SAMPLE"]
            for m in metric_cols:
                data[m] = {}
        for row in reader:
            key = tuple(row[c] for c in site_key_cols)
            sites.add(key)
            sample = row["SAMPLE"]
            samples.add(sample)
            for m in metric_cols:
                data[m].setdefault(key, {})[sample] = row[m]

sorted_sites = sorted(sites, key=lambda k: (k[0], int(k[1])))
sorted_samples = sorted(samples)

for m in metric_cols:
    out_path = f"~{output_basename}.{m}.tsv"
    with open(out_path, "w") as out:
        out.write("\t".join(site_key_cols + sorted_samples) + "\n")
        for key in sorted_sites:
            row_vals = list(key) + [data[m].get(key, {}).get(s, ".") for s in sorted_samples]
            out.write("\t".join(row_vals) + "\n")
    print(f"Wrote {out_path}: {len(sorted_sites)} sites x {len(sorted_samples)} samples")
CODE
  >>>

  output {
    Array[File] wide_matrices = glob("~{output_basename}.*.tsv")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 20 HDD"
    preemptible: 2
  }
}

task PrepareSampleSelection {
  input {
    File all_samples_file
    Array[File] all_tables
    File germline_file
    File mosaic_file
    Int n_ref_samples
    Int random_seed
    String docker
  }

  command <<<
    set -euo pipefail

    python3 <<CODE
import random
import shutil

with open("~{all_samples_file}") as f:
    all_samples = [line.strip() for line in f if line.strip()]
with open("~{germline_file}") as f:
    germline = [line.strip() for line in f if line.strip()]
with open("~{mosaic_file}") as f:
    mosaic = [line.strip() for line in f if line.strip()]

all_tables = "~{sep=',' all_tables}".split(",")
sample_to_table = dict(zip(all_samples, all_tables))

carrier_set = set(germline) | set(mosaic)
ref_pool = [s for s in all_samples if s not in carrier_set]

random.seed(~{random_seed})
ref_samples = random.sample(ref_pool, min(~{n_ref_samples}, len(ref_pool)))

selected = (
    [(s, "germline") for s in germline]
    + [(s, "mosaic") for s in mosaic]
    + [(s, "ref") for s in ref_samples]
)

with open("selected_samples.txt", "w") as fs, open("selected_groups.txt", "w") as fg:
    for i, (sample, group) in enumerate(selected):
        shutil.copy(sample_to_table[sample], f"{i:04d}.{sample}.metrics.tsv")
        fs.write(sample + "\n")
        fg.write(group + "\n")
CODE
  >>>

  output {
    Array[String] selected_samples = read_lines("selected_samples.txt")
    Array[String] selected_groups = read_lines("selected_groups.txt")
    Array[File] selected_tables = glob("*.metrics.tsv")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 20 HDD"
    preemptible: 2
  }
}

task PlotSampleBAFCN {
  input {
    File metrics_tsv
    String sample
    String region
    Array[Int] breakpoints
    Int smoothing_window
    String docker
  }

  command <<<
    set -euo pipefail
    pip install --quiet --no-cache-dir matplotlib numpy

    python3 <<CODE
import csv
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

region = "~{region}"
core_start, core_end = (int(x) for x in re.match(r"chr\w+:(\d+)-(\d+)", region).groups())
breakpoints = [~{sep=',' breakpoints}]
window = ~{smoothing_window}

BLUE = "#2a78d6"
ORANGE = "#eb6834"
GREY = "#8a8980"


def zone(x):
    return "core" if core_start <= x <= core_end else "flank"


def rolling_mean(arr, w):
    return np.array([arr[max(0, i - w // 2): i + w // 2 + 1].mean() for i in range(len(arr))])


def split_by_zone(xs, ys, zones):
    segs = []
    cz, cx, cy = zones[0], [xs[0]], [ys[0]]
    for x, y, z in zip(xs[1:], ys[1:], zones[1:]):
        if z != cz:
            cx.append(x)
            cy.append(y)
            segs.append((cz, cx, cy))
            cz, cx, cy = z, [x], [y]
        else:
            cx.append(x)
            cy.append(y)
    segs.append((cz, cx, cy))
    return segs


rows = []
with open("~{metrics_tsv}") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for r in reader:
        if r["BAF"] == "." or r["LRR"] == ".":
            continue
        pos = int(r["POS"])
        baf = float(r["BAF"])
        cn = 2 * (2 ** float(r["LRR"]))
        rows.append((pos, baf, cn, zone(pos)))
rows.sort()

xs = [r[0] for r in rows]
zones = [r[3] for r in rows]
cn_smoothed = rolling_mean(np.array([r[2] for r in rows]), window)

fig, (ax_baf, ax_cn) = plt.subplots(2, 1, figsize=(13, 7), dpi=150, sharex=True)

for z in ("flank", "core"):
    color = BLUE if z == "core" else ORANGE
    zx = [r[0] for r in rows if r[3] == z]
    zbaf = [r[1] for r in rows if r[3] == z]
    zcn = [r[2] for r in rows if r[3] == z]
    ax_baf.scatter(zx, zbaf, s=6, color=color, alpha=0.6, zorder=3)
    ax_cn.scatter(zx, zcn, s=6, color=color, alpha=0.35, zorder=2)

for z, sx, sy in split_by_zone(xs, list(cn_smoothed), zones):
    ax_cn.plot(sx, sy, color=BLUE if z == "core" else ORANGE, linewidth=2, zorder=3)

for bp in breakpoints:
    ax_baf.axvline(bp, color=GREY, linewidth=1, zorder=2)
    ax_cn.axvline(bp, color=GREY, linewidth=1, zorder=2)
ax_cn.axhline(2, color=GREY, linewidth=1, alpha=0.5, zorder=1)

ax_baf.set_ylabel("BAF")
ax_baf.set_ylim(-0.05, 1.05)
ax_baf.grid(axis="y", color="#e3e2dc", linewidth=1, zorder=0)
ax_baf.spines[["top", "right"]].set_visible(False)
ax_baf.set_title(f"BAF and estimated copy number across {region} — sample ~{sample}")

ax_cn.set_xlabel(f"Position on {region.split(':')[0]}")
ax_cn.set_ylabel("Estimated CN (2×2^LRR)")
ax_cn.set_ylim(0, 3)
ax_cn.set_xlim(min(xs), max(xs))
ax_cn.grid(axis="y", color="#e3e2dc", linewidth=1, zorder=0)
ax_cn.spines[["top", "right"]].set_visible(False)

fig.tight_layout()
fig.savefig("~{sample}.BAF_CN.png")
CODE
  >>>

  output {
    File plot_png = "~{sample}.BAF_CN.png"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}

task PlotGroupCN {
  input {
    Array[File] tables
    Array[String] samples
    Array[String] groups
    String region
    Array[Int] breakpoints
    Int smoothing_window
    String output_basename
    String docker
  }

  command <<<
    set -euo pipefail
    pip install --quiet --no-cache-dir matplotlib numpy

    python3 <<CODE
import csv
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

region = "~{region}"
core_start, core_end = (int(x) for x in re.match(r"chr\w+:(\d+)-(\d+)", region).groups())
breakpoints = [~{sep=',' breakpoints}]
window = ~{smoothing_window}

tables = "~{sep=',' tables}".split(",")
samples = "~{sep=',' samples}".split(",")
groups = "~{sep=',' groups}".split(",")

GREY = "#8a8980"
LIGHTBLUE = "#7fb2e8"
BLUE = "#2a78d6"
CORE_SHADE = "#eef3fa"
GROUP_COLOR = {"ref": GREY, "mosaic": LIGHTBLUE, "germline": BLUE}
GROUP_LABEL = {"ref": "Ref", "mosaic": "Mosaic", "germline": "Germline deletion"}
GROUP_ZORDER = {"ref": 2, "mosaic": 3, "germline": 4}


def rolling_mean(arr, w):
    return np.array([arr[max(0, i - w // 2): i + w // 2 + 1].mean() for i in range(len(arr))])


fig, ax = plt.subplots(figsize=(13, 5), dpi=150)
all_x = []
counts = {"ref": 0, "mosaic": 0, "germline": 0}

for table, sample, group in zip(tables, samples, groups):
    rows = []
    with open(table) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            if r["LRR"] == ".":
                continue
            rows.append((int(r["POS"]), 2 * (2 ** float(r["LRR"]))))
    rows.sort()
    xs = [r[0] for r in rows]
    all_x.extend(xs)
    cn_smoothed = rolling_mean(np.array([r[1] for r in rows]), window)
    counts[group] += 1
    ax.plot(xs, cn_smoothed, color=GROUP_COLOR[group], linewidth=0.6, alpha=0.5, zorder=GROUP_ZORDER[group])

ax.axvspan(core_start, core_end, color=CORE_SHADE, zorder=0)
for bp in breakpoints:
    ax.axvline(bp, color="#52514e", linewidth=1, zorder=5)
ax.axhline(2, color="#52514e", linewidth=1, alpha=0.4, zorder=1)

ax.set_xlabel(f"Position on {region.split(':')[0]}")
ax.set_ylabel("Estimated CN (2×2^LRR)")
ax.set_ylim(0, 3)
ax.set_xlim(min(all_x), max(all_x))
ax.grid(axis="y", color="#e3e2dc", linewidth=1, zorder=0)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title(
    f"Estimated CN, smoothed ({window}-probe rolling mean) per sample across {region} — "
    f"ref (n={counts['ref']}), mosaic (n={counts['mosaic']}), germline (n={counts['germline']})"
)

legend_handles = [
    Line2D([0], [0], color=GROUP_COLOR[g], linewidth=2, label=GROUP_LABEL[g])
    for g in ("ref", "mosaic", "germline")
]
legend_handles.append(Line2D([0], [0], color="#52514e", linewidth=1, label="Breakpoint"))
ax.legend(handles=legend_handles, frameon=False, loc="upper right")

fig.tight_layout()
fig.savefig("~{output_basename}.group_CN.png")
CODE
  >>>

  output {
    File plot_png = "~{output_basename}.group_CN.png"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 10 HDD"
    preemptible: 2
  }
}
