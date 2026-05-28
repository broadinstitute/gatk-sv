#!/usr/bin/env python3
"""
Figure: Paired-end RNA-seq reads aligned against gene exons.
Shows exon structure, read pairs, junction-spanning reads, and coverage.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np

OUT = "/Users/xzhao/Downloads/github/gatk-sv/module19_LR_benchmark/scripts/paralog_expression_method.png"

# ── palette ────────────────────────────────────────────────────────────────
EXON_CLR   = "#2E75B6"
INTRON_CLR = "#7F7F7F"
READ1_CLR  = "#E07B39"
READ2_CLR  = "#5BA85B"
JUNC_CLR   = "#C00000"
INSERT_CLR = "#CCCCCC"
COV_CLR    = "#BDD7EE"
COV_EDG    = "#2E75B6"
BG         = "#FFFFFF"
DARK       = "#1A1A2E"
LABEL_CLR  = "#444444"

EXONS      = [(5, 15), (28, 40), (52, 65), (78, 95)]
EXON_NAMES = ["Exon 1", "Exon 2", "Exon 3", "Exon 4"]

def exon_ax(x, ax_l=0.08, ax_r=0.92):
    return ax_l + (x - 0) / 100 * (ax_r - ax_l)

def is_junction_spanning(r_start, r_end):
    juncs = []
    for i in range(len(EXONS) - 1):
        ex_end, next_start = EXONS[i][1], EXONS[i+1][0]
        if r_start < ex_end and r_end > next_start:
            juncs.append((ex_end, next_start))
    return juncs

def draw_read(ax, r_start, r_end, row_y, color, rh=0.038, zorder=4):
    juncs = is_junction_spanning(r_start, r_end)
    mid_y = row_y + rh / 2
    if not juncs:
        x0, x1 = exon_ax(r_start), exon_ax(r_end)
        ax.add_patch(mpatches.FancyArrow(
            x0, mid_y, x1 - x0, 0,
            width=rh * 0.65, head_width=rh * 0.8, head_length=0.008,
            fc=color, ec=color, alpha=0.85, length_includes_head=True, zorder=zorder))
    else:
        segments, seg_s = [], r_start
        for (ex_end, next_start) in juncs:
            segments.append((seg_s, min(r_end, ex_end)))
            seg_s = max(r_start, next_start)
        segments.append((seg_s, r_end))
        for k, (s, e) in enumerate(segments):
            if e <= s:
                continue
            x0, x1 = exon_ax(s), exon_ax(e)
            # last segment gets arrow head
            if k == len(segments) - 1:
                ax.add_patch(mpatches.FancyArrow(
                    x0, mid_y, x1 - x0, 0,
                    width=rh * 0.65, head_width=rh * 0.8, head_length=0.008,
                    fc=color, ec=color, alpha=0.85, length_includes_head=True, zorder=zorder))
            else:
                ax.add_patch(FancyBboxPatch(
                    (x0, row_y + rh * 0.1), x1 - x0, rh * 0.8,
                    boxstyle="square,pad=0", linewidth=0,
                    facecolor=color, alpha=0.85, zorder=zorder))
        for (ex_end, next_start) in juncs:
            ax.annotate("", xy=(exon_ax(next_start), mid_y),
                        xytext=(exon_ax(ex_end), mid_y),
                        arrowprops=dict(arrowstyle="-|>", color=JUNC_CLR,
                                        lw=1.5, mutation_scale=8,
                                        connectionstyle="arc3,rad=-0.45"),
                        zorder=zorder + 1)

READ_PAIRS = [
    # (r1_start, r1_end, r2_start, r2_end, label, note)
    (6.0,  11.5, 29.5, 36.0,  "Pair A", "Mates in adjacent exons"),
    (12.0, 28.5, 30.5, 37.5,  "Pair B", "Read 1 spans Exon 1\u20132 junction"),
    (30.5, 36.5, 38.5, 53.5,  "Pair C", "Read 2 spans Exon 2\u20133 junction"),
    (53.5, 62.0, 79.5, 90.0,  "Pair D", "Mates in non-adjacent exons"),
]
READ_ROWS = [0.572, 0.514, 0.455, 0.396]

# ─────────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(14, 9), dpi=150, facecolor=BG)
ax  = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off"); ax.set_facecolor(BG)

ax.text(0.5, 0.97, "Paired-end RNA-seq read alignment to gene exons",
        ha="center", va="top", fontsize=13, fontweight="bold", color=DARK)

for ly in [0.64, 0.29]:
    ax.axhline(ly, xmin=0.02, xmax=0.98, color="#DDDDDD", lw=1)

# ── Panel A: Gene model ────────────────────────────────────────────────────
ax.text(0.03, 0.93, "A   Gene model", fontsize=10, fontweight="bold", color=DARK)
gene_y, gene_h = 0.79, 0.030
intron_y = gene_y + gene_h / 2
ax.plot([exon_ax(0), exon_ax(100)], [intron_y, intron_y],
        color=INTRON_CLR, lw=1.5, zorder=1)
for i in range(len(EXONS) - 1):
    ex_end, ex_next = EXONS[i][1], EXONS[i+1][0]
    for t in range(max(1, int((ex_next - ex_end) / 5))):
        tx = exon_ax(ex_end + (ex_next - ex_end) * (t + 0.5) / max(1, int((ex_next - ex_end)/5)))
        ax.annotate("", xy=(tx + 0.009, intron_y), xytext=(tx, intron_y),
                    arrowprops=dict(arrowstyle="-|>", color=INTRON_CLR,
                                   lw=0.8, mutation_scale=5), zorder=2)
for (es, ee), ename in zip(EXONS, EXON_NAMES):
    x0, x1 = exon_ax(es), exon_ax(ee)
    ax.add_patch(FancyBboxPatch((x0, gene_y), x1-x0, gene_h,
                                boxstyle="square,pad=0", linewidth=1.5,
                                edgecolor=EXON_CLR, facecolor=EXON_CLR, zorder=3))
    ax.text((x0+x1)/2, gene_y + gene_h + 0.012, ename,
            ha="center", va="bottom", fontsize=8.5, color=EXON_CLR, fontweight="bold")
ax.text(exon_ax(100) + 0.01, intron_y, "5′ → 3′", fontsize=8, color=INTRON_CLR, va="center")

# ── Panel B: Read alignment ────────────────────────────────────────────────
ax.text(0.03, 0.628, "B   Paired-end read alignment", fontsize=10, fontweight="bold", color=DARK)
for (es, ee) in EXONS:
    for c in [es, ee]:
        ax.axvline(exon_ax(c), ymin=0.300, ymax=0.630, color=EXON_CLR, lw=0.5, alpha=0.18)
for (r1s, r1e, r2s, r2e, label, note), row_y in zip(READ_PAIRS, READ_ROWS):
    mid_y = row_y + 0.019
    ax.plot([exon_ax(r1e), exon_ax(r2s)], [mid_y, mid_y],
            color=INSERT_CLR, lw=1.2, linestyle="--", zorder=3)
    draw_read(ax, r1s, r1e, row_y, READ1_CLR)
    draw_read(ax, r2s, r2e, row_y, READ2_CLR)
    ax.text(0.015, mid_y, label, va="center", ha="left", fontsize=7.5,
            fontweight="bold", color=LABEL_CLR)
    ax.text(0.08, mid_y - 0.015, note, va="top", ha="left", fontsize=7,
            color=LABEL_CLR, style="italic")

# legend
lx, ly = 0.45, 0.298
for ldx, lc, lt in [(0, READ1_CLR, "Read 1"), (0.12, READ2_CLR, "Read 2")]:
    ax.add_patch(FancyBboxPatch((lx+ldx, ly-0.010), 0.030, 0.019,
                                boxstyle="square,pad=0", facecolor=lc, linewidth=0, alpha=0.85))
    ax.text(lx+ldx+0.035, ly+0.0, lt, va="center", fontsize=8, color=DARK)
ax.plot([lx+0.24, lx+0.27], [ly, ly], color=INSERT_CLR, lw=1.5, linestyle="--")
ax.text(lx+0.275, ly, "Insert", va="center", fontsize=8, color=DARK)
alx = lx + 0.35
ax.annotate("", xy=(alx+0.040, ly), xytext=(alx, ly),
            arrowprops=dict(arrowstyle="-|>", color=JUNC_CLR, lw=1.5,
                            mutation_scale=8, connectionstyle="arc3,rad=-0.5"))
ax.text(alx+0.045, ly, "Spliced junction", va="center", fontsize=8, color=DARK)

# ── Panel C: Per-exon coverage ────────────────────────────────────────────
ax.text(0.03, 0.270, "C   Per-exon coverage", fontsize=10, fontweight="bold", color=DARK)
cov_base, cov_max_h = 0.075, 0.160
exon_depths = [5, 8, 7, 4]
for (es, ee), depth, ename in zip(EXONS, exon_depths, EXON_NAMES):
    x0, x1 = exon_ax(es), exon_ax(ee)
    bh = depth / max(exon_depths) * cov_max_h
    ax.add_patch(FancyBboxPatch((x0+0.003, cov_base), x1-x0-0.006, bh,
                                boxstyle="square,pad=0", linewidth=1.2,
                                edgecolor=COV_EDG, facecolor=COV_CLR, zorder=3))
    ax.text((x0+x1)/2, cov_base-0.013, ename, ha="center", va="top",
            fontsize=8.5, color=EXON_CLR, fontweight="bold")
    ax.text((x0+x1)/2, cov_base+bh+0.007, f"{depth}×",
            ha="center", va="bottom", fontsize=9, fontweight="bold", color=EXON_CLR)
for i in range(len(EXONS)-1):
    ax.plot([exon_ax(EXONS[i][1]), exon_ax(EXONS[i+1][0])], [cov_base, cov_base],
            color=COV_EDG, lw=0.8, alpha=0.4, linestyle=":")
ax.text(0.055, cov_base + cov_max_h/2, "Coverage\n(reads)",
        ha="center", va="center", fontsize=8, color=DARK, rotation=90)

plt.savefig(OUT, dpi=150, bbox_inches="tight", facecolor=BG)
print(f"Saved: {OUT}")
