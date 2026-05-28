#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Load data ────────────────────────────────────────────────────────────────
DATA = ("/Users/xzhao/zxf0419@gmail.com - Google Drive/My Drive/"
        "Talkowski_Lab/gnomAD_LR/annotation/aou_I.counts_sites.tsv")
df = pd.read_csv(DATA, sep='\t')

# ── Extract target rows ──────────────────────────────────────────────────────
VAR_COLS = ['SNV', 'INS 1-49bp', 'DEL 1-49bp',
            'INS 50-499bp', 'DEL 50-499bp',
            'INS >499bp', 'DEL >499bp']

# All row: all four key columns == "All"
row_all = df[
    (df['category'] == 'All') &
    (df['sub_category'] == 'All') &
    (df['tr_status'] == 'All') &
    (df['region'] == 'All')
].iloc[0]

# TR row: category="All", sub_category="All", tr_status="TR", region="All"
row_tr = df[
    (df['category'] == 'All') &
    (df['sub_category'] == 'All') &
    (df['tr_status'] == 'TR') &
    (df['region'] == 'All')
].iloc[0]

counts_all = row_all[VAR_COLS].values.astype(float)
counts_tr  = row_tr[VAR_COLS].values.astype(float)

# ── Colors ───────────────────────────────────────────────────────────────────
COLORS = {
    'SNV':         '#7B2D8B',   # purple
    'INS 1-49bp':  '#FFB6C1',   # pink (light)
    'DEL 1-49bp':  '#CC0000',   # red
    'INS 50-499bp':'#FF69B4',   # hot pink
    'DEL 50-499bp':'#FF4444',   # lighter red
    'INS >499bp':  '#FF1493',   # deep pink
    'DEL >499bp':  '#FF6666',   # salmon red
}
bar_colors = [COLORS[c] for c in VAR_COLS]

# ── Compute TR strip extents (in log10-space) ────────────────────────────────
# Strip height in log10-space = log10(count_all) * (count_tr / count_all)
# Strip sits on top of bar:
#   bar top  in log10-space : log10(count_all)
#   strip top in log10-space: log10(count_all) + log10(count_all) * proportion
#                           = log10(count_all) * (1 + proportion)
log10_all = np.log10(counts_all)
proportions = counts_tr / counts_all
strip_log_top  = log10_all * (1 + proportions)   # in log10 units
strip_data_top = 10 ** strip_log_top             # back to data units

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(11, 7))

x = np.arange(len(VAR_COLS))
width = 0.6

bars = ax.bar(x, counts_all, width=width, color=bar_colors,
              edgecolor='white', linewidth=0.8, zorder=3)

# TR proportion strips (darker shade, hatched)
for i, (col, c_all, c_tr, s_top) in enumerate(
        zip(VAR_COLS, counts_all, counts_tr, strip_data_top)):
    strip_bottom = c_all
    strip_height = s_top - strip_bottom
    # choose a darker version of the bar colour for the strip
    darker = {
        'SNV':         '#4B006E',
        'INS 1-49bp':  '#E75480',
        'DEL 1-49bp':  '#660000',
        'INS 50-499bp':'#C71585',
        'DEL 50-499bp':'#990000',
        'INS >499bp':  '#8B0057',
        'DEL >499bp':  '#993333',
    }[col]
    ax.bar(i, strip_height, width=width,
           bottom=strip_bottom, color=darker,
           edgecolor='white', linewidth=0.5,
           alpha=0.85, zorder=4,
           hatch='///')
    # proportion label inside strip (centred vertically in log-space)
    label_y = 10 ** ((np.log10(strip_bottom) + np.log10(s_top)) / 2)
    pct = proportions[i] * 100
    ax.text(i, label_y, f'{pct:.1f}%',
            ha='center', va='center', fontsize=7.5,
            color='white', fontweight='bold', zorder=5)

# Count labels above each bar (above the strip)
for i, (c_all, s_top) in enumerate(zip(counts_all, strip_data_top)):
    ax.text(i, s_top * 1.15,
            f'{int(c_all):,}',
            ha='center', va='bottom', fontsize=8.5,
            fontweight='bold', color='#333333', zorder=5)

# ── Axes formatting ──────────────────────────────────────────────────────────
ax.set_yscale('log')
ax.set_xticks(x)
ax.set_xticklabels(VAR_COLS, rotation=25, ha='right', fontsize=11)

# Y-axis: actual numbers
from matplotlib.ticker import LogLocator, FuncFormatter

def fmt_num(val, _):
    if val >= 1e12:
        return f'{val/1e12:.0f}T'
    if val >= 1e9:
        v = val / 1e9
        return f'{v:.0f}B' if v >= 10 else f'{v:.1f}B'
    if val >= 1e6:
        v = val / 1e6
        return f'{v:.0f}M' if v >= 10 else f'{v:.1f}M'
    if val >= 1e3:
        return f'{val/1e3:.0f}K'
    return f'{int(val)}'

ax.yaxis.set_major_formatter(FuncFormatter(fmt_num))
ax.yaxis.set_minor_locator(LogLocator(subs='all'))
ax.set_ylabel('Number of variant sites', fontsize=12)
ax.set_xlabel('Variant class', fontsize=12)
ax.set_title('AoU variant site counts\n(All categories, strip = TR proportion)',
             fontsize=13, fontweight='bold')

ax.set_ylim(bottom=1e3)
ax.grid(axis='y', which='major', linestyle='--', alpha=0.4, zorder=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
legend_handles = [
    mpatches.Patch(color='#7B2D8B', label='SNV'),
    mpatches.Patch(color='#FF69B4', label='INS (pink shades)'),
    mpatches.Patch(color='#CC0000', label='DEL (red shades)'),
    mpatches.Patch(facecolor='grey', hatch='///', alpha=0.85,
                   edgecolor='white', label='TR proportion (strip)'),
]
ax.legend(handles=legend_handles, loc='upper right', fontsize=9,
          framealpha=0.85)

plt.tight_layout()

OUT = ("/Users/xzhao/zxf0419@gmail.com - Google Drive/My Drive/"
       "Talkowski_Lab/gnomAD_LR/annotation/aou_I.counts_sites.barplot.pdf")
plt.savefig(OUT, dpi=200, bbox_inches='tight')
print(f"Saved: {OUT}")

OUT_PNG = OUT.replace('.pdf', '.png')
plt.savefig(OUT_PNG, dpi=150, bbox_inches='tight')
print(f"Saved: {OUT_PNG}")
