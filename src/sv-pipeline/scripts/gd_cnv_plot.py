#!/bin/python

"""
GD CNV Plotting

Generate visualisation plots for GD copy-number variant calls produced by
gd_cnv_call.py (or the combined plot_gd_cnv_output.py).  Shows depth
profiles at GD loci with annotations for segmental duplications, gene
transcripts, breakpoint intervals, and copy-number calls.

Flanking regions are compressed in display space by ``--flank-scale``
(default 0.10) so the locus body dominates the figure width.

Usage:
    python gd_cnv_plot.py \
        --calls gd_cnv_calls.tsv.gz \
        --cn-posteriors cn_posteriors.tsv.gz \
        --gd-table gd_table.tsv \
        --gtf genes.gtf.gz \
        --segdup-bed segdups.bed.gz \
        --output-dir plots/
"""

import argparse
import gzip
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Import GDLocus and GDTable from gd_cnv_pyro.py to avoid duplication
from gd_cnv_pyro import GDLocus, GDTable


# =============================================================================
# Logging
# =============================================================================


class TeeStream:
    """Write to both the original stream and a log file."""

    def __init__(self, original_stream, log_file):
        self.original_stream = original_stream
        self.log_file = log_file

    def write(self, message: str):
        self.original_stream.write(message)
        self.log_file.write(message)

    def flush(self):
        self.original_stream.flush()
        self.log_file.flush()

    def __getattr__(self, name):
        return getattr(self.original_stream, name)


def setup_logging(output_dir: str, filename: str = "plot_log.txt"):
    """Redirect stdout and stderr to both the console and a log file."""
    log_path = os.path.join(output_dir, filename)
    log_fh = open(log_path, "w")
    sys.stdout = TeeStream(sys.__stdout__, log_fh)
    sys.stderr = TeeStream(sys.__stderr__, log_fh)
    return log_fh


# =============================================================================
# Annotation Loaders
# =============================================================================


class GTFParser:
    """Parser for GTF gene annotation files."""

    def __init__(self, filepath: str):
        """
        Load gene annotations from GTF file.

        Args:
            filepath: Path to GTF file (can be gzipped)
        """
        self.genes = []
        self.transcripts = []
        self._parse_gtf(filepath)
        self._build_index()
        print(f"Loaded {len(self.genes)} genes and {len(self.transcripts)} transcripts")

    def _parse_gtf(self, filepath: str):
        """Parse GTF file and extract gene/transcript information."""
        opener = gzip.open if filepath.endswith(".gz") else open
        mode = "rt" if filepath.endswith(".gz") else "r"

        with opener(filepath, mode) as f:
            for line in f:
                if line.startswith("#"):
                    continue

                fields = line.strip().split("\t")
                if len(fields) < 9:
                    continue

                chrom, source, feature, start, end, score, strand, frame, attributes = fields

                if feature not in ["gene", "transcript"]:
                    continue

                # Parse attributes
                attr_dict = self._parse_attributes(attributes)

                record = {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "strand": strand,
                    "gene_id": attr_dict.get("gene_id", ""),
                    "gene_name": attr_dict.get("gene_name", ""),
                    "gene_type": attr_dict.get("gene_type", ""),
                }

                if feature == "gene":
                    self.genes.append(record)
                elif feature == "transcript":
                    record["transcript_id"] = attr_dict.get("transcript_id", "")
                    record["transcript_type"] = attr_dict.get("transcript_type", "")
                    self.transcripts.append(record)

    def _parse_attributes(self, attr_string: str) -> dict:
        """Parse GTF attribute string into dictionary."""
        attrs = {}
        pattern = r'(\w+)\s+"([^"]+)"'
        for match in re.finditer(pattern, attr_string):
            key, value = match.groups()
            attrs[key] = value
        return attrs

    def _build_index(self):
        """Build per-chromosome sorted arrays for fast region queries."""
        from collections import defaultdict

        genes_by_chrom = defaultdict(list)
        for g in self.genes:
            genes_by_chrom[g["chrom"]].append(g)
        self._gene_index = {}
        for chrom, genes in genes_by_chrom.items():
            genes.sort(key=lambda x: x["start"])
            self._gene_index[chrom] = (
                np.array([g["start"] for g in genes]),
                np.array([g["end"] for g in genes]),
                genes,
            )

        tx_by_chrom = defaultdict(list)
        for t in self.transcripts:
            tx_by_chrom[t["chrom"]].append(t)
        self._tx_index = {}
        for chrom, txs in tx_by_chrom.items():
            txs.sort(key=lambda x: x["start"])
            self._tx_index[chrom] = (
                np.array([t["start"] for t in txs]),
                np.array([t["end"] for t in txs]),
                txs,
            )

    def get_genes_in_region(self, chrom: str, start: int, end: int,
                            gene_types: Optional[List[str]] = None) -> List[dict]:
        """Get genes overlapping a genomic region."""
        if chrom not in self._gene_index:
            return []
        starts, ends, records = self._gene_index[chrom]
        mask = (starts < end) & (ends > start)
        indices = np.where(mask)[0]
        if gene_types:
            gene_types_set = set(gene_types)
            return [records[i] for i in indices if records[i]["gene_type"] in gene_types_set]
        return [records[i] for i in indices]

    def get_transcripts_in_region(self, chrom: str, start: int, end: int,
                                  transcript_types: Optional[List[str]] = None) -> List[dict]:
        """Get transcripts overlapping a genomic region."""
        if chrom not in self._tx_index:
            return []
        starts, ends, records = self._tx_index[chrom]
        mask = (starts < end) & (ends > start)
        indices = np.where(mask)[0]
        if transcript_types:
            tx_types_set = set(transcript_types)
            return [records[i] for i in indices if records[i]["transcript_type"] in tx_types_set]
        return [records[i] for i in indices]


class SegDupAnnotation:
    """Handler for segmental duplication annotations."""

    def __init__(self, filepath: str):
        """Load segmental duplication regions from BED file."""
        self.df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=["chr", "start", "end", "name", "score", "strand"],
            compression="gzip" if filepath.endswith(".gz") else None,
        )
        self._build_index()
        print(f"Loaded {len(self.df)} segmental duplication regions")

    def _build_index(self):
        """Build index by chromosome."""
        self.by_chrom = {}
        for chrom, group in self.df.groupby("chr"):
            self.by_chrom[chrom] = group[["start", "end"]].values

    def get_regions_in_range(self, chrom: str, start: int, end: int) -> List[Tuple[int, int]]:
        """Get segdup regions overlapping a genomic range."""
        if chrom not in self.by_chrom:
            return []
        regions = self.by_chrom[chrom]
        overlaps = (regions[:, 0] < end) & (regions[:, 1] > start)
        return [(int(r[0]), int(r[1])) for r in regions[overlaps]]


class GapsAnnotation:
    """Handler for reference assembly gap annotations."""

    def __init__(self, filepath: str):
        """Load gap regions from a BED file (optionally gzipped)."""
        self.df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            usecols=[0, 1, 2],
            names=["chr", "start", "end"],
            comment="#",
            compression="gzip" if filepath.endswith(".gz") else None,
        )
        self._build_index()
        print(f"Loaded {len(self.df)} reference gap regions")

    def _build_index(self):
        """Build index by chromosome."""
        self.by_chrom = {}
        for chrom, group in self.df.groupby("chr"):
            self.by_chrom[chrom] = group[["start", "end"]].values

    def get_regions_in_range(self, chrom: str, start: int, end: int) -> List[Tuple[int, int]]:
        """Get gap regions overlapping a genomic range."""
        if chrom not in self.by_chrom:
            return []
        regions = self.by_chrom[chrom]
        overlaps = (regions[:, 0] < end) & (regions[:, 1] > start)
        return [(int(r[0]), int(r[1])) for r in regions[overlaps]]


# =============================================================================
# Coordinate Transform
# =============================================================================


class FlankCompressor:
    """Piecewise-linear coordinate transform that compresses flanking regions.

    Maps genomic positions to *display* positions where each flanking
    region outside ``[locus_start, locus_end]`` occupies exactly
    *flank_scale* of the total display width.  The locus body fills
    the remainder.  The display range is ``[0, d_end]``.

    When *flank_scale* is ``None`` no compression occurs: each section's
    display width is proportional to its genomic width (but coordinates
    are still shifted to start at 0; use :meth:`format_genomic_ticks`
    to label axes with real genomic positions).
    """

    def __init__(
        self,
        region_start: float,
        region_end: float,
        locus_start: float,
        locus_end: float,
        flank_scale: Optional[float] = 0.20,
    ):
        self.region_start = region_start
        self.region_end = region_end
        self.locus_start = locus_start
        self.locus_end = locus_end
        self.flank_scale = flank_scale

        left_width = max(locus_start - region_start, 0)
        body_width = max(locus_end - locus_start, 0)
        right_width = max(region_end - locus_end, 0)
        total_genomic = left_width + body_width + right_width

        has_left = left_width > 0
        has_right = right_width > 0

        # D: total display width (arbitrary; keep in genomic scale for
        # comfortable numbers).
        D = max(total_genomic, 1)

        if flank_scale is None:
            # No compression – display proportional to genomic.
            self.d_left = left_width
            self.d_body = body_width
            self.d_right = right_width
        else:
            self.d_left = flank_scale * D if has_left else 0.0
            self.d_right = flank_scale * D if has_right else 0.0
            self.d_body = D - self.d_left - self.d_right
            # Guard: body must be positive.
            if self.d_body <= 0:
                n_flanks = int(has_left) + int(has_right)
                self.d_body = 0.01 * D
                per_flank = (D - self.d_body) / max(n_flanks, 1)
                self.d_left = per_flank if has_left else 0.0
                self.d_right = per_flank if has_right else 0.0

        self.d_body_start = self.d_left
        self.d_body_end = self.d_left + self.d_body
        self.d_end = self.d_body_end + self.d_right

        # Per-section scale factors (display units per genomic unit).
        self._left_scale = self.d_left / left_width if left_width > 0 else 0.0
        self._body_scale = self.d_body / body_width if body_width > 0 else 0.0
        self._right_scale = self.d_right / right_width if right_width > 0 else 0.0

    # -- forward: genomic → display ------------------------------------------

    def __call__(self, x):
        """Transform genomic coordinate(s) to display.  NaN-safe."""
        scalar = np.ndim(x) == 0
        x = np.atleast_1d(np.asarray(x, dtype=float))
        out = np.empty_like(x)
        nan = np.isnan(x)
        left = (~nan) & (x < self.locus_start)
        body = (~nan) & (x >= self.locus_start) & (x <= self.locus_end)
        right = (~nan) & (x > self.locus_end)
        out[nan] = np.nan
        out[left] = (x[left] - self.region_start) * self._left_scale
        out[body] = self.d_body_start + (x[body] - self.locus_start) * self._body_scale
        out[right] = self.d_body_end + (x[right] - self.locus_end) * self._right_scale
        return float(out[0]) if scalar else out

    # -- inverse: display → genomic ------------------------------------------

    def inverse(self, d):
        """Transform display coordinate(s) back to genomic."""
        scalar = np.ndim(d) == 0
        d = np.atleast_1d(np.asarray(d, dtype=float))
        out = np.empty_like(d)
        nan = np.isnan(d)
        left = (~nan) & (d < self.d_body_start)
        body = (~nan) & (d >= self.d_body_start) & (d <= self.d_body_end)
        right = (~nan) & (d > self.d_body_end)
        out[nan] = np.nan
        if self._left_scale > 0:
            out[left] = self.region_start + d[left] / self._left_scale
        else:
            out[left] = self.region_start
        if self._body_scale > 0:
            out[body] = self.locus_start + (d[body] - self.d_body_start) / self._body_scale
        else:
            out[body] = self.locus_start
        if self._right_scale > 0:
            out[right] = self.locus_end + (d[right] - self.d_body_end) / self._right_scale
        else:
            out[right] = self.locus_end
        return float(out[0]) if scalar else out

    # -- tick formatting helper ----------------------------------------------

    def format_genomic_ticks(self, ax, breakpoints=None, label_x=True):
        """Set up major and minor x-axis ticks with genomic-position labels.

        Major ticks (labeled) at the left / right plot edges and at each
        breakpoint boundary.  Minor ticks (unlabeled) every 10 % of the
        GD body width, extending across the flanks so the compression
        is visible.

        Args:
            ax: Matplotlib axis.
            breakpoints: List of ``(start, end)`` genomic coordinate
                pairs for breakpoint intervals.  When *None*, only edge
                ticks are drawn.
        """
        inv = self.inverse

        # -- Major ticks: plot edges + breakpoint boundaries ---------------
        major_display = [0.0, self.d_end]
        if breakpoints:
            for bp_start, bp_end in breakpoints:
                major_display.append(self(bp_start))
                major_display.append(self(bp_end))
        # De-duplicate (round to avoid float noise) and sort
        major_display = sorted({round(d, 6) for d in major_display})

        ax.set_xticks(major_display)
        if label_x:
            major_labels = [f"{inv(d):,.0f}" for d in major_display]
            ax.set_xticklabels(major_labels, fontsize=7, rotation=45, ha="right", color='red')
        else:
            ax.set_xticklabels([])
        ax.tick_params(axis='x', which='major', length=5, width=1, color='red', labelcolor='red')

        # -- Minor ticks: every 10 % of body width (genomic) across full
        #    display range so the compression is visible ----------------
        body_genomic = self.locus_end - self.locus_start
        if body_genomic > 0:
            step = body_genomic * 0.10
            minor_genomic = []
            # Left from locus_start into left flank
            pos = self.locus_start - step
            while pos >= self.region_start:
                minor_genomic.append(pos)
                pos -= step
            # Right from locus_start through body and into right flank
            pos = self.locus_start + step
            while pos <= self.region_end:
                minor_genomic.append(pos)
                pos += step
            minor_display = sorted(
                float(d) for d in self(np.array(minor_genomic)))
            # Drop any that land on (or very near) a major tick
            major_set = {round(d, 2) for d in major_display}
            minor_display = [d for d in minor_display
                            if round(d, 2) not in major_set]
            ax.set_xticks(minor_display, minor=True)
            ax.tick_params(axis='x', which='minor', length=3, width=0.5)


# =============================================================================
# Plotting Helpers
# =============================================================================


def get_sample_columns(df: pd.DataFrame) -> list:
    """Extract sample column names from DataFrame."""
    metadata_cols = ["Cluster", "Chr", "Start", "End", "source_file", "Bin"]
    return [col for col in df.columns if col not in metadata_cols]


def draw_annotations_panel(
    ax,
    locus: GDLocus,
    region_start: int,
    region_end: int,
    chrom: str,
    title: str,
    gtf: Optional[GTFParser] = None,
    segdup: Optional[SegDupAnnotation] = None,
    gaps: Optional[GapsAnnotation] = None,
    show_gd_entries: bool = True,
    min_gene_label_spacing: float = 0.05,
    xform: Optional[FlankCompressor] = None,
):
    """
    Draw annotations panel with genes, segdups, reference gaps, and breakpoints.

    All x-coordinates are mapped through *xform* so that flanking regions
    can be visually compressed.

    Args:
        ax: Matplotlib axis to draw on
        locus: GDLocus object
        region_start: Region start position (genomic)
        region_end: Region end position (genomic)
        chrom: Chromosome name
        title: Panel title
        gtf: Optional GTFParser for gene annotations
        segdup: Optional SegDupAnnotation
        gaps: Optional GapsAnnotation for reference assembly gaps
        show_gd_entries: Whether to show GD entry regions and breakpoint intervals
        min_gene_label_spacing: Minimum distance between gene label centres as a
            fraction of the plot width.
        xform: Coordinate transform (genomic → display).  When *None* a
            default identity-like transform (no compression) is created.
    """
    if xform is None:
        xform = FlankCompressor(region_start, region_end,
                                locus.start, locus.end, flank_scale=None)

    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0, 0.25)
    ax.set_ylabel("Annotations")
    ax.set_title(title)

    # Shade flanking regions
    flank_defs = []
    if region_start < locus.start:
        flank_defs.append((region_start, locus.start, "left flank"))
    if locus.end < region_end:
        flank_defs.append((locus.end, region_end, "right flank"))
    for flank_start, flank_end, flank_name in flank_defs:
        d_start = xform(flank_start)
        d_end = xform(flank_end)
        ax.axvspan(d_start, d_end, alpha=0.08, color="black", zorder=0)
        ax.text((d_start + d_end) / 2, 0.225,
                flank_name, ha="center", va="center",
                fontsize=7, color="gray", style="italic")

    # Draw breakpoint intervals
    if show_gd_entries:
        intervals = locus.get_intervals()
        colors = plt.cm.Set3(np.linspace(0, 1, len(intervals)))
        for i, (start, end, name) in enumerate(intervals):
            ds, de = xform(start), xform(end)
            ax.axvspan(ds, de, alpha=0.15, color=colors[i], label=name, zorder=0)
            ax.text((ds + de) / 2, 0.225, name, ha="center", va="center", fontsize=7)

    # Draw breakpoint ranges as shaded regions
    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ds, de = xform(bp_start), xform(bp_end)
        ax.axvspan(ds, de, ymin=0.04, ymax=0.15, alpha=0.3, color="red", zorder=1)
        bp_name = locus.breakpoint_names[i] if i < len(locus.breakpoint_names) else str(i + 1)
        ax.text((ds + de) / 2, 0.02, f"BP {bp_name}", ha="center", va="center",
                fontsize=7, fontweight="bold", color="darkred", zorder=10)

    # Draw segdup regions
    if segdup:
        sd_regions = segdup.get_regions_in_range(chrom, region_start, region_end)
        for sd_start, sd_end in sd_regions:
            ax.axvspan(xform(max(sd_start, region_start)),
                       xform(min(sd_end, region_end)),
                       ymin=0.2, ymax=0.5, alpha=0.2, color="orange", zorder=2)

    # Draw reference gap regions
    if gaps:
        gap_regions = gaps.get_regions_in_range(chrom, region_start, region_end)
        for gap_start, gap_end in gap_regions:
            ax.axvspan(
                xform(max(gap_start, region_start)),
                xform(min(gap_end, region_end)),
                alpha=0.25, color="lightgray", hatch="//", edgecolor="gray",
                linewidth=0.5, zorder=2,
            )

    # Draw genes
    if gtf:
        genes = gtf.get_genes_in_region(chrom, region_start, region_end,
                                        gene_types=["protein_coding"])
        y_pos = 0.2
        display_width = xform.d_end
        min_spacing = min_gene_label_spacing * display_width
        last_labeled_center = -np.inf
        for gene in genes:
            gene_start = max(gene["start"], region_start)
            gene_end = min(gene["end"], region_end)
            gene_center = (gene_start + gene_end) / 2
            d_gs, d_ge = xform(gene_start), xform(gene_end)
            d_gc = xform(gene_center)
            ax.hlines(y_pos, d_gs, d_ge, colors="blue", linewidth=4, zorder=3)
            if d_gc - last_labeled_center >= min_spacing:
                ax.text(d_gc, y_pos - 0.05,
                        gene["gene_name"], ha="center", va="bottom", fontsize=7, zorder=3)
                last_labeled_center = d_gc

    ax.set_yticks([])
    xform.format_genomic_ticks(ax, locus.breakpoints, label_x=False)


def _build_ploidy_lookup(
    ploidy_df: Optional[pd.DataFrame],
) -> Dict[Tuple[str, str], int]:
    """Build (sample, contig) -> ploidy mapping from a ploidy DataFrame."""
    lookup: Dict[Tuple[str, str], int] = {}
    if ploidy_df is not None:
        for _, row in ploidy_df.iterrows():
            lookup[(str(row["sample"]), str(row["contig"]))] = int(row["ploidy"])
    return lookup


def _sort_samples_by_ploidy(
    sample_cols: List[str],
    chrom: str,
    ploidy_lookup: Dict[Tuple[str, str], int],
) -> List[str]:
    """Return *sample_cols* sorted by ploidy (ascending), then by name."""
    return sorted(sample_cols, key=lambda s: (ploidy_lookup.get((s, chrom), 2), s))


def _build_gap_positions(
    bin_mids: np.ndarray,
    bin_starts: np.ndarray,
    bin_ends: np.ndarray,
) -> np.ndarray:
    """Build a positions array with NaN sentinels at bin gaps."""
    n = len(bin_starts)
    if n == 0:
        return np.array([], dtype=float)
    if n == 1:
        return bin_mids.copy()
    bin_spacing = np.median(bin_starts[1:] - bin_ends[:-1])
    gap_threshold = 3 * bin_spacing
    is_gap = (bin_starts[1:] - bin_ends[:-1]) > gap_threshold
    n_gaps = int(is_gap.sum())
    result = np.full(n + n_gaps, np.nan)
    cum_gaps = np.empty(n, dtype=int)
    cum_gaps[0] = 0
    cum_gaps[1:] = np.cumsum(is_gap)
    result[np.arange(n) + cum_gaps] = bin_mids
    return result


def _depths_with_gaps(
    raw_depths: np.ndarray,
    plot_positions: np.ndarray,
) -> np.ndarray:
    """Map per-bin depths onto the gap-aware position array."""
    result = np.full_like(plot_positions, np.nan)
    non_nan = ~np.isnan(plot_positions)
    result[non_nan] = raw_depths
    return result


def _build_heatmap_matrix(
    region_df: pd.DataFrame,
    sample_cols: List[str],
    region_start: int,
    region_end: int,
    xform: Optional[FlankCompressor] = None,
    n_viz_bins: int = 1000,
) -> np.ndarray:
    """Build a dense (n_samples, n_viz_bins) matrix for heatmap display.

    When *xform* is provided, pixel columns are spaced evenly in display
    coordinates and then inverse-transformed to genomic space for the
    searchsorted lookup so that compressed flanks occupy fewer pixels.
    """
    n_samples = len(sample_cols)
    viz_matrix = np.full((n_samples, n_viz_bins), np.nan)

    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    data = region_df[sample_cols].values  # shape (n_bins, n_samples)

    if len(bin_starts) == 0:
        return viz_matrix

    if xform is not None:
        display_span = max(xform.d_end, 1)
        pixel_display = (np.arange(n_viz_bins) + 0.5) * (display_span / n_viz_bins)
        pixel_genomic = xform.inverse(pixel_display)
    else:
        region_span = max(region_end - region_start, 1)
        pixel_genomic = (
            (np.arange(n_viz_bins) + 0.5) * (region_span / n_viz_bins)
            + region_start
        )

    bin_idx = np.searchsorted(bin_starts, pixel_genomic, side='right') - 1
    safe_idx = np.clip(bin_idx, 0, len(bin_starts) - 1)
    valid = (bin_idx >= 0) & (pixel_genomic < bin_ends[safe_idx])
    viz_matrix[:, valid] = data[bin_idx[valid], :].T

    return viz_matrix


# =============================================================================
# Overview Column Drawing
# =============================================================================


def _draw_overview_column(
    axes_col: List,
    region_df: pd.DataFrame,
    locus: GDLocus,
    calls_df: pd.DataFrame,
    region_start: int,
    region_end: int,
    carrier_cols: List[str],
    non_carrier_cols: List[str],
    all_ploidies: List[int],
    ploidy_lookup: Dict[Tuple[str, str], int],
    sample_cols: List[str],
    carriers,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    col_title: str,
    min_gene_label_spacing: float = 0.05,
    show_colorbar: bool = True,
    gaps: Optional[GapsAnnotation] = None,
    xform: Optional[FlankCompressor] = None,
):
    """Draw annotation + carrier heatmap + non-carrier heatmap + mean-depth
    panels into a list of axes (one column of the overview figure)."""
    if xform is None:
        xform = FlankCompressor(region_start, region_end,
                                locus.start, locus.end, flank_scale=None)

    chrom = locus.chrom
    n_carriers = len(carrier_cols)
    n_non_carriers = len(non_carrier_cols)
    n_ploidy_panels = len(all_ploidies)

    d_start = 0.0
    d_end = xform.d_end

    # Panel 1: Gene annotations and segdup regions
    ax = axes_col[0]
    draw_annotations_panel(
        ax, locus, region_start, region_end, chrom, col_title,
        gtf=gtf, segdup=segdup, gaps=gaps, show_gd_entries=True,
        min_gene_label_spacing=min_gene_label_spacing, xform=xform,
    )

    # Panel 2: Carriers heatmap
    ax = axes_col[1]
    if n_carriers > 0:
        viz_matrix = _build_heatmap_matrix(
            region_df, carrier_cols, region_start, region_end, xform=xform)
        cmap = plt.cm.RdBu_r.copy()
        cmap.set_bad(color='white')
        ax.imshow(viz_matrix, aspect="auto", cmap=cmap,
                  vmin=0, vmax=4, interpolation="nearest",
                  extent=[d_start, d_end, n_carriers, 0])
        ax.set_ylabel(f"Carriers (n={n_carriers})")
        ax.set_yticks([])
        ax.set_xlim(d_start, d_end)
        ax.set_ylim(0, n_carriers)
        xform.format_genomic_ticks(ax, locus.breakpoints, label_x=False)
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'No carriers', ha='center', va='center',
                transform=ax.transAxes)

    # Panel 3: Non-carriers heatmap
    ax = axes_col[2]
    if n_non_carriers > 0:
        viz_matrix = _build_heatmap_matrix(
            region_df, non_carrier_cols, region_start, region_end, xform=xform)
        cmap = plt.cm.RdBu_r.copy()
        cmap.set_bad(color='white')
        im = ax.imshow(viz_matrix, aspect="auto", cmap=cmap,
                       vmin=0, vmax=4, interpolation="nearest",
                       extent=[d_start, d_end, n_non_carriers, 0])
        ax.set_ylabel(f"Non-carriers (n={n_non_carriers})")
        ax.set_yticks([])
        ax.set_xlim(d_start, d_end)
        ax.set_ylim(0, n_non_carriers)
        xform.format_genomic_ticks(ax, locus.breakpoints, label_x=False)
        if show_colorbar:
            cbar = plt.colorbar(im, ax=ax, orientation='horizontal',
                                pad=0.08, aspect=40)
            cbar.set_label("Normalized read depth")
    else:
        ax.axis('off')
        ax.text(0.5, 0.5, 'No non-carriers', ha='center', va='center',
                transform=ax.transAxes)

    # ---- Mean depth panels: one per ploidy state ----
    bin_mids = (region_df["Start"].values + region_df["End"].values) / 2
    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    # Build gap-aware positions in genomic space, then transform to display
    plot_positions_genomic = _build_gap_positions(bin_mids, bin_starts, bin_ends)
    display_positions = xform(plot_positions_genomic)

    # Pre-transform breakpoint coords
    bp_display = [(xform(s), xform(e)) for s, e in locus.breakpoints]

    for panel_idx, ploidy_val in enumerate(all_ploidies):
        ax = axes_col[3 + panel_idx]

        ploidy_samples = [s for s in sample_cols
                          if ploidy_lookup.get((s, chrom), 2) == ploidy_val]
        ploidy_carriers = [s for s in ploidy_samples if s in carriers]
        ploidy_non_carriers = [s for s in ploidy_samples if s not in carriers]

        if len(ploidy_carriers) > 0 and len(ploidy_non_carriers) > 0:
            nc_depths = _depths_with_gaps(
                region_df[ploidy_non_carriers].mean(axis=1).values, plot_positions_genomic)
            ax.plot(display_positions, nc_depths, "b-", linewidth=1, alpha=0.7,
                    label=f"Non-carriers (n={len(ploidy_non_carriers)})")

            carrier_calls = calls_df[
                (calls_df["cluster"] == locus.cluster) &
                (calls_df["is_carrier"]) &
                (calls_df["is_best_match"]) &
                (calls_df["sample"].isin(ploidy_carriers))
            ].copy()

            def get_bp_labels(row):
                for entry in locus.gd_entries:
                    if (entry["start_GRCh38"] == row["start"] and
                            entry["end_GRCh38"] == row["end"] and
                            entry["svtype"] == row["svtype"]):
                        return f"{entry['BP1']}-{entry['BP2']}"
                for entry in locus.gd_entries:
                    if entry["start_GRCh38"] == row["start"] and entry["end_GRCh38"] == row["end"]:
                        return f"{entry['BP1']}-{entry['BP2']}"
                return "unknown"

            if len(carrier_calls) > 0:
                carrier_calls["bp_interval"] = carrier_calls.apply(
                    get_bp_labels, axis=1)
                del_colors = ['#FF6B6B', '#FF4757', '#EE5A6F', '#C23B4E']
                dup_colors = ['#6B9BD1', '#4169E1', '#5080D0', '#3A5BB8']
                color_idx = 0
                plotted_any = False
                for (svtype, bp_interval), group in carrier_calls.groupby(
                        ["svtype", "bp_interval"]):
                    valid = [s for s in group["sample"].unique()
                             if s in ploidy_carriers]
                    if len(valid) > 0:
                        gd = _depths_with_gaps(
                            region_df[valid].mean(axis=1).values,
                            plot_positions_genomic)
                        if not all(np.isnan(gd)):
                            color = (del_colors if svtype == "DEL"
                                     else dup_colors)[color_idx % 4]
                            ax.plot(display_positions, gd, "-", linewidth=1.5,
                                    alpha=0.8, color=color,
                                    label=f"{svtype} {bp_interval} (n={len(valid)})")
                            color_idx += 1
                            plotted_any = True
                if plotted_any:
                    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.35),
                              ncol=6, fontsize=7, frameon=True, fancybox=True)
        else:
            if len(ploidy_samples) > 0:
                md = _depths_with_gaps(
                    region_df[ploidy_samples].mean(axis=1).values,
                    plot_positions_genomic)
                ax.plot(display_positions, md, "b-", linewidth=1.5,
                        label=f"All (n={len(ploidy_samples)})")
                ax.legend(fontsize=7)

        for d_bp_start, d_bp_end in bp_display:
            ax.axvline(d_bp_start, color="red", linestyle="-", alpha=1, zorder=5)
            ax.axvline(d_bp_end, color="red", linestyle="-", alpha=1, zorder=5)
            ax.axvspan(d_bp_start, d_bp_end, alpha=0.1, color="red", zorder=0)

        ax.axhline(float(ploidy_val), color="gray", linestyle=":", alpha=0.5)
        ax.set_xlim(d_start, d_end)
        ax.set_ylim(0, max(4, ploidy_val + 2))
        ax.set_ylabel(f"Depth (P={ploidy_val})")
        ax.grid(True, alpha=0.3, axis="y")

        if panel_idx == n_ploidy_panels - 1:
            ax.set_xlabel(f"Position on {chrom}")
            xform.format_genomic_ticks(ax, locus.breakpoints)
        else:
            ax.set_xticks([])


# =============================================================================
# Top-level Plot Functions
# =============================================================================


def plot_locus_overview(
    locus: GDLocus,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
    ploidy_df: Optional[pd.DataFrame] = None,
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
    raw_sample_medians: Optional[Dict[str, float]] = None,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
):
    """
    Create overview plot for a GD locus showing all samples.

    When *raw_counts_df* is provided the figure has two columns:
      - **Left**: normalised raw depth (before filtering/rebinning)
      - **Right**: processed rebinned counts from the model
    """
    if not locus.breakpoints:
        print(f"  Warning: No breakpoints defined for locus {locus.cluster}, skipping")
        return

    chrom = locus.chrom
    mask = (
        (depth_df["Cluster"] == locus.cluster) &
        (depth_df["Chr"] == chrom)
    )
    region_df = depth_df[mask].copy().sort_values("Start")

    if len(region_df) == 0:
        print(f"  Warning: No depth data found for locus {locus.cluster} "
              f"({chrom}:{locus.start:,}-{locus.end:,}), skipping plot")
        return

    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())

    # Build coordinate transform
    xform = FlankCompressor(region_start, region_end,
                            locus.start, locus.end, flank_scale=flank_scale)

    sample_cols = get_sample_columns(region_df)
    ploidy_lookup = _build_ploidy_lookup(ploidy_df)

    carriers = calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["is_carrier"])
    ]["sample"].unique()

    carrier_cols = _sort_samples_by_ploidy(
        [s for s in sample_cols if s in carriers], chrom, ploidy_lookup)
    non_carrier_cols = _sort_samples_by_ploidy(
        [s for s in sample_cols if s not in carriers], chrom, ploidy_lookup)
    n_carriers = len(carrier_cols)
    n_non_carriers = len(non_carrier_cols)

    all_ploidies = sorted(set(
        ploidy_lookup.get((s, chrom), 2) for s in sample_cols
    ))
    n_ploidy_panels = len(all_ploidies)

    # ---- Build raw-normalised DataFrame for the left column (if available) ----
    have_raw = (raw_counts_df is not None and raw_sample_medians is not None)
    raw_region_df: Optional[pd.DataFrame] = None
    if have_raw:
        raw_mask = (
            (raw_counts_df["Chr"] == chrom)
            & (raw_counts_df["End"] > region_start)
            & (raw_counts_df["Start"] < region_end)
        )
        raw_region = raw_counts_df[raw_mask].copy().sort_values("Start")
        if len(raw_region) > 0:
            raw_norm_data = {}
            usable_samples = []
            for s in sample_cols:
                if s in raw_region.columns and s in raw_sample_medians:
                    med = raw_sample_medians[s]
                    if med > 0:
                        raw_norm_data[s] = 2.0 * raw_region[s].values.astype(float) / med
                        usable_samples.append(s)
            if usable_samples:
                raw_region_df = pd.concat(
                    [raw_region[["Chr", "Start", "End"]].reset_index(drop=True),
                     pd.DataFrame(raw_norm_data, columns=usable_samples)],
                    axis=1,
                )
        if raw_region_df is None or len(raw_region_df) == 0:
            have_raw = False
            raw_region_df = None

    # ---- figure layout ----
    n_cols = 2 if have_raw else 1
    carrier_height = min(6, max(1, n_carriers * 0.15)) if n_carriers > 0 else 0.5
    non_carrier_height = min(6, max(1, n_non_carriers * 0.15)) if n_non_carriers > 0 else 0.5
    depth_panel_height = 1.0
    fig_height = 3 + carrier_height + non_carrier_height + n_ploidy_panels * depth_panel_height

    height_ratios = (
        [1, carrier_height, non_carrier_height]
        + [depth_panel_height] * n_ploidy_panels
    )
    n_rows = 3 + n_ploidy_panels
    fig_width = 16 * n_cols

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(fig_width, fig_height),
        gridspec_kw={"height_ratios": height_ratios},
        squeeze=False,
    )

    locus_label = f"{locus.cluster} ({chrom}:{locus.start:,}-{locus.end:,})"

    if have_raw:
        _draw_overview_column(
            [axes[r, 0] for r in range(n_rows)],
            raw_region_df, locus, calls_df,
            region_start, region_end,
            [s for s in carrier_cols if s in raw_region_df.columns],
            [s for s in non_carrier_cols if s in raw_region_df.columns],
            all_ploidies, ploidy_lookup,
            [s for s in sample_cols if s in raw_region_df.columns],
            carriers, gtf, segdup,
            col_title=f"Raw normalised — {locus_label}",
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=False,
            gaps=gaps,
            xform=xform,
        )
        _draw_overview_column(
            [axes[r, 1] for r in range(n_rows)],
            region_df, locus, calls_df,
            region_start, region_end,
            carrier_cols, non_carrier_cols,
            all_ploidies, ploidy_lookup,
            sample_cols, carriers, gtf, segdup,
            col_title=f"Processed — {locus_label}",
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=False,
            gaps=gaps,
            xform=xform,
        )
    else:
        _draw_overview_column(
            [axes[r, 0] for r in range(n_rows)],
            region_df, locus, calls_df,
            region_start, region_end,
            carrier_cols, non_carrier_cols,
            all_ploidies, ploidy_lookup,
            sample_cols, carriers, gtf, segdup,
            col_title=locus_label,
            min_gene_label_spacing=min_gene_label_spacing,
            show_colorbar=False,
            gaps=gaps,
            xform=xform,
        )

    plt.tight_layout()

    locus_dir = os.path.join(output_dir, "locus_plots")
    os.makedirs(locus_dir, exist_ok=True)
    filename = f"{locus.cluster.replace('/', '_')}_overview.png"
    plt.savefig(os.path.join(locus_dir, filename), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Created: locus_plots/{filename}")


def plot_sample_at_locus(
    sample_id: str,
    locus: GDLocus,
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
    min_gene_label_spacing: float = 0.05,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
):
    """Create detailed plot for a specific sample at a GD locus."""
    if not locus.breakpoints:
        return

    chrom = locus.chrom
    mask = (
        (depth_df["Cluster"] == locus.cluster) &
        (depth_df["Chr"] == chrom)
    )
    region_df = depth_df[mask].copy().sort_values("Start")

    if len(region_df) == 0 or sample_id not in region_df.columns:
        return

    region_start = int(region_df["Start"].min())
    region_end = int(region_df["End"].max())

    xform = FlankCompressor(region_start, region_end,
                            locus.start, locus.end, flank_scale=flank_scale)

    sample_calls = calls_df[
        (calls_df["cluster"] == locus.cluster) &
        (calls_df["sample"] == sample_id)
    ]

    fig, axes = plt.subplots(3, 1, figsize=(14, 8),
                             gridspec_kw={"height_ratios": [1, 2, 1]})

    bin_starts = region_df["Start"].values
    bin_ends = region_df["End"].values
    bin_mids = (bin_starts + bin_ends) / 2
    sample_depth = region_df[sample_id].values

    # Transform bar geometry
    d_bin_mids = xform(bin_mids)
    d_bar_widths = xform(bin_ends) - xform(bin_starts)

    # Panel 1: Annotations
    is_carrier = sample_calls["is_carrier"].any() if len(sample_calls) > 0 else False
    carrier_str = " [CARRIER]" if is_carrier else ""
    title = f"{sample_id} at {locus.cluster}{carrier_str}"
    draw_annotations_panel(axes[0], locus, region_start, region_end, chrom, title, gtf, segdup,
                           gaps=gaps, show_gd_entries=True,
                           min_gene_label_spacing=min_gene_label_spacing, xform=xform)
    axes[0].set_ylabel("Annotations")

    # Panel 2: Depth profile
    ax = axes[1]
    ax.bar(d_bin_mids, sample_depth, width=d_bar_widths * 0.9, alpha=0.7,
           color="steelblue", edgecolor="none")

    intervals = locus.get_intervals()
    for i, (start, end, name) in enumerate(intervals):
        ax.axvspan(xform(start), xform(end), alpha=0.1,
                   color=plt.cm.Set3(i / max(len(intervals), 1)))

    ax.axhline(2.0, color="green", linestyle="-", alpha=0.5, linewidth=1, label="CN=2")
    ax.axhline(1.0, color="orange", linestyle="--", alpha=0.5, linewidth=1, label="CN=1")
    ax.axhline(3.0, color="purple", linestyle="--", alpha=0.5, linewidth=1, label="CN=3")

    for i, (bp_start, bp_end) in enumerate(locus.breakpoints):
        ax.axvspan(xform(bp_start), xform(bp_end), alpha=0.2, color="red", zorder=0)

    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0, 5)
    ax.set_ylabel("Normalized Depth")
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3, axis="y")
    xform.format_genomic_ticks(ax, locus.breakpoints)

    # Panel 3: Call summary
    ax = axes[2]
    ax.set_xlim(0.0, xform.d_end)
    ax.set_ylim(0, 1)

    y_pos = 0.5
    for _, call in sample_calls.iterrows():
        if call["is_carrier"]:
            color = "red" if call["svtype"] == "DEL" else "blue"
            alpha = 0.8 if call["is_best_match"] else 0.4
            ds, de = xform(call["start"]), xform(call["end"])
            ax.axvspan(ds, de, alpha=alpha, color=color)
            label = f"{call['svtype']} CN={call['mean_cn']:.2f}"
            if call["is_best_match"]:
                label += " (best)"
            ax.text((ds + de) / 2, y_pos, label,
                    ha="center", va="center", fontsize=9,
                    fontweight="bold" if call["is_best_match"] else "normal")

    ax.set_yticks([])
    ax.set_xlabel(f"Position on {chrom}")
    ax.set_ylabel("Calls")
    xform.format_genomic_ticks(ax, locus.breakpoints)

    plt.tight_layout()

    sample_dir = os.path.join(output_dir, "sample_plots", locus.cluster.replace("/", "_"))
    os.makedirs(sample_dir, exist_ok=True)
    filename = f"{sample_id.replace('/', '_')}.png"
    plt.savefig(os.path.join(sample_dir, filename), dpi=150, bbox_inches="tight")
    plt.close()


def plot_carrier_summary(calls_df: pd.DataFrame, output_dir: str):
    """Create summary plots of carrier calls across all loci."""
    carriers = calls_df[calls_df["is_carrier"]].copy()

    if len(carriers) == 0:
        print("No carriers to plot.")
        return

    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    ax = axes[0]
    locus_counts = carriers.groupby(["cluster", "svtype"])["sample"].nunique().unstack(fill_value=0)
    locus_counts.plot(kind="barh", ax=ax, color={"DEL": "red", "DUP": "blue"}, alpha=0.7)
    ax.set_xlabel("Number of Carriers")
    ax.set_ylabel("Locus")
    ax.set_title("Carriers per GD Locus")
    ax.legend(title="SV Type")

    ax = axes[1]
    all_scores = carriers["log_prob_score"].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 21)
    for svtype, color in [("DEL", "red"), ("DUP", "blue")]:
        sv_data = carriers[carriers["svtype"] == svtype]
        if len(sv_data) > 0:
            ax.hist(sv_data["log_prob_score"], bins=bins, alpha=0.6, color=color,
                    label=f"{svtype} (n={len(sv_data)})", edgecolor="black")
    ax.set_xlabel("Log Probability Score")
    ax.set_ylabel("Count")
    ax.set_title("Log Probability Score Distribution in Carriers")
    ax.legend()

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "carrier_summary.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("Created: carrier_summary.png")


def plot_confidence_distribution(calls_df: pd.DataFrame, output_dir: str):
    """Plot distribution of confidence scores."""
    fig, axes = plt.subplots(1, 2, figsize=(8, 4))

    ax = axes[0]
    carriers = calls_df[calls_df["is_carrier"]]
    non_carriers = calls_df[~calls_df["is_carrier"]]
    all_scores = calls_df["log_prob_score"].dropna()
    bins = np.linspace(all_scores.min(), all_scores.max(), 31)
    ax.hist(non_carriers["log_prob_score"], bins=bins, alpha=0.6, label="Non-carriers",
            color="gray", edgecolor="black")
    ax.hist(carriers["log_prob_score"], bins=bins, alpha=0.6, label="Carriers",
            color="green", edgecolor="black")
    ax.set_xlabel("Log Probability Score")
    ax.set_ylabel("Count")
    ax.set_yscale("log")
    ax.legend()

    ax = axes[1]
    non_carriers_mask = ~calls_df["is_carrier"]
    ax.scatter(calls_df.loc[non_carriers_mask, "mean_depth"],
               calls_df.loc[non_carriers_mask, "log_prob_score"],
               c="gray", alpha=0.1, s=10, label="Non-carriers")
    carriers_mask = calls_df["is_carrier"]
    ax.scatter(calls_df.loc[carriers_mask, "mean_depth"],
               calls_df.loc[carriers_mask, "log_prob_score"],
               c="green", alpha=0.5, s=10, label="Carriers")
    ax.set_xlabel("Mean Depth (normalized)")
    ax.set_ylabel("Confidence")

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "confidence_distribution.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("Created: confidence_distribution.png")


def create_carrier_pdf(
    calls_df: pd.DataFrame,
    depth_df: pd.DataFrame,
    gd_table: GDTable,
    gtf: Optional[GTFParser],
    segdup: Optional[SegDupAnnotation],
    output_dir: str,
    padding: int = 50000,
    min_gene_label_spacing: float = 0.05,
    raw_counts_df: Optional[pd.DataFrame] = None,
    gaps: Optional[GapsAnnotation] = None,
    flank_scale: float = 0.20,
):
    """Create a PDF with plots for all carrier samples."""
    carriers = calls_df[calls_df["is_carrier"]]

    if len(carriers) == 0:
        print("No carriers to include in PDF.")
        return

    # Precompute per-sample genome-wide autosomal median raw counts
    raw_sample_medians: Dict[str, float] = {}
    if raw_counts_df is not None:
        autosomal = [str(c) for c in range(1, 23)] + [f"chr{c}" for c in range(1, 23)]
        raw_auto = raw_counts_df[raw_counts_df["Chr"].isin(autosomal)]
        raw_sample_cols = [c for c in raw_counts_df.columns
                          if c not in ("Chr", "Start", "End", "source_file", "Bin")]
        for s in raw_sample_cols:
            med = raw_auto[s].median()
            if med > 0:
                raw_sample_medians[s] = med
        print(f"  Computed raw autosomal medians for {len(raw_sample_medians)} samples")

    pdf_path = os.path.join(output_dir, "carrier_plots.pdf")
    print(f"Creating carrier PDF: {pdf_path}")

    with PdfPages(pdf_path) as pdf:
        for cluster in sorted(carriers["cluster"].unique()):
            locus = gd_table.loci.get(cluster)
            if locus is None:
                continue

            if not locus.breakpoints:
                print(f"  Warning: No breakpoints for {cluster}, skipping")
                continue

            cluster_carriers = carriers[carriers["cluster"] == cluster]["sample"].unique()
            print(f"  Adding {len(cluster_carriers)} carriers for {cluster}")

            chrom = locus.chrom
            locus_mask = (
                (depth_df["Cluster"] == cluster) &
                (depth_df["Chr"] == chrom)
            )
            region_df = depth_df[locus_mask].sort_values("Start")
            if len(region_df) == 0:
                continue

            region_start = int(region_df["Start"].min())
            region_end = int(region_df["End"].max())

            # Build coordinate transform for this locus
            xform = FlankCompressor(region_start, region_end,
                                    locus.start, locus.end,
                                    flank_scale=flank_scale)

            cluster_calls_df = calls_df[calls_df["cluster"] == cluster]
            _raw_region = None
            if raw_counts_df is not None:
                _raw_mask = (
                    (raw_counts_df["Chr"] == chrom)
                    & (raw_counts_df["End"] > region_start)
                    & (raw_counts_df["Start"] < region_end)
                )
                _rr = raw_counts_df[_raw_mask].sort_values("Start")
                if len(_rr) > 0:
                    _raw_region = _rr

            for sample_id in sorted(cluster_carriers):
                if sample_id not in region_df.columns:
                    continue

                sample_calls = cluster_calls_df[cluster_calls_df["sample"] == sample_id]

                fig, axes = plt.subplots(2, 1, figsize=(12, 4),
                                         gridspec_kw={"height_ratios": [1, 2]})

                bin_mids = (region_df["Start"].values + region_df["End"].values) / 2
                sample_depth = region_df[sample_id].values

                # Transform bar geometry
                d_bin_mids = xform(bin_mids)
                d_bar_widths = xform(region_df["End"].values) - xform(region_df["Start"].values)

                # Panel 1: Annotations
                carrier_call = (
                    sample_calls[sample_calls["is_carrier"]].iloc[0]
                    if len(sample_calls[sample_calls["is_carrier"]]) > 0
                    else None
                )
                call_info = (
                    f" - {carrier_call['svtype']} confidence={carrier_call['log_prob_score']:.2f}"
                    if carrier_call is not None else ""
                )
                title = f"{sample_id} at {cluster}{call_info}"
                draw_annotations_panel(
                    axes[0], locus, region_start, region_end, chrom, title,
                    gtf, segdup, gaps=gaps, show_gd_entries=False,
                    min_gene_label_spacing=min_gene_label_spacing, xform=xform,
                )

                # Panel 2: Depth
                ax = axes[1]

                best_call = None
                if len(sample_calls) > 0:
                    carrier_calls_local = sample_calls[sample_calls["is_carrier"]]
                    if len(carrier_calls_local) > 0:
                        best_matches = carrier_calls_local[carrier_calls_local["is_best_match"]]
                        if len(best_matches) > 0:
                            best_call = best_matches.iloc[0]
                        else:
                            best_call = carrier_calls_local.loc[
                                carrier_calls_local["log_prob_score"].idxmax()]

                if best_call is not None:
                    interval_start = best_call["start"]
                    interval_end = best_call["end"]
                    mean_depth = best_call["mean_depth"]
                    svtype = best_call["svtype"]
                    color = "#FF6B6B" if svtype == "DEL" else "#6B9BD1"
                    d_is, d_ie = xform(interval_start), xform(interval_end)
                    ax.axvspan(d_is, d_ie, alpha=0.2, color=color, zorder=1,
                               label=f"{svtype} region")
                    ax.hlines(mean_depth, d_is, d_ie,
                              colors="black", linewidth=2.5, alpha=0.8, zorder=2,
                              label=f"Mean depth={mean_depth:.2f}")

                ax.bar(d_bin_mids, sample_depth, width=d_bar_widths * 0.9, alpha=0.6,
                       color="steelblue", edgecolor="none", zorder=3)

                # Overlay raw depth trace if available
                if _raw_region is not None and sample_id in raw_sample_medians:
                    if sample_id in _raw_region.columns:
                        raw_mids = (_raw_region["Start"].values + _raw_region["End"].values) / 2
                        raw_depth = _raw_region[sample_id].values.astype(float)
                        norm_raw = 2.0 * raw_depth / raw_sample_medians[sample_id]
                        ax.plot(xform(raw_mids), norm_raw, color="darkorange", linewidth=0.6,
                                alpha=0.7, zorder=4, label="Raw depth")

                ax.axhline(2.0, color="green", linestyle="-", alpha=0.4, linewidth=1.5,
                           label="Reference CN=2", zorder=0)
                ax.axhline(1.0, color="orange", linestyle=":", alpha=0.4, linewidth=1, zorder=0)
                ax.axhline(3.0, color="purple", linestyle=":", alpha=0.4, linewidth=1, zorder=0)

                ax.set_xlim(0.0, xform.d_end)
                ax.set_ylim(0, 5)
                ax.set_xlabel(f"Position on {chrom}")
                ax.set_ylabel("Normalized Depth")
                ax.grid(True, alpha=0.3, axis="y", zorder=0)
                ax.legend(loc="upper right", fontsize=8)
                xform.format_genomic_ticks(ax, locus.breakpoints)

                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()

    print(f"  Saved carrier PDF: {pdf_path}")


# =============================================================================
# Main
# =============================================================================


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot GD CNV detection results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--calls", "-c",
        required=True,
        help="GD CNV calls file (gd_cnv_calls.tsv.gz) produced by gd_cnv_call.py",
    )
    parser.add_argument(
        "--cn-posteriors",
        required=True,
        help="CN posteriors file (cn_posteriors.tsv.gz) with depth values",
    )
    parser.add_argument(
        "--gd-table", "-g",
        required=True,
        help="GD locus definition table (TSV)",
    )
    parser.add_argument(
        "--ploidy-table",
        required=False,
        help="Ploidy estimates TSV (ploidy_estimates.tsv) from gd_cnv_pyro.py. "
             "Columns: sample, contig, median_depth, ploidy. "
             "If not provided, ploidy=2 is assumed for all sample/contig pairs.",
    )
    parser.add_argument(
        "--gtf",
        required=False,
        help="Gene annotation GTF file (can be gzipped)",
    )
    parser.add_argument(
        "--segdup-bed",
        required=False,
        help="Segmental duplication BED file (can be gzipped)",
    )
    parser.add_argument(
        "--gaps-bed",
        required=False,
        help="Reference assembly gaps BED file (can be gzipped).  Regions "
             "are drawn as hatched gray spans on annotation panels.",
    )
    parser.add_argument(
        "--raw-counts",
        required=False,
        help="Raw read-count matrix TSV (can be gzipped).  When provided, "
             "the raw depth profile is superimposed on carrier plots.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        required=True,
        help="Output directory for plots",
    )
    parser.add_argument(
        "--padding",
        type=int,
        default=50000,
        help="Padding around locus boundaries (bp)",
    )
    parser.add_argument(
        "--plot-all-samples",
        action="store_true",
        help="Plot individual sample plots for all samples (not just carriers)",
    )
    parser.add_argument(
        "--sample",
        type=str,
        help="Plot only this specific sample",
    )
    parser.add_argument(
        "--min-gene-label-spacing",
        type=float,
        default=0.05,
        help="Minimum distance between adjacent gene label centres as a fraction "
             "of the plot width (default: 0.05 = 5%%).",
    )
    parser.add_argument(
        "--locus",
        action="append",
        dest="loci",
        metavar="LOCUS",
        help="Restrict plotting to this locus (cluster name).  May be given "
             "multiple times to select several loci.  When omitted all loci "
             "are plotted.",
    )
    parser.add_argument(
        "--flank-scale",
        type=float,
        default=0.20,
        help="Fraction of the total plot width that each flanking region "
             "should occupy.  0.20 means each flank gets 20%% of the "
             "display width, leaving 60%% for the locus body (when both "
             "flanks are present).  Smaller values compress flanks more.  "
             "Use 0 to hide flanks entirely.",
    )
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    log_fh = setup_logging(args.output_dir)

    print(f"Output directory: {args.output_dir}")
    print(f"Flank scale: {args.flank_scale}")

    # Load data
    print("\nLoading data...")

    print(f"  Loading calls: {args.calls}")
    calls_df = pd.read_csv(args.calls, sep="\t", compression="infer")
    print(f"    {len(calls_df)} call records")

    print(f"  Loading CN posteriors: {args.cn_posteriors}")
    cn_posteriors_df = pd.read_csv(args.cn_posteriors, sep="\t", compression="infer")
    print(f"    {len(cn_posteriors_df)} bin-sample records")

    print(f"  Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"    {len(gd_table.loci)} loci")

    # Load ploidy table if provided
    ploidy_df = None
    if args.ploidy_table:
        print(f"  Loading ploidy table: {args.ploidy_table}")
        ploidy_df = pd.read_csv(args.ploidy_table, sep="\t")
        print(f"    {len(ploidy_df)} sample/contig ploidy records")

    # Convert cn_posteriors to depth_df format (bins x samples matrix)
    print("\n  Converting posteriors to depth matrix format...")
    depth_df = cn_posteriors_df.pivot(
        index=["cluster", "chr", "start", "end"],
        columns="sample",
        values="depth",
    ).reset_index()
    depth_df = depth_df.rename(columns={
        "cluster": "Cluster", "chr": "Chr", "start": "Start", "end": "End",
    })
    sample_cols_all = [c for c in depth_df.columns
                       if c not in ["Cluster", "Chr", "Start", "End"]]
    print(f"    {len(depth_df)} bin-rows x {len(sample_cols_all)} samples (across all loci)")

    # Optional annotations
    gtf = None
    if args.gtf:
        print(f"\n  Loading GTF: {args.gtf}")
        gtf = GTFParser(args.gtf)

    segdup = None
    if args.segdup_bed:
        print(f"  Loading segdup BED: {args.segdup_bed}")
        segdup = SegDupAnnotation(args.segdup_bed)

    gaps = None
    if args.gaps_bed:
        print(f"  Loading gaps BED: {args.gaps_bed}")
        gaps = GapsAnnotation(args.gaps_bed)

    # Load raw counts matrix if provided
    raw_counts_df = None
    raw_sample_medians: Dict[str, float] = {}
    if args.raw_counts:
        print(f"\n  Loading raw counts: {args.raw_counts}")
        raw_counts_df = pd.read_csv(args.raw_counts, sep="\t", compression="infer")
        if "#Chr" in raw_counts_df.columns:
            raw_counts_df.rename(columns={"#Chr": "Chr"}, inplace=True)
        raw_counts_df["Start"] = raw_counts_df["Start"].astype(int)
        raw_counts_df["End"] = raw_counts_df["End"].astype(int)
        print(f"    {len(raw_counts_df)} bins")
        autosomal = [str(c) for c in range(1, 23)] + [f"chr{c}" for c in range(1, 23)]
        raw_auto = raw_counts_df[raw_counts_df["Chr"].isin(autosomal)]
        raw_sample_cols = [c for c in raw_counts_df.columns
                          if c not in ("Chr", "Start", "End", "source_file", "Bin")]
        for s in raw_sample_cols:
            med = raw_auto[s].median()
            if med > 0:
                raw_sample_medians[s] = med
        print(f"    Computed raw autosomal medians for {len(raw_sample_medians)} samples")

    # ---- Determine which loci to plot ----
    if args.loci:
        loci_to_plot = {}
        for name in args.loci:
            if name in gd_table.loci:
                loci_to_plot[name] = gd_table.loci[name]
            else:
                print(f"  WARNING: --locus '{name}' not found in GD table, skipping")
        if not loci_to_plot:
            print("ERROR: none of the requested loci were found in the GD table")
            sys.exit(1)
        print(f"  Restricting to {len(loci_to_plot)} of {len(gd_table.loci)} loci")
    else:
        loci_to_plot = gd_table.loci

    # Optionally subset calls to selected loci so summary plots reflect the
    # same scope as the per-locus plots.
    if args.loci:
        plot_calls_df = calls_df[calls_df["cluster"].isin(loci_to_plot)].copy()
    else:
        plot_calls_df = calls_df

    # Build a filtered GDTable-like object for create_carrier_pdf
    class _FilteredGDTable:
        def __init__(self, loci_dict):
            self.loci = loci_dict

    filtered_gd_table = _FilteredGDTable(loci_to_plot)

    # Create summary plots
    print("\nCreating summary plots...")
    plot_carrier_summary(plot_calls_df, args.output_dir)
    plot_confidence_distribution(plot_calls_df, args.output_dir)

    # Create locus overview plots
    print("\nCreating locus overview plots...")
    for cluster, locus in loci_to_plot.items():
        print(f"  Processing {cluster}...")
        plot_locus_overview(
            locus, calls_df, depth_df, gtf, segdup,
            args.output_dir, padding=args.padding,
            ploidy_df=ploidy_df,
            min_gene_label_spacing=args.min_gene_label_spacing,
            raw_counts_df=raw_counts_df,
            raw_sample_medians=raw_sample_medians if raw_sample_medians else None,
            gaps=gaps,
            flank_scale=args.flank_scale,
        )

    # Create individual sample plots
    if args.plot_all_samples or args.sample:
        print("\nCreating individual sample plots...")
        sample_cols = [c for c in depth_df.columns if c not in ['Chr', 'Start', 'End']]
        carriers = calls_df[calls_df["is_carrier"]]

        for cluster, locus in loci_to_plot.items():
            if args.sample:
                samples_to_plot = [args.sample]
            elif args.plot_all_samples:
                samples_to_plot = sample_cols
            else:
                samples_to_plot = carriers[carriers["cluster"] == cluster]["sample"].unique()

            for sample_id in samples_to_plot:
                plot_sample_at_locus(
                    sample_id, locus, calls_df, depth_df, gtf, segdup,
                    args.output_dir, padding=args.padding,
                    min_gene_label_spacing=args.min_gene_label_spacing,
                    gaps=gaps,
                    flank_scale=args.flank_scale,
                )

    # Create carrier PDF
    print("\nCreating carrier PDF...")
    create_carrier_pdf(
        plot_calls_df, depth_df, filtered_gd_table, gtf, segdup,
        args.output_dir, padding=args.padding,
        min_gene_label_spacing=args.min_gene_label_spacing,
        raw_counts_df=raw_counts_df,
        gaps=gaps,
        flank_scale=args.flank_scale,
    )

    print("\n" + "=" * 80)
    print("Plotting complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
