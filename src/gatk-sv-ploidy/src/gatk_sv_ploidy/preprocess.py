"""
Preprocess subcommand — read, normalise, and filter depth data.

Reads a raw depth TSV, normalises each sample by its autosomal median so that
diploid depth ≈ 2.0, filters low-quality bins by median/MAD thresholds, and
writes the preprocessed matrix to disk.
"""

from __future__ import annotations

import argparse
import logging
import os
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from gatk_sv_ploidy._util import get_sample_columns
from gatk_sv_ploidy.data import read_depth_tsv

logger = logging.getLogger(__name__)


# ── normalisation ───────────────────────────────────────────────────────────


def normalise_depth(df: pd.DataFrame) -> pd.DataFrame:
    """Normalise per-sample depth so that autosomal CN = 2 corresponds to 2.0.

    For each sample the median depth across autosomal bins is computed and all
    bins (including sex chromosomes) are scaled by ``2 / median``.

    Args:
        df: DataFrame with ``Chr`` metadata column and per-sample depth columns.

    Returns:
        Copy of *df* with normalised depth values.
    """
    sample_cols = get_sample_columns(df)
    autosome_mask = ~df["Chr"].isin(["chrX", "chrY"])
    medians = np.median(df.loc[autosome_mask, sample_cols].values, axis=0)

    logger.info(
        "Autosomal medians: min=%.3f, max=%.3f, mean=%.3f",
        medians.min(),
        medians.max(),
        medians.mean(),
    )

    df = df.copy()
    df[sample_cols] = 2.0 * df[sample_cols].values / medians[np.newaxis, :]
    return df


# ── bin quality filtering ───────────────────────────────────────────────────


def filter_low_quality_bins(
    df: pd.DataFrame,
    *,
    autosome_median_min: float = 1.0,
    autosome_median_max: float = 3.0,
    autosome_mad_max: float = 2.0,
    chrX_median_min: float = 0.0,
    chrX_median_max: float = 3.0,
    chrX_mad_max: float = 2.0,
    chrY_median_min: float = 0.0,
    chrY_median_max: float = 3.0,
    chrY_mad_max: float = 2.0,
    min_bins_per_chr: int = 10,
) -> pd.DataFrame:
    """Filter bins whose cross-sample statistics fall outside thresholds.

    Thresholds are applied separately for autosomes, chrX, and chrY because
    expected depth distributions differ across chromosome types.

    Args:
        df: Normalised depth DataFrame.
        autosome_median_min/max: Median depth range for autosomal bins.
        autosome_mad_max: Maximum MAD for autosomal bins.
        chrX_median_min/max: Median depth range for chrX bins.
        chrX_mad_max: Maximum MAD for chrX bins.
        chrY_median_min/max: Median depth range for chrY bins.
        chrY_mad_max: Maximum MAD for chrY bins.
        min_bins_per_chr: Raise if any chromosome has fewer bins after
            filtering.

    Returns:
        Filtered DataFrame.

    Raises:
        ValueError: If a chromosome has fewer than *min_bins_per_chr* bins.
    """
    sample_cols = get_sample_columns(df)
    depths = df[sample_cols].values
    medians = np.median(depths, axis=1)
    mads = stats.median_abs_deviation(depths, axis=1)

    logger.info("Starting bins: %d", len(df))

    keep = np.ones(len(df), dtype=bool)

    # -- per-chromosome-type thresholds --
    _thresholds = {
        "autosome": (
            ~df["Chr"].isin(["chrX", "chrY"]),
            autosome_median_min,
            autosome_median_max,
            autosome_mad_max,
        ),
        "chrX": (
            df["Chr"] == "chrX",
            chrX_median_min,
            chrX_median_max,
            chrX_mad_max,
        ),
        "chrY": (
            df["Chr"] == "chrY",
            chrY_median_min,
            chrY_median_max,
            chrY_mad_max,
        ),
    }

    for label, (mask, med_min, med_max, mad_max) in _thresholds.items():
        if not mask.any():
            continue
        ok = (medians >= med_min) & (medians <= med_max) & (mads <= mad_max)
        n_before = int(mask.sum())
        keep[mask.values] &= ok[mask.values]
        n_after = int((mask & keep).sum())
        logger.info(
            "%s: median [%.1f, %.1f], MAD ≤ %.1f → %d / %d bins kept",
            label,
            med_min,
            med_max,
            mad_max,
            n_after,
            n_before,
        )

    df_out = df[keep].copy()
    logger.info("Bins after filtering: %d (removed %d)", len(df_out), len(df) - len(df_out))

    # Verify sufficient bins per chromosome
    counts = df_out.groupby("Chr").size()
    bad = counts[counts < min_bins_per_chr]
    if len(bad) > 0:
        msg = "; ".join(f"{c}={n}" for c, n in bad.items())
        raise ValueError(
            f"Insufficient bins after filtering ({msg}). "
            "Consider --skip-bin-filter or relaxing thresholds."
        )

    return df_out

# ── allele fraction preprocessing (per-site joint model) ─────────────────────


def read_site_depth_tsv(
    path: str, stride: int = 1, position_stride: int = 0,
) -> pd.DataFrame:
    """Read a GATK CollectSVEvidence site-depth file.

    SD files are headerless, tab-delimited, optionally bgzipped, with columns:

    1. contig (str)
    2. position (int, 0-based)
    3. sample (str)
    4. A depth (int)
    5. C depth (int)
    6. G depth (int)
    7. T depth (int)

    Args:
        path: Path to ``.sd.txt.gz`` or plain SD file.
        stride: Keep every *stride*-th row by file-row index (1 = keep all).
            Ignored when *position_stride* > 0.
        position_stride: Keep positions where ``position % position_stride == 0``.
            This is slower than *stride* (the full file must be decompressed)
            but guarantees that every sample keeps the **same** set of genomic
            positions, which is critical for multi-sample allele-fraction
            analysis.  When > 0, *stride* is ignored.

    Returns:
        DataFrame with columns ``contig``, ``position``, ``sample``,
        ``A``, ``C``, ``G``, ``T``.
    """
    _dtypes = {
        "contig": str,
        "position": np.int64,
        "sample": str,
        "A": np.int32,
        "C": np.int32,
        "G": np.int32,
        "T": np.int32,
    }
    _names = ["contig", "position", "sample", "A", "C", "G", "T"]

    position_stride = max(0, int(position_stride))
    stride = max(1, int(stride))

    if position_stride > 1:
        # Position-based striding: read in chunks and keep only rows whose
        # genomic position is divisible by position_stride.  This ensures
        # cross-sample consistency at the cost of a full file scan.
        logger.info(
            "Reading site depth file: %s (position_stride=%d)", path, position_stride,
        )
        chunks: list[pd.DataFrame] = []
        reader = pd.read_csv(
            path, sep="\t", compression="infer", header=None,
            names=_names, dtype=_dtypes, chunksize=2_000_000,
        )
        for chunk in reader:
            filtered = chunk[chunk["position"] % position_stride == 0]
            if len(filtered) > 0:
                chunks.append(filtered)
        df = (
            pd.concat(chunks, ignore_index=True)
            if chunks
            else pd.DataFrame(columns=_names)
        )
    else:
        logger.info("Reading site depth file: %s (stride=%d)", path, stride)
        skiprows = (lambda i: i % stride != 0) if stride > 1 else None
        df = pd.read_csv(
            path, sep="\t", compression="infer", header=None,
            names=_names, dtype=_dtypes, skiprows=skiprows,
        )

    logger.info(
        "  %d sites for %d sample(s)",
        len(df),
        df["sample"].nunique(),
    )
    return df


def _site_alt_total(
    a: np.ndarray, c: np.ndarray, g: np.ndarray, t: np.ndarray,
    ref_base: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute alt allele count and total depth given reference base identity.

    For each position the reference allele is extracted based on *ref_base*
    (one of ``'A'``, ``'C'``, ``'G'``, ``'T'``).  The alt count is
    ``total - ref_count``.

    Args:
        a, c, g, t: Per-base read counts (1-D integer arrays).
        ref_base: 1-D array of single-character strings indicating the
            reference allele at each position.

    Returns:
        ``(alt, total)`` integer arrays with the same length as the inputs.
    """
    counts = np.stack([a, c, g, t], axis=-1)       # (N, 4)
    total = counts.sum(axis=-1)                      # (N,)

    base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}
    ref_idx = np.array([base_to_idx.get(b, 0) for b in ref_base], dtype=np.intp)
    ref_count = counts[np.arange(len(counts)), ref_idx]
    alt = total - ref_count
    return alt, total


def _site_minor_total(a: np.ndarray, c: np.ndarray, g: np.ndarray, t: np.ndarray):
    """Compute minor allele count and total depth from per-base counts.

    The major allele is defined as the base with the most reads.  The minor
    allele count is ``total - major``.

    Returns:
        ``(minor, total)`` integer arrays with the same length as the inputs.
    """
    counts = np.stack([a, c, g, t], axis=-1)
    total = counts.sum(axis=-1)
    major = counts.max(axis=-1)
    minor = total - major
    return minor, total


def read_known_sites(path: str) -> pd.DataFrame:
    """Read a known-sites file with population allele frequencies.

    Accepts either a VCF (``#CHROM`` header) or a simple 4-column TSV
    (``contig``, ``position``, ``ref``, ``pop_af``).  For VCFs the ``AF``
    INFO field is extracted; for multi-allelic records only the first ALT AF
    is used.  Positions are stored as 0-based.

    Args:
        path: Path to known-sites file (optionally gzipped).

    Returns:
        DataFrame with columns ``contig``, ``position`` (0-based), ``ref``,
        ``pop_af``.
    """
    logger.info("Reading known sites: %s", path)
    # Peek at the first line to detect format
    import gzip
    _open = gzip.open if path.endswith(".gz") else open
    with _open(path, "rt") as fh:
        first_line = ""
        for line in fh:
            if not line.startswith("##"):
                first_line = line
                break

    if first_line.startswith("#CHROM") or first_line.startswith("CHROM"):
        # VCF format — parse manually for speed
        rows = []
        with _open(path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t", 8)
                chrom = parts[0]
                pos = int(parts[1]) - 1   # VCF is 1-based → 0-based
                ref = parts[3]
                info = parts[7]
                # Extract AF from INFO
                af_val = 0.5
                for field in info.split(";"):
                    if field.startswith("AF="):
                        af_str = field[3:].split(",")[0]
                        af_val = float(af_str)
                        break
                rows.append((chrom, pos, ref, af_val))
        df = pd.DataFrame(rows, columns=["contig", "position", "ref", "pop_af"])
    else:
        # Simple TSV format
        df = pd.read_csv(
            path, sep="\t", compression="infer",
            names=["contig", "position", "ref", "pop_af"],
            dtype={"contig": str, "position": np.int64,
                   "ref": str, "pop_af": np.float64},
        )

    logger.info("  %d known sites loaded", len(df))
    return df


def build_per_site_data(
    sd_paths: List[str],
    bins_df: pd.DataFrame,
    known_sites_df: Optional[pd.DataFrame] = None,
    *,
    min_site_depth: int = 10,
    max_sites_per_bin: int = 50,
    sd_stride: int = 100,
    seed: int = 42,
) -> dict:
    """Build per-site allele-count arrays for the joint genotype/CN model.

    For every known SNP that falls inside a genomic bin, the alt and total
    allele counts are stored per-sample.  Sites are grouped by bin and padded
    to *max_sites_per_bin*.

    When *known_sites_df* is ``None``, falls back to a naive mode where every
    position in the SD files is treated as a known site with ``pop_af=0.5``
    and the alt count is defined as ``total - major_allele``.

    Args:
        sd_paths: Paths to per-sample ``.sd.txt.gz`` files.
        bins_df: Preprocessed depth DataFrame with ``Chr``, ``Start``, ``End``
            metadata columns.
        known_sites_df: DataFrame with ``contig``, ``position`` (0-based),
            ``ref``, ``pop_af`` (from :func:`read_known_sites`).  If ``None``,
            every observed position is used with ``pop_af = 0.5``.
        min_site_depth: Minimum total depth at a site to include.
        max_sites_per_bin: Pad/subsample to this many sites per bin.
        sd_stride: Keep every *sd_stride*-th row when reading SD files
            (1 = keep all, 100 = 100× downsample).  Ignored when
            *known_sites_df* is provided (all matching sites are kept).
        seed: Random seed for reproducible subsampling.

    Returns:
        Dictionary with the following NumPy arrays (ready for
        :func:`numpy.savez`):

        - ``site_alt``: ``int32 (n_bins, max_sites, n_samples)`` alt counts
        - ``site_total``: ``int32 (n_bins, max_sites, n_samples)`` total counts
        - ``site_pop_af``: ``float32 (n_bins, max_sites)`` population AFs
        - ``site_mask``: ``bool (n_bins, max_sites, n_samples)`` valid entries
        - ``sample_ids``: ``str (n_samples,)``
        - ``bin_chr``: ``str (n_bins,)``
        - ``bin_start``: ``int64 (n_bins,)``
        - ``bin_end``: ``int64 (n_bins,)``
    """
    sample_cols = get_sample_columns(bins_df)
    n_bins = len(bins_df)
    n_samples = len(sample_cols)
    sample_to_idx = {s: i for i, s in enumerate(sample_cols)}

    have_known = known_sites_df is not None

    # Build bin-lookup per chromosome
    bin_lookup: dict = {}
    for chrom in bins_df["Chr"].unique():
        chrom_mask = (bins_df["Chr"] == chrom).values
        chrom_indices = np.where(chrom_mask)[0]
        starts = bins_df["Start"].values[chrom_mask]
        ends = bins_df["End"].values[chrom_mask]
        bin_lookup[chrom] = (chrom_indices, starts, ends)

    # Build known-site lookup per chromosome as sorted arrays for
    # vectorised searchsorted membership tests (no per-element dict lookups).
    # known_arrays[chrom] = (sorted_positions, ref_bases, pop_afs)
    known_arrays: dict = {}
    if have_known:
        for chrom, grp in known_sites_df.groupby("contig"):
            order = np.argsort(grp["position"].values)
            known_arrays[chrom] = (
                grp["position"].values[order].astype(np.int64),
                grp["ref"].str.upper().values[order],
                grp["pop_af"].values[order].astype(np.float32),
            )
        n_known = sum(len(v[0]) for v in known_arrays.values())
        logger.info(
            "Known-site lookup: %d chromosomes, %d total sites",
            len(known_arrays), n_known,
        )

    # Intermediate: collect per-bin lists of (site_pos, pop_af, per-sample alt, per-sample total)
    # For each bin we accumulate site-level data across all SD files
    # site_data[global_bin_idx] = {pos: {"pop_af": float, "alt": [n_samples], "total": [n_samples]}}
    site_data: List[dict] = [{} for _ in range(n_bins)]

    total_sites_used = 0

    # When known sites are provided the searchsorted filter is fast
    # enough; only downsample in fallback (no known sites) mode.
    effective_stride = 1 if have_known else sd_stride

    for sd_path in sd_paths:
        if have_known:
            # Known-sites mode: read at full resolution (fast searchsorted
            # filter handles the reduction).
            sd_df = read_site_depth_tsv(sd_path, stride=1)
        else:
            # No known sites: use position-based striding so that every
            # sample keeps the same set of genomic positions.  This is
            # slower than file-row stride but critical for multi-sample
            # allele-fraction analysis where disjoint positions produce
            # meaningless results.
            sd_df = read_site_depth_tsv(sd_path, position_stride=effective_stride)

        for sample_id in sd_df["sample"].unique():
            if sample_id not in sample_to_idx:
                logger.warning(
                    "Sample %s in SD file not in depth data, skipping",
                    sample_id,
                )
                continue

            si = sample_to_idx[sample_id]
            smask = sd_df["sample"].values == sample_id
            s_contigs = sd_df["contig"].values[smask]
            s_positions = sd_df["position"].values[smask].astype(np.int64)
            s_a = sd_df["A"].values[smask]
            s_c = sd_df["C"].values[smask]
            s_g = sd_df["G"].values[smask]
            s_t = sd_df["T"].values[smask]

            for chrom, (chrom_bin_idx, starts, ends) in bin_lookup.items():
                cmask = s_contigs == chrom
                if not cmask.any():
                    continue

                c_pos = s_positions[cmask]
                c_a = s_a[cmask]
                c_c = s_c[cmask]
                c_g = s_g[cmask]
                c_t = s_t[cmask]

                # ── Assign sites to bins (vectorised) ──
                bi = np.searchsorted(starts, c_pos, side="right") - 1
                bi_safe = np.clip(bi, 0, len(starts) - 1)
                in_bin = (
                    (bi >= 0) & (bi < len(starts))
                    & (c_pos < ends[bi_safe])
                )
                c_pos, c_a, c_c, c_g, c_t = (
                    c_pos[in_bin], c_a[in_bin], c_c[in_bin],
                    c_g[in_bin], c_t[in_bin],
                )
                g_bi_arr = chrom_bin_idx[bi[in_bin]]
                if len(c_pos) == 0:
                    continue

                # ── Vectorised depth filter ──
                totals = (c_a + c_c + c_g + c_t).astype(np.int64)
                depth_ok = totals >= min_site_depth
                c_pos = c_pos[depth_ok]
                c_a, c_c, c_g, c_t = (
                    c_a[depth_ok], c_c[depth_ok],
                    c_g[depth_ok], c_t[depth_ok],
                )
                g_bi_arr = g_bi_arr[depth_ok]
                totals = totals[depth_ok]
                if len(c_pos) == 0:
                    continue

                # ── Known-sites filter / alt-count computation (vectorised) ──
                if have_known:
                    kp, kr, kaf = known_arrays.get(
                        chrom, (np.empty(0, np.int64),
                                np.empty(0, dtype="U1"),
                                np.empty(0, np.float32)),
                    )
                    if len(kp) == 0:
                        continue
                    ins = np.searchsorted(kp, c_pos)
                    ins_safe = np.clip(ins, 0, max(len(kp) - 1, 0))
                    match = kp[ins_safe] == c_pos

                    c_pos = c_pos[match]
                    c_a, c_c, c_g, c_t = (
                        c_a[match], c_c[match], c_g[match], c_t[match],
                    )
                    g_bi_arr = g_bi_arr[match]
                    totals = totals[match]
                    ref_bases = kr[ins_safe[match]]
                    pop_afs = kaf[ins_safe[match]]

                    ref_counts = np.where(
                        ref_bases == "A", c_a,
                        np.where(ref_bases == "C", c_c,
                                 np.where(ref_bases == "G", c_g, c_t)),
                    )
                    alts = totals - ref_counts
                else:
                    major = np.maximum(
                        np.maximum(c_a, c_c), np.maximum(c_g, c_t),
                    )
                    alts = totals - major
                    pop_afs = np.full(len(c_pos), 0.5, dtype=np.float32)

                if len(c_pos) == 0:
                    continue

                # ── Accumulate into per-bin dicts ──
                n_sites_chrom = len(c_pos)
                total_sites_used += n_sites_chrom
                for i in range(n_sites_chrom):
                    pos = int(c_pos[i])
                    g_bi = int(g_bi_arr[i])
                    entry = site_data[g_bi]
                    if pos not in entry:
                        entry[pos] = {
                            "pop_af": float(pop_afs[i]),
                            "alt": np.zeros(n_samples, dtype=np.int32),
                            "total": np.zeros(n_samples, dtype=np.int32),
                        }
                    entry[pos]["alt"][si] = int(alts[i])
                    entry[pos]["total"][si] = int(totals[i])

                if n_sites_chrom > 0:
                    logger.debug(
                        "  %s %s: %d sites accumulated",
                        sample_id, chrom, n_sites_chrom,
                    )

    # ── pack into padded arrays ──────────────────────────────────────────
    site_alt = np.zeros((n_bins, max_sites_per_bin, n_samples), dtype=np.int32)
    site_total = np.zeros((n_bins, max_sites_per_bin, n_samples), dtype=np.int32)
    site_pop_af = np.zeros((n_bins, max_sites_per_bin), dtype=np.float32)
    site_mask = np.zeros((n_bins, max_sites_per_bin, n_samples), dtype=bool)

    bins_with_data = 0
    for g_bi in range(n_bins):
        entries = site_data[g_bi]
        if not entries:
            continue
        bins_with_data += 1

        positions = sorted(entries.keys())
        if len(positions) > max_sites_per_bin:
            # Prioritise positions with data for the most samples.
            coverage = np.array(
                [int(np.sum(entries[p]["total"] > 0)) for p in positions],
            )
            # Stable argsort descending; among ties, keep original order
            order = np.argsort(-coverage, kind="mergesort")
            positions = sorted([positions[i] for i in order[:max_sites_per_bin]])

        for slot, pos in enumerate(positions):
            e = entries[pos]
            site_alt[g_bi, slot, :] = e["alt"]
            site_total[g_bi, slot, :] = e["total"]
            site_pop_af[g_bi, slot] = e["pop_af"]
            site_mask[g_bi, slot, :] = e["total"] > 0

    logger.info(
        "Per-site data: %d total site-sample entries, %d / %d bins with data, "
        "padded to %d sites/bin",
        total_sites_used, bins_with_data, n_bins, max_sites_per_bin,
    )

    return {
        "site_alt": site_alt,
        "site_total": site_total,
        "site_pop_af": site_pop_af,
        "site_mask": site_mask,
        "sample_ids": np.array(sample_cols, dtype=object),
        "bin_chr": bins_df["Chr"].values.astype(object),
        "bin_start": bins_df["Start"].values.astype(np.int64),
        "bin_end": bins_df["End"].values.astype(np.int64),
    }

# ── CLI ─────────────────────────────────────────────────────────────────────


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the preprocess subcommand."""
    p = argparse.ArgumentParser(
        description="Preprocess depth data for aneuploidy inference",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument(
        "-i", "--input", required=True,
        help="Input TSV file with raw read-depth (bins × samples)",
    )
    p.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory",
    )
    p.add_argument(
        "--viable-only", action="store_true", default=False,
        help="Subset to chr13, chr18, chr21, chrX, chrY only",
    )
    p.add_argument(
        "--skip-bin-filter", action="store_true", default=False,
        help="Skip bin quality filtering",
    )

    # Filter thresholds
    g = p.add_argument_group("bin-filter thresholds")
    g.add_argument("--autosome-median-min", type=float, default=1.0)
    g.add_argument("--autosome-median-max", type=float, default=3.0)
    g.add_argument("--autosome-mad-max", type=float, default=2.0)
    g.add_argument("--chrX-median-min", type=float, default=0.0)
    g.add_argument("--chrX-median-max", type=float, default=3.0)
    g.add_argument("--chrX-mad-max", type=float, default=2.0)
    g.add_argument("--chrY-median-min", type=float, default=0.0)
    g.add_argument("--chrY-median-max", type=float, default=3.0)
    g.add_argument("--chrY-mad-max", type=float, default=2.0)

    # Allele fraction arguments
    a = p.add_argument_group("allele fraction (site depth)")
    a.add_argument(
        "--site-depth-list", default=None,
        help="Text file listing paths to per-sample SD files (one per line)",
    )
    a.add_argument(
        "--site-depth", nargs="*", default=None,
        help="One or more per-sample SD file paths directly on the command line",
    )
    a.add_argument(
        "--known-sites", default=None,
        help="VCF or TSV of known SNP sites with population AFs.  "
             "If omitted, all observed positions are used with pop_af=0.5",
    )
    a.add_argument("--min-site-depth", type=int, default=10,
                   help="Minimum total depth at a site to include")
    a.add_argument("--max-sites-per-bin", type=int, default=50,
                   help="Maximum sites per bin (pad/subsample to this)")
    a.add_argument("--sd-stride", type=int, default=1,
                   help="Keep every Nth row from SD files (100 = 100x "
                        "downsample).  Ignored when --known-sites is given.")

    return p.parse_args()


def main() -> None:
    """Entry point for ``gatk-sv-ploidy preprocess``."""
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Read
    df = read_depth_tsv(args.input)

    # 2. Optionally subset to viable trisomy chromosomes
    if args.viable_only:
        viable = {"chr13", "chr18", "chr21", "chrX", "chrY"}
        n_before = len(df)
        df = df[df["Chr"].isin(viable)]
        logger.info(
            "Viable-only filter: %d → %d bins (%s)",
            n_before,
            len(df),
            sorted(df["Chr"].unique()),
        )

    # 3. Normalise
    df = normalise_depth(df)

    # 4. Filter
    if args.skip_bin_filter:
        logger.info("Skipping bin quality filtering (--skip-bin-filter)")
    else:
        df = filter_low_quality_bins(
            df,
            autosome_median_min=args.autosome_median_min,
            autosome_median_max=args.autosome_median_max,
            autosome_mad_max=args.autosome_mad_max,
            chrX_median_min=args.chrX_median_min,
            chrX_median_max=args.chrX_median_max,
            chrX_mad_max=args.chrX_mad_max,
            chrY_median_min=args.chrY_median_min,
            chrY_median_max=args.chrY_median_max,
            chrY_mad_max=args.chrY_mad_max,
        )

    # 5. Write preprocessed depth
    out_path = os.path.join(args.output_dir, "preprocessed_depth.tsv")
    df.to_csv(out_path, sep="\t")
    logger.info("Preprocessed depth written to %s", out_path)

    # 6. Build per-site allele data from site-depth files (optional)
    sd_paths: List[str] = []
    if args.site_depth_list:
        with open(args.site_depth_list) as fh:
            sd_paths.extend(line.strip() for line in fh if line.strip())
    if args.site_depth:
        sd_paths.extend(args.site_depth)

    if sd_paths:
        known_sites_df = None
        if args.known_sites:
            known_sites_df = read_known_sites(args.known_sites)

        logger.info("Building per-site allele data from %d SD file(s) …", len(sd_paths))
        site_arrays = build_per_site_data(
            sd_paths,
            df,
            known_sites_df=known_sites_df,
            min_site_depth=args.min_site_depth,
            max_sites_per_bin=args.max_sites_per_bin,
            sd_stride=args.sd_stride,
        )
        site_data_path = os.path.join(args.output_dir, "site_data.npz")
        np.savez_compressed(site_data_path, **site_arrays)
        logger.info("Per-site allele data written to %s", site_data_path)
    else:
        logger.info("No site-depth files provided; skipping allele data")


if __name__ == "__main__":
    main()
