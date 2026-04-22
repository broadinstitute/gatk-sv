"""
Data containers and I/O for depth matrices.

Provides :class:`DepthData` (a torch-tensor wrapper around a bins × samples
depth matrix) and helpers for reading / writing TSV depth files.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import torch

from gatk_sv_ploidy._util import (
    AUTOSOME_NAMES,
    CHR_ORDER,
    DEPTH_SPACES,
    get_sample_columns,
)

logger = logging.getLogger(__name__)


class DepthData:
    """Torch-backed container for a bins × samples normalised depth matrix.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns ``Chr``, ``Start``, ``End``, and one column per
        sample containing normalised read-depth values.
    device : str
        Torch device (``'cpu'`` or ``'cuda'``).
    dtype : torch.dtype
        Floating-point type for the depth tensor.
    subsample_bins : int, optional
        If given, randomly subsample this many bins.
    subsample_samples : int, optional
        If given, randomly subsample this many samples.
    seed : int
        Random seed used when subsampling.
    clamp_threshold : float, optional
        Values above this are clamped.  Set to ``None`` to disable.
    site_data : dict, optional
        Per-site allele count data from :func:`load_site_data`.  Keys are
        ``site_alt``, ``site_total``, ``site_pop_af``, ``site_mask``.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        device: str = "cpu",
        dtype: torch.dtype = torch.float32,
        subsample_bins: Optional[int] = None,
        subsample_samples: Optional[int] = None,
        seed: int = 42,
        clamp_threshold: float = 5.0,
        site_data: Optional[Dict[str, np.ndarray]] = None,
        depth_space: str = "normalized",
    ) -> None:
        sample_cols = get_sample_columns(df)

        # ── optional subsampling ────────────────────────────────────────
        if subsample_bins is not None or subsample_samples is not None:
            rng = np.random.RandomState(seed)
            if subsample_bins is not None and subsample_bins < len(df):
                logger.info("Subsampling %d / %d bins", subsample_bins, len(df))
                idx = np.sort(rng.choice(len(df), subsample_bins, replace=False))
                df = df.iloc[idx].copy()
                if site_data is not None:
                    site_data = {
                        "site_alt": site_data["site_alt"][idx],
                        "site_total": site_data["site_total"][idx],
                        "site_pop_af": site_data["site_pop_af"][idx],
                        "site_mask": site_data["site_mask"][idx],
                    }
            if subsample_samples is not None and subsample_samples < len(sample_cols):
                logger.info(
                    "Subsampling %d / %d samples",
                    subsample_samples,
                    len(sample_cols),
                )
                sel = rng.choice(len(sample_cols), subsample_samples, replace=False)
                sample_cols = [sample_cols[i] for i in sel]
                if site_data is not None:
                    site_data = {
                        "site_alt": site_data["site_alt"][:, :, sel],
                        "site_total": site_data["site_total"][:, :, sel],
                        "site_pop_af": site_data["site_pop_af"],
                        "site_mask": site_data["site_mask"][:, :, sel],
                    }

        # ── sort rows by chromosome order ───────────────────────────────
        df = df.copy()
        df["_order"] = df["Chr"].map(CHR_ORDER)
        df["_orig_pos"] = np.arange(len(df))
        df = df.sort_values(["_order", "Start"])
        sort_perm = df["_orig_pos"].values
        df = df.drop(["_order", "_orig_pos"], axis=1)

        # Reindex site_data to match the sorted row order
        if site_data is not None:
            reindex_keys = {
                "site_alt", "site_total", "site_pop_af", "site_mask",
                "bin_chr", "bin_start", "bin_end",
            }
            site_data = {
                k: site_data[k][sort_perm] if k in reindex_keys
                else site_data[k]
                for k in site_data
            }

            # Cross-check bin coordinates when NPZ includes metadata
            if "bin_chr" in site_data:
                npz_chr = site_data["bin_chr"]
                npz_start = site_data["bin_start"]
                if not np.array_equal(npz_chr, df["Chr"].values):
                    mismatches = np.where(npz_chr != df["Chr"].values)[0]
                    raise ValueError(
                        f"site_data bin_chr does not match depth Chr "
                        f"after sort ({len(mismatches)} mismatches, "
                        f"first at index {mismatches[0]}). "
                        "Regenerate site_data.npz."
                    )
                if not np.array_equal(npz_start, df["Start"].values):
                    mismatches = np.where(
                        npz_start != df["Start"].values
                    )[0]
                    raise ValueError(
                        f"site_data bin_start does not match depth Start "
                        f"after sort ({len(mismatches)} mismatches, "
                        f"first at index {mismatches[0]}). "
                        "Regenerate site_data.npz."
                    )
                logger.info(
                    "Site-data bin coordinates validated against depth bins."
                )

        # ── store metadata arrays ───────────────────────────────────────
        self.chr: np.ndarray = df["Chr"].values
        self.start: np.ndarray = df["Start"].values
        self.end: np.ndarray = df["End"].values
        self.depth_space = str(depth_space).strip().lower()
        if self.depth_space not in DEPTH_SPACES:
            raise ValueError(
                f"Unknown depth_space: {depth_space!r}. Choose one of {DEPTH_SPACES}."
            )
        if "BinLengthBp" in df.columns:
            self.bin_length_bp = np.maximum(
                df["BinLengthBp"].to_numpy(dtype=np.int64),
                1,
            )
        else:
            self.bin_length_bp = np.maximum(
                self.end - self.start,
                1,
            ).astype(np.int64, copy=False)
        self.bin_length_kb = torch.tensor(
            self.bin_length_bp.astype(np.float32) / 1000.0,
            dtype=dtype,
            device=device,
        )
        self.sample_ids: List[str] = sample_cols

        # ── chromosome type tensor (0=autosome, 1=chrX, 2=chrY) ────────
        _chr_type_map = {c: 0 for c in AUTOSOME_NAMES}
        _chr_type_map["chrX"] = 1
        _chr_type_map["chrY"] = 2
        self.chr_type: torch.Tensor = torch.tensor(
            np.array(
                [_chr_type_map.get(c, 0) for c in self.chr], dtype=np.int64,
            ),
            device=device,
        )

        # ── build depth tensor ──────────────────────────────────────────
        depth_values = df[sample_cols].to_numpy(dtype=np.float64)
        if self.depth_space == "raw":
            if np.any(depth_values < 0):
                raise ValueError("Raw count depth input must be non-negative.")
            if not np.allclose(depth_values, np.rint(depth_values)):
                raise ValueError(
                    "Raw count depth input must be integer-valued. "
                    "Run preprocess with --output-space raw or use normalized input."
                )

        depth_matrix = depth_values.astype(np.float32)

        if clamp_threshold is not None:
            n_clamped = int(np.sum(depth_matrix > clamp_threshold))
            if n_clamped > 0:
                logger.info(
                    "Clamping %d values above %.1f", n_clamped, clamp_threshold
                )
                depth_matrix = np.clip(depth_matrix, None, clamp_threshold)

        self.depth = torch.tensor(depth_matrix, dtype=dtype, device=device)
        self.n_bins: int = self.depth.shape[0]
        self.n_samples: int = self.depth.shape[1]

        logger.info("Loaded data: %d bins × %d samples", self.n_bins, self.n_samples)
        logger.info(
            "Depth range: [%.3f, %.3f], mean=%.3f",
            self.depth.min().item(),
            self.depth.max().item(),
            self.depth.mean().item(),
        )

        # ── optional per-site allele data ───────────────────────────────
        if site_data is not None:
            sa = site_data["site_alt"]     # (n_bins, max_sites, n_samples)
            st = site_data["site_total"]   # (n_bins, max_sites, n_samples)
            sp = site_data["site_pop_af"]  # (n_bins, max_sites)
            sm = site_data["site_mask"]    # (n_bins, max_sites, n_samples)

            # ── validate dimensions ─────────────────────────────────────
            if sa.shape[0] != self.n_bins:
                raise ValueError(
                    f"site_data bin count ({sa.shape[0]}) does not match "
                    f"depth bin count ({self.n_bins}). "
                    "Regenerate site_data.npz from the same preprocessed "
                    "depth file."
                )
            if sa.shape[2] != self.n_samples:
                raise ValueError(
                    f"site_data sample count ({sa.shape[2]}) does not match "
                    f"depth sample count ({self.n_samples}). "
                    "Regenerate site_data.npz from the same preprocessed "
                    "depth file."
                )

            self.site_alt = torch.tensor(
                sa.astype(np.float32), dtype=dtype, device=device,
            )
            self.site_total = torch.tensor(
                st.astype(np.float32), dtype=dtype, device=device,
            )
            self.site_pop_af = torch.tensor(
                sp.astype(np.float32), dtype=dtype, device=device,
            )
            self.site_mask = torch.tensor(sm, dtype=torch.bool, device=device)
            self.max_sites: int = sa.shape[1]

            n_informative = int(self.site_mask.any(dim=1).any(dim=1).sum().item())
            total_sites = int(self.site_mask.any(dim=2).sum().item())
            logger.info(
                "Per-site allele data: %d bins with data, %d total sites, "
                "max %d sites/bin",
                n_informative, total_sites, self.max_sites,
            )
        else:
            self.site_alt: Optional[torch.Tensor] = None
            self.site_total: Optional[torch.Tensor] = None
            self.site_pop_af: Optional[torch.Tensor] = None
            self.site_mask: Optional[torch.Tensor] = None
            self.max_sites: int = 0


# ── reading / writing helpers ───────────────────────────────────────────────


def read_depth_tsv(path: str) -> pd.DataFrame:
    """Read a single depth TSV file (optionally gzipped).

    The file is expected to have columns ``#Chr`` (or ``Chr``), ``Start``,
    ``End``, followed by one column per sample.  A ``Bin`` index column and a
    ``source_file`` column are added.

    Args:
        path: Path to the TSV file.

    Returns:
        DataFrame with metadata and per-sample depth columns.

    Raises:
        Exception: Propagated from :func:`pandas.read_csv`.
    """
    logger.info("Loading: %s", path)
    df = pd.read_csv(path, sep="\t", compression="infer")

    # Normalise chromosome column name
    if "#Chr" in df.columns:
        df = df.rename(columns={"#Chr": "Chr"})

    # Add provenance info
    bin_chr = df["Chr"].astype(str)
    bin_start = df["Start"].astype(str)
    bin_end = df["End"].astype(str)
    df["Bin"] = bin_chr.str.cat(bin_start, sep=":").str.cat(bin_end, sep="-")
    df = df.set_index("Bin")
    df["source_file"] = str(Path(path).parent.parent.name)
    return df


def load_site_data(path: str) -> Dict[str, np.ndarray]:
    """Load per-site allele data from an ``.npz`` archive.

    The archive is expected to contain the keys ``site_alt``, ``site_total``,
    ``site_pop_af``, ``site_mask`` (as produced by
    :func:`~gatk_sv_ploidy.preprocess.build_per_site_data`).  Bin-coordinate
    metadata (``bin_chr``, ``bin_start``, ``bin_end``) is also loaded when
    present; :class:`DepthData` uses these to validate alignment after
    sorting.

    Args:
        path: Path to the ``.npz`` file.

    Returns:
        Dictionary with NumPy arrays ready for :class:`DepthData`.
    """
    logger.info("Loading per-site allele data: %s", path)
    npz = np.load(path, allow_pickle=True)
    result = {
        "site_alt": npz["site_alt"],
        "site_total": npz["site_total"],
        "site_pop_af": npz["site_pop_af"],
        "site_mask": npz["site_mask"],
    }
    # Include bin coordinate metadata when available (for alignment
    # validation inside DepthData).
    for key in ("bin_chr", "bin_start", "bin_end"):
        if key in npz:
            result[key] = npz[key]
    logger.info(
        "  shape: bins=%d, max_sites=%d, samples=%d",
        result["site_alt"].shape[0],
        result["site_alt"].shape[1],
        result["site_alt"].shape[2],
    )
    return result
