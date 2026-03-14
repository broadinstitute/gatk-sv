"""
Synthesize GD CNVs into read-depth matrices.

Spikes simulated genomic-disorder CNV signal into high- and low-resolution
binned read-count matrices.  For each sample selected to carry a GD event
(controlled by ``--gd-probability``), a single randomly chosen GD entry is
assigned and the depth counts in bins overlapping the GD locus are scaled.
When a preprocess-generated ploidy table is provided, the copy-step is made
relative to the sample's baseline contig ploidy:

* **Deletions** – counts multiplied by $\frac{p - 1}{p}$
* **Duplications** – counts multiplied by $\frac{p + 1}{p}$

where $p$ is the baseline ploidy for that sample/contig.  Without a ploidy
table, synthesize falls back to the historical diploid assumptions
(``0.5`` for deletions and ``1.5`` for duplications). Counts are rounded to
the nearest integer after scaling.

A truth table mapping sample IDs to their assigned GD IDs is written
alongside the modified matrices.

Usage::

    gatk-sv-gd synthesize \\
        --lo-res-counts lowres.rd.txt.gz \\
        --hi-res-counts highres.rd.txt.gz \\
        --gd-table gd_loci.tsv \\
        --output-dir out/ \\
        --seed 42

CLI entry-point: :func:`main`.
"""

import argparse
import csv
import gzip
import multiprocessing
import os
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pysam
from tqdm import tqdm

from gatk_sv_gd.models import GDLocus, GDTable


# ---------------------------------------------------------------------------
# Region parsing (shared with preprocess)
# ---------------------------------------------------------------------------

def _parse_region(region_str: str) -> Tuple[str, Optional[int], Optional[int]]:
    """Parse ``chr1:3000-4000`` or ``chr1`` into (chrom, start|None, end|None)."""
    if ":" in region_str:
        chrom, coords = region_str.split(":", 1)
        parts = coords.replace(",", "").split("-")
        if len(parts) != 2:
            raise ValueError(
                f"Invalid region format '{region_str}': expected chrom:start-end"
            )
        return chrom, int(parts[0]), int(parts[1])
    return region_str, None, None


def _gd_entry_overlaps_regions(
    entry: dict,
    chrom: str,
    regions: List[Tuple[str, Optional[int], Optional[int]]],
) -> bool:
    """Return True if a GD entry overlaps any parsed region."""
    for r_chrom, r_start, r_end in regions:
        if chrom != r_chrom:
            continue
        if r_start is None:
            return True
        if entry["start_GRCh38"] < r_end and entry["end_GRCh38"] > r_start:
            return True
    return False


def _open_textfile(filepath: str, mode: str = "rt"):
    """Open plain-text or gzipped text files transparently."""
    return gzip.open(filepath, mode) if filepath.endswith(".gz") else open(filepath, mode)


def _detect_baf_columns(filepath: str) -> Tuple[Optional[int], List[str]]:
    """Detect whether a BAF file has a header and return logical columns."""
    with _open_textfile(filepath, "rt") as handle:
        first_line = handle.readline().strip()

    if not first_line:
        return None, ["Chr", "Pos", "BAF", "Sample"]

    parts = first_line.split("\t")
    if len(parts) < 4:
        raise ValueError(
            f"BAF table {filepath} must have at least 4 tab-delimited columns"
        )

    try:
        int(parts[1])
        float(parts[2])
    except ValueError:
        return 0, parts[:4]

    return None, ["Chr", "Pos", "BAF", "Sample"]


def _read_count_sample_ids(filepath: str) -> List[str]:
    """Read sample column names from a count matrix header."""
    with _open_textfile(filepath, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
    meta_cols = {"#Chr", "Chr", "Start", "End", "source_file", "Bin"}
    return [col for col in header if col not in meta_cols]


def _build_gd_lookup(gd_table: GDTable) -> Dict[str, Tuple[str, dict]]:
    """Map ``GD_ID`` to its source ``(chrom, entry)`` pair."""
    lookup: Dict[str, Tuple[str, dict]] = {}
    for locus in gd_table.get_all_loci().values():
        for entry in locus.gd_entries:
            gd_id = str(entry["GD_ID"])
            if gd_id in lookup:
                raise ValueError(f"Duplicate GD_ID in GD table: {gd_id}")
            lookup[gd_id] = (locus.chrom, entry)
    return lookup


def _load_truth_assignments(
    truth_table_path: str,
    gd_lookup: Dict[str, Tuple[str, dict]],
) -> Dict[str, Tuple[str, dict]]:
    """Load synthesize-format truth assignments from an existing table."""
    with _open_textfile(truth_table_path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "GD_ID"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Truth table must contain synthesize-format columns: sample_id, GD_ID"
            )

        assignments: Dict[str, Tuple[str, dict]] = {}
        for row in reader:
            sample_id = str(row["sample_id"]).strip()
            gd_id = str(row["GD_ID"]).strip()
            if not sample_id or not gd_id:
                continue
            if gd_id not in gd_lookup:
                raise ValueError(f"GD_ID '{gd_id}' from truth table not found in GD table")
            existing = assignments.get(sample_id)
            if existing is not None and str(existing[1]["GD_ID"]) != gd_id:
                raise ValueError(
                    f"Sample '{sample_id}' has multiple GD assignments in truth table"
                )
            assignments[sample_id] = gd_lookup[gd_id]
    return assignments


def _read_baf_sample_ids(filepath: str) -> List[str]:
    """Scan a BAF table and return the distinct sample IDs it contains."""
    header, column_names = _detect_baf_columns(filepath)
    sample_idx = column_names.index("Sample")
    sample_ids: Set[str] = set()
    with _open_textfile(filepath, "rt") as handle:
        for line_no, line in enumerate(handle):
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= sample_idx:
                continue
            if header == 0 and line_no == 0:
                continue
            sample_id = str(fields[sample_idx]).strip()
            if sample_id:
                sample_ids.add(sample_id)
    return sorted(sample_ids)


def _discover_sample_ids(
    lo_res_counts: Optional[str],
    hi_res_counts: Optional[str],
    baf_table: Optional[str],
) -> List[str]:
    """Discover the batch sample set from the first available input table."""
    if lo_res_counts:
        return _read_count_sample_ids(lo_res_counts)
    if hi_res_counts:
        return _read_count_sample_ids(hi_res_counts)
    if baf_table:
        return _read_baf_sample_ids(baf_table)
    return []


def _load_ploidy_lookup(ploidy_table_path: str) -> Dict[Tuple[str, str], int]:
    """Load a ``{(sample, contig): ploidy}`` lookup from preprocess output."""
    with open(ploidy_table_path, "r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "contig", "ploidy"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                "Ploidy table must contain columns: sample, contig, ploidy"
            )

        lookup: Dict[Tuple[str, str], int] = {}
        for row in reader:
            sample_id = str(row["sample"]).strip()
            contig = str(row["contig"]).strip()
            ploidy_str = str(row["ploidy"]).strip()
            if not sample_id or not contig or not ploidy_str:
                continue
            lookup[(sample_id, contig)] = int(float(ploidy_str))
    return lookup


def _resolve_sample_contig_ploidy(
    sample_id: str,
    chrom: str,
    ploidy_lookup: Optional[Dict[Tuple[str, str], int]],
    default_ploidy: int = 2,
) -> int:
    """Return baseline ploidy for a sample/contig, tolerating chr-prefix mismatches."""
    if not ploidy_lookup:
        return default_ploidy

    keys_to_try = [(str(sample_id), str(chrom))]
    if str(chrom).startswith("chr"):
        keys_to_try.append((str(sample_id), str(chrom)[3:]))
    else:
        keys_to_try.append((str(sample_id), f"chr{chrom}"))

    for key in keys_to_try:
        if key in ploidy_lookup:
            return max(1, int(ploidy_lookup[key]))
    return default_ploidy


def _resolve_event_multiplier(
    svtype: str,
    baseline_ploidy: int,
    fallback_del_multiplier: float,
    fallback_dup_multiplier: float,
) -> float:
    """Return a per-event count multiplier using ploidy-aware copy steps when possible."""
    if baseline_ploidy <= 0:
        return fallback_del_multiplier if svtype == "DEL" else fallback_dup_multiplier
    if svtype == "DEL":
        return max(0.0, float(baseline_ploidy - 1) / float(baseline_ploidy))
    if svtype == "DUP":
        return float(baseline_ploidy + 1) / float(baseline_ploidy)
    return 1.0


# ---------------------------------------------------------------------------
# Sample ↔ GD assignment
# ---------------------------------------------------------------------------

def assign_gd_to_samples(
    sample_ids: List[str],
    eligible_entries: List[Tuple[str, dict]],
    rng: np.random.Generator,
    gd_probability: float,
) -> Dict[str, Tuple[str, dict]]:
    """Randomly assign at most one GD entry to each sample.

    Args:
        sample_ids: All sample column names.
        eligible_entries: List of ``(chrom, gd_entry_dict)`` pairs eligible
            for spiking.
        rng: Numpy random generator (for reproducibility).
        gd_probability: Probability that a given sample receives a GD.

    Returns:
        Dict mapping ``sample_id`` → ``(chrom, gd_entry)`` for affected
        samples only.
    """
    if not eligible_entries:
        return {}

    assignments: Dict[str, Tuple[str, dict]] = {}
    n_entries = len(eligible_entries)
    for sid in sample_ids:
        if rng.random() < gd_probability:
            idx = rng.integers(0, n_entries)
            assignments[sid] = eligible_entries[idx]
    return assignments


def _select_sample_subset(
    sample_ids: List[str],
    fraction: float,
    rng: np.random.Generator,
    minimum: int = 0,
) -> List[str]:
    """Select an exact-size random sample subset from a target fraction."""
    if not sample_ids or fraction <= 0:
        return []
    target_n = max(minimum, int(round(len(sample_ids) * fraction)))
    target_n = min(len(sample_ids), target_n)
    if target_n <= 0:
        return []
    chosen = rng.choice(sample_ids, size=target_n, replace=False)
    return [str(sample_id) for sample_id in chosen.tolist()]


def _make_synth_event(
    sample_id: str,
    chrom: str,
    start: int,
    end: int,
    svtype: str,
    multiplier: float,
    baseline_ploidy: int,
    event_id: str,
    source: str,
    cluster: Optional[str] = None,
    gd_id: Optional[str] = None,
    extra: Optional[Dict[str, object]] = None,
) -> dict:
    """Build a normalized synthetic event record used for matrix rewriting."""
    record = {
        "sample_id": str(sample_id),
        "chrom": str(chrom),
        "start": int(start),
        "end": int(end),
        "svtype": str(svtype),
        "multiplier": float(multiplier),
        "baseline_ploidy": int(baseline_ploidy),
        "event_id": str(event_id),
        "source": str(source),
        "cluster": cluster if cluster is not None else "",
        "GD_ID": gd_id if gd_id is not None else "",
    }
    if extra:
        record.update(extra)
    return record


def _build_assignment_events(
    assignments: Dict[str, Tuple[str, dict]],
    del_multiplier: float,
    dup_multiplier: float,
    ploidy_lookup: Optional[Dict[Tuple[str, str], int]] = None,
) -> List[dict]:
    """Convert primary GD assignments to generic synth event records."""
    events: List[dict] = []
    for sample_id, (chrom, entry) in assignments.items():
        svtype = str(entry["svtype"])
        baseline_ploidy = _resolve_sample_contig_ploidy(sample_id, chrom, ploidy_lookup)
        events.append(
            _make_synth_event(
                sample_id=sample_id,
                chrom=chrom,
                start=int(entry["start_GRCh38"]),
                end=int(entry["end_GRCh38"]),
                svtype=svtype,
                multiplier=_resolve_event_multiplier(
                    svtype,
                    baseline_ploidy,
                    del_multiplier,
                    dup_multiplier,
                ),
                baseline_ploidy=baseline_ploidy,
                event_id=str(entry["GD_ID"]),
                source="gd",
                cluster=str(entry.get("cluster", "")),
                gd_id=str(entry["GD_ID"]),
            )
        )
    return events


def _intervals_overlap(start1: int, end1: int, start2: int, end2: int) -> bool:
    """Return True when two half-open genomic intervals overlap."""
    return int(start1) < int(end2) and int(start2) < int(end1)


def _background_conflicts_with_primary(
    sample_id: str,
    chrom: str,
    start: int,
    end: int,
    assignments: Optional[Dict[str, Tuple[str, dict]]],
) -> bool:
    """Return True if a background event would overlap the sample's primary truth event."""
    if not assignments or sample_id not in assignments:
        return False

    primary_chrom, primary_entry = assignments[sample_id]
    if str(primary_chrom) != str(chrom):
        return False

    return _intervals_overlap(
        start,
        end,
        int(primary_entry["start_GRCh38"]),
        int(primary_entry["end_GRCh38"]),
    )


def _matches_canonical_gd_interval(
    locus: GDLocus,
    start: int,
    end: int,
) -> bool:
    """Return True if an interval exactly matches a labeled GD body in the locus."""
    for entry in locus.gd_entries:
        if int(entry["start_GRCh38"]) == int(start) and int(entry["end_GRCh38"]) == int(end):
            return True
    return False


def _matches_canonical_breakpoint_pair(
    locus: GDLocus,
    bp1_name: str,
    bp2_name: str,
) -> bool:
    """Return True if a breakpoint pair already defines a labeled GD entry."""
    left = str(bp1_name)
    right = str(bp2_name)
    for entry in locus.gd_entries:
        if str(entry.get("BP1", "")) == left and str(entry.get("BP2", "")) == right:
            return True
    return False


def generate_salted_flank_bleed_events(
    sample_ids: List[str],
    eligible_loci: List[GDLocus],
    rng: np.random.Generator,
    probability: float,
    del_multiplier: float,
    dup_multiplier: float,
    ploidy_lookup: Optional[Dict[Tuple[str, str], int]] = None,
    assignments: Optional[Dict[str, Tuple[str, dict]]] = None,
) -> List[dict]:
    """Generate nuisance DEL/DUP events that bleed from GD bodies into flanks."""
    selected_samples = _select_sample_subset(sample_ids, probability, rng)
    if not selected_samples or not eligible_loci:
        return []

    events: List[dict] = []
    flank_modes = ["left", "right", "both"]

    for event_idx, sample_id in enumerate(selected_samples, start=1):
        candidate_event = None
        max_attempts = max(1, len(eligible_loci) * 3)
        for _ in range(max_attempts):
            locus = eligible_loci[int(rng.integers(0, len(eligible_loci)))]
            breakpoint_names = list(locus.breakpoint_names)
            if len(breakpoint_names) < 2:
                continue

            left_idx = int(rng.integers(0, len(breakpoint_names) - 1))
            right_idx = int(rng.integers(left_idx + 1, len(breakpoint_names)))
            covered_intervals = locus.get_intervals_between(
                breakpoint_names[left_idx],
                breakpoint_names[right_idx],
            )
            if not covered_intervals:
                covered_intervals = locus.get_intervals()
            if not covered_intervals:
                continue
            if _matches_canonical_breakpoint_pair(
                locus,
                breakpoint_names[left_idx],
                breakpoint_names[right_idx],
            ):
                continue

            body_start = int(covered_intervals[0][0])
            body_end = int(covered_intervals[-1][1])
            if _matches_canonical_gd_interval(locus, body_start, body_end):
                continue
            body_size = max(1, body_end - body_start)
            flank_mode = flank_modes[int(rng.integers(0, len(flank_modes)))]

            left_bleed = 0
            right_bleed = 0
            if flank_mode in {"left", "both"}:
                left_bleed = int(round(body_size * rng.uniform(0.5, 1.0)))
            if flank_mode in {"right", "both"}:
                right_bleed = int(round(body_size * rng.uniform(0.5, 1.0)))

            event_start = max(0, body_start - left_bleed)
            event_end = body_end + right_bleed
            if _background_conflicts_with_primary(
                sample_id,
                locus.chrom,
                event_start,
                event_end,
                assignments,
            ):
                continue

            svtype = "DEL" if rng.random() < 0.5 else "DUP"
            baseline_ploidy = _resolve_sample_contig_ploidy(
                sample_id,
                locus.chrom,
                ploidy_lookup,
            )
            candidate_event = _make_synth_event(
                sample_id=sample_id,
                chrom=locus.chrom,
                start=event_start,
                end=event_end,
                svtype=svtype,
                multiplier=_resolve_event_multiplier(
                    svtype,
                    baseline_ploidy,
                    del_multiplier,
                    dup_multiplier,
                ),
                baseline_ploidy=baseline_ploidy,
                event_id=f"salt_{svtype.lower()}_{event_idx:04d}",
                source="salt",
                cluster=locus.cluster,
                extra={
                    "body_start": body_start,
                    "body_end": body_end,
                    "flank_mode": flank_mode,
                    "body_bp1": breakpoint_names[left_idx],
                    "body_bp2": breakpoint_names[right_idx],
                },
            )
            break

        if candidate_event is not None:
            events.append(candidate_event)

    return events


def _infer_chrom_name(example_chroms: List[str], base_name: str) -> str:
    """Infer chromosome naming style from existing contig names."""
    if any(str(chrom).startswith("chr") for chrom in example_chroms):
        return f"chr{base_name}"
    return base_name


def generate_viable_trisomy_events(
    sample_ids: List[str],
    chrom_examples: List[str],
    rng: np.random.Generator,
    probability: float,
    dup_multiplier: float,
    ploidy_lookup: Optional[Dict[Tuple[str, str], int]] = None,
    assignments: Optional[Dict[str, Tuple[str, dict]]] = None,
) -> List[dict]:
    """Generate whole-chromosome viable trisomy / YY spike-ins."""
    viable_types = [
        ("21", "TRISOMY_21"),
        ("18", "TRISOMY_18"),
        ("13", "TRISOMY_13"),
        ("X", "TRISOMY_X"),
        ("Y", "YY"),
    ]
    minimum = min(len(sample_ids), len(viable_types))
    selected_samples = _select_sample_subset(
        sample_ids,
        probability,
        rng,
        minimum=minimum,
    )
    if not selected_samples:
        return []

    events: List[dict] = []
    chrom_end = 2_000_000_000
    viable_by_sample: List[Tuple[str, Tuple[str, str]]] = []
    fallback_cycle = list(viable_types)
    while len(fallback_cycle) < len(selected_samples):
        fallback_cycle.extend(viable_types)
    rng.shuffle(fallback_cycle)

    for sample_idx, sample_id in enumerate(selected_samples):
        blocked_chrom = None
        if assignments and sample_id in assignments:
            blocked_chrom = str(assignments[sample_id][0])

        candidate_types = [
            item for item in viable_types
            if _infer_chrom_name(chrom_examples, item[0]) != blocked_chrom
        ]
        if not candidate_types:
            candidate_types = viable_types
        chosen_type = candidate_types[int(rng.integers(0, len(candidate_types)))]
        viable_by_sample.append((sample_id, chosen_type))

    for event_idx, (sample_id, (base_chrom, label)) in enumerate(
        viable_by_sample,
        start=1,
    ):
        chrom = _infer_chrom_name(chrom_examples, base_chrom)
        baseline_ploidy = _resolve_sample_contig_ploidy(sample_id, chrom, ploidy_lookup)
        events.append(
            _make_synth_event(
                sample_id=sample_id,
                chrom=chrom,
                start=0,
                end=chrom_end,
                svtype="DUP",
                multiplier=_resolve_event_multiplier(
                    "DUP",
                    baseline_ploidy,
                    dup_multiplier,
                    dup_multiplier,
                ),
                baseline_ploidy=baseline_ploidy,
                event_id=f"aneuploid_{label.lower()}_{event_idx:04d}",
                source="aneuploidy",
                cluster=label,
            )
        )
    return events


def _build_interval_map_from_events(
    events: List[dict],
) -> Dict[Tuple[str, int, int], List[Tuple[str, float]]]:
    """Build a mapping from genomic interval to sample-specific multipliers."""
    interval_map: Dict[Tuple[str, int, int], List[Tuple[str, float]]] = {}
    for event in events:
        key = (str(event["chrom"]), int(event["start"]), int(event["end"]))
        interval_map.setdefault(key, []).append(
            (str(event["sample_id"]), float(event["multiplier"]))
        )
    return interval_map


# ---------------------------------------------------------------------------
# Build a fast lookup: for each (chrom, GD_ID) → set of samples + multiplier
# ---------------------------------------------------------------------------

def _build_spike_index(
    assignments: Dict[str, Tuple[str, dict]],
) -> Dict[Tuple[str, int, int], List[Tuple[int, float]]]:
    """Build a mapping from genomic interval → [(sample_col_index, multiplier)].

    The returned dict is keyed by ``(chrom, start, end)`` of the GD entry.
    Values are lists of ``(column_index_in_sample_array, multiplier)``.

    Column indices are set to -1 here as a placeholder; the caller fills
    them in once the header is known.
    """
    # Group by (chrom, start, end)
    interval_map: Dict[Tuple[str, int, int], List[Tuple[str, float]]] = {}
    for sid, (chrom, entry) in assignments.items():
        key = (chrom, int(entry["start_GRCh38"]), int(entry["end_GRCh38"]))
        mult = 0.5 if entry["svtype"] == "DEL" else 1.5
        interval_map.setdefault(key, []).append((sid, mult))
    return interval_map


# ---------------------------------------------------------------------------
# BGZF / tabix helpers
# ---------------------------------------------------------------------------

# Standard 28-byte empty BGZF block that marks end-of-file (SAM spec §4.1.2).
_BGZF_EOF = bytes([
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00,
    0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00,
    0x00, 0x00, 0x00, 0x00,
])


def _ensure_tabix_index(path: str) -> None:
    """Create a ``.tbi`` tabix index alongside *path* if one is missing.

    The file must be block-gzipped (BGZF) with BED-like columns
    ``chrom, start, end, …`` and a header line beginning with ``#``.
    """
    tbi = path + ".tbi"
    if not os.path.exists(tbi):
        print(f"  Creating tabix index for {os.path.basename(path)} …")
        pysam.tabix_index(
            path, seq_col=0, start_col=1, end_col=2,
            meta_char="#", zerobased=True, force=True,
        )


def _ensure_baf_tabix_index(path: str) -> None:
    """Create a tabix index for a point-based BAF table if needed."""
    tbi = path + ".tbi"
    if os.path.exists(tbi):
        return

    header, _ = _detect_baf_columns(path)
    print(f"  Creating tabix index for {os.path.basename(path)} …")
    pysam.tabix_index(
        path,
        seq_col=0,
        start_col=1,
        end_col=1,
        meta_char="#",
        zerobased=False,
        line_skip=1 if header == 0 else 0,
        force=True,
    )


def _resolve_intervals(
    interval_map: Dict[Tuple[str, int, int], List[Tuple[str, float]]],
    col_index: Dict[str, int],
) -> Dict[str, List[Tuple[int, int, List[Tuple[int, float]]]]]:
    """Map sample-name intervals to column-index intervals.

    Returns a dict keyed by chromosome whose values are lists of
    ``(gd_start, gd_end, [(col_idx, multiplier), …])`` sorted by
    start position.
    """
    intervals_by_chrom: Dict[
        str, List[Tuple[int, int, List[Tuple[str, float]]]]
    ] = {}
    for (chrom, start, end), sample_mults in interval_map.items():
        intervals_by_chrom.setdefault(chrom, []).append(
            (start, end, sample_mults)
        )

    resolved: Dict[str, List[Tuple[int, int, List[Tuple[int, float]]]]] = {}
    for chrom, intervals in intervals_by_chrom.items():
        resolved_intervals = []
        for start, end, sample_mults in intervals:
            resolved_mults = [
                (col_index[sid], mult)
                for sid, mult in sample_mults
                if sid in col_index
            ]
            if resolved_mults:
                resolved_intervals.append((start, end, resolved_mults))
        if resolved_intervals:
            resolved_intervals.sort(key=lambda x: x[0])
            resolved[chrom] = resolved_intervals
    return resolved


# ---------------------------------------------------------------------------
# Per-contig worker (runs in a child process)
# ---------------------------------------------------------------------------


def _process_contig_group(
    input_path: str,
    contig_temp_pairs: List[Tuple[str, str]],
    resolved: Dict[str, List[Tuple[int, int, List[Tuple[int, float]]]]],
    start_col: int,
    end_col: int,
) -> Tuple[int, int]:
    """Fetch assigned contigs via tabix, spike GD signal, write per-contig
    BGZF temp files.

    Each worker opens the block-gzipped input independently (parallel
    decompression) and writes each contig to its own BGZF temp file
    (parallel compression).  Because every contig appears in exactly one
    worker, bins are mutually exclusive across workers.

    Args:
        input_path: BGZF counts file with a ``.tbi`` index.
        contig_temp_pairs: ``[(contig, temp_path), …]`` — contigs this
            worker is responsible for and the per-contig output paths.
        resolved: GD interval lookup — ``chrom →
            [(gd_start, gd_end, [(col_idx, mult), …]), …]``.
        start_col: Column index of the bin-start field.
        end_col: Column index of the bin-end field.

    Returns:
        ``(total_rows, total_modified)``.
    """
    tbx = pysam.TabixFile(input_path)
    total_rows = 0
    total_modified = 0

    for contig, temp_path in contig_temp_pairs:
        chrom_intervals = resolved.get(contig)
        with pysam.BGZFile(temp_path, "w") as fout:
            try:
                for row_str in tbx.fetch(contig):
                    total_rows += 1
                    if chrom_intervals is None:
                        fout.write((row_str + "\n").encode())
                        continue

                    fields = row_str.split("\t")
                    bin_start = int(fields[start_col])
                    bin_end = int(fields[end_col])

                    row_modified = False
                    for gd_start, gd_end, col_mults in chrom_intervals:
                        if bin_end <= gd_start or bin_start >= gd_end:
                            continue
                        for col_idx, mult in col_mults:
                            original = float(fields[col_idx])
                            fields[col_idx] = str(int(round(original * mult)))
                        row_modified = True

                    if row_modified:
                        total_modified += 1
                        fout.write(("\t".join(fields) + "\n").encode())
                    else:
                        fout.write((row_str + "\n").encode())
            except ValueError:
                pass  # contig not present in the file

    tbx.close()
    return total_rows, total_modified


# ---------------------------------------------------------------------------
# BGZF concatenation
# ---------------------------------------------------------------------------


def _concat_bgzf_parts(
    output_path: str,
    part_paths: List[str],
    label: str = "",
) -> None:
    """Concatenate BGZF part files into a single valid BGZF file.

    Intermediate parts have their 28-byte EOF marker stripped; the final
    part keeps its EOF so the result is a well-formed BGZF file.
    """
    n_parts = len(part_paths)
    with open(output_path, "wb") as fout:
        for i, path in enumerate(
            tqdm(part_paths, desc=f"  [{label}] Concatenating",
                 unit=" parts", dynamic_ncols=True)
        ):
            fsize = os.path.getsize(path)
            is_last = (i == n_parts - 1)
            with open(path, "rb") as fin:
                if is_last:
                    shutil.copyfileobj(fin, fout)
                else:
                    # Strip the 28-byte BGZF EOF marker.
                    to_copy = fsize - len(_BGZF_EOF)
                    if to_copy <= 0:
                        continue  # empty part (only EOF)
                    copied = 0
                    while copied < to_copy:
                        chunk = fin.read(min(to_copy - copied, 1 << 20))
                        if not chunk:
                            break
                        fout.write(chunk)
                        copied += len(chunk)


# ---------------------------------------------------------------------------
# Parallel tabix-based TSV rewriter
# ---------------------------------------------------------------------------


def _rewrite_counts_file(
    input_path: str,
    output_path: str,
    interval_map: Dict[Tuple[str, int, int], List[Tuple[str, float]]],
    label: str = "",
    n_workers: int = 1,
) -> int:
    """Rewrite a block-gzipped counts file, spiking GD signal in parallel.

    Uses **tabix** for random-access chunking by contig — each worker
    independently decompresses only its assigned contigs (parallel
    decompression) and writes per-contig BGZF temp files (parallel
    compression).  The temp parts are then binary-concatenated in the
    original contig order and tabix-indexed.

    Memory usage is proportional to a single BGZF block per worker, not
    the full file size.

    Args:
        input_path: Block-gzipped input (``.rd.txt.gz``).
        output_path: Block-gzipped output path.
        interval_map: ``(chrom, start, end) → [(sample_id, mult), …]``.
        label: Human-readable label for progress bars.
        n_workers: Parallel worker processes.

    Returns:
        Number of data rows that were modified.
    """
    # ── Ensure tabix index exists ────────────────────────────────────
    _ensure_tabix_index(input_path)

    file_size = os.path.getsize(input_path)
    print(f"  [{label}] {input_path} "
          f"({file_size / 1e9:.2f} GB compressed)")

    tbx = pysam.TabixFile(input_path)
    contigs = list(tbx.contigs)
    header_lines = list(tbx.header)
    tbx.close()

    if not header_lines:
        raise ValueError(f"No header found in {input_path}")

    # ── Parse header → column index map ──────────────────────────────
    header_text = header_lines[-1]  # last header line has column names
    header_cols = header_text.split("\t")
    if header_cols[0].startswith("#"):
        header_cols[0] = header_cols[0].lstrip("#")
    col_index: Dict[str, int] = {
        name: i for i, name in enumerate(header_cols)
    }
    start_col = col_index.get("Start", 1)
    end_col = col_index.get("End", 2)

    resolved = _resolve_intervals(interval_map, col_index)

    # ── Partition contigs across workers (round-robin for balance) ───
    effective_workers = max(1, min(n_workers, len(contigs)))

    with tempfile.TemporaryDirectory() as tmpdir:
        # Per-contig temp paths, kept in the original contig order so
        # concatenation preserves the sorted genome order.
        contig_temp_paths: List[Tuple[str, str]] = [
            (ctg, os.path.join(tmpdir, f"contig_{i:04d}_{ctg}.bgzf"))
            for i, ctg in enumerate(contigs)
        ]

        # Round-robin assignment — balances load across workers while
        # keeping per-worker groups small enough to avoid memory spikes.
        groups: List[List[Tuple[str, str]]] = [
            [] for _ in range(effective_workers)
        ]
        for i, ct_pair in enumerate(contig_temp_paths):
            groups[i % effective_workers].append(ct_pair)
        groups = [g for g in groups if g]
        n_groups = len(groups)

        print(f"  [{label}] {len(contigs)} contigs → {n_groups} worker "
              f"groups ({effective_workers} workers)")

        # ── Write header as a BGZF part ──────────────────────────────
        header_path = os.path.join(tmpdir, "header.bgzf")
        with pysam.BGZFile(header_path, "w") as hf:
            for h in header_lines:
                hf.write((h + "\n").encode())

        # ── Process contigs in parallel ──────────────────────────────
        total_rows = 0
        total_modified = 0

        if effective_workers <= 1:
            for group in tqdm(
                groups, desc=f"  [{label}] Processing",
                unit=" groups", dynamic_ncols=True,
            ):
                n_rows, n_mod = _process_contig_group(
                    input_path, group, resolved, start_col, end_col,
                )
                total_rows += n_rows
                total_modified += n_mod
        else:
            with ProcessPoolExecutor(
                max_workers=effective_workers,
            ) as pool:
                futures = {
                    pool.submit(
                        _process_contig_group,
                        input_path, group, resolved,
                        start_col, end_col,
                    ): i
                    for i, group in enumerate(groups)
                }
                with tqdm(
                    total=n_groups,
                    desc=f"  [{label}] Processing",
                    unit=" groups",
                    dynamic_ncols=True,
                ) as pbar:
                    for future in as_completed(futures):
                        n_rows, n_mod = future.result()
                        total_rows += n_rows
                        total_modified += n_mod
                        pbar.update(1)

        # ── Concatenate BGZF parts in contig order ───────────────────
        ordered_parts = [header_path] + [tp for _, tp in contig_temp_paths]
        _concat_bgzf_parts(output_path, ordered_parts, label=label)

    # ── Create tabix index on the output ─────────────────────────────
    print(f"  [{label}] Indexing {os.path.basename(output_path)} …")
    pysam.tabix_index(
        output_path, seq_col=0, start_col=1, end_col=2,
        meta_char="#", zerobased=True, force=True,
    )

    print(f"  [{label}] Done: {total_rows:,} rows, {total_modified:,} modified")
    return total_modified


def _spike_baf_value(original_baf: float, svtype: str, baseline_ploidy: int) -> float:
    """Return a synthetic BAF value consistent with the spiked CN state."""
    if svtype == "DEL":
        return 1.0 if original_baf >= 0.5 else 0.0
    if svtype == "DUP":
        if baseline_ploidy <= 1:
            return 1.0 if original_baf >= 0.5 else 0.0
        return (2.0 / 3.0) if original_baf >= 0.5 else (1.0 / 3.0)
    return original_baf


def _build_baf_interval_map(
    events: List[dict],
) -> Dict[str, List[Tuple[int, int, Dict[str, Tuple[str, int]]]]]:
    """Build a chromosome-indexed lookup for BAF spike-ins."""
    by_chrom: Dict[str, Dict[Tuple[int, int], Dict[str, Tuple[str, int]]]] = {}
    for event in events:
        interval_key = (int(event["start"]), int(event["end"]))
        chrom = str(event["chrom"])
        sample_id = str(event["sample_id"])
        chrom_map = by_chrom.setdefault(chrom, {})
        chrom_map.setdefault(interval_key, {})[sample_id] = (
            str(event["svtype"]),
            int(event.get("baseline_ploidy", 2)),
        )

    resolved: Dict[str, List[Tuple[int, int, Dict[str, Tuple[str, int]]]]] = {}
    for chrom, interval_map in by_chrom.items():
        resolved[chrom] = [
            (start, end, sample_map)
            for (start, end), sample_map in sorted(interval_map.items())
        ]
    return resolved


def _process_baf_contig_group(
    input_path: str,
    contig_temp_pairs: List[Tuple[str, str]],
    resolved: Dict[str, List[Tuple[int, int, Dict[str, str]]]],
    pos_col: int,
    baf_col: int,
    sample_col: int,
) -> Tuple[int, int]:
    """Fetch point BAF records by contig, spike values, write BGZF parts."""
    tbx = pysam.TabixFile(input_path)
    total_rows = 0
    total_modified = 0

    for contig, temp_path in contig_temp_pairs:
        chrom_intervals = resolved.get(contig)
        with pysam.BGZFile(temp_path, "w") as fout:
            try:
                for row_str in tbx.fetch(contig):
                    total_rows += 1
                    if chrom_intervals is None:
                        fout.write((row_str + "\n").encode())
                        continue

                    fields = row_str.split("\t")
                    if len(fields) <= max(pos_col, baf_col, sample_col):
                        fout.write((row_str + "\n").encode())
                        continue

                    try:
                        pos = int(fields[pos_col])
                        baf = float(fields[baf_col])
                    except ValueError:
                        fout.write((row_str + "\n").encode())
                        continue

                    sample_id = str(fields[sample_col]).strip()
                    row_modified = False
                    for gd_start, gd_end, sample_map in chrom_intervals:
                        if pos < gd_start or pos >= gd_end:
                            continue
                        event_signature = sample_map.get(sample_id)
                        if event_signature is None:
                            continue
                        svtype, baseline_ploidy = event_signature
                        fields[baf_col] = (
                            f"{_spike_baf_value(baf, svtype, baseline_ploidy):.6f}"
                        )
                        row_modified = True
                        break

                    if row_modified:
                        total_modified += 1
                        fout.write(("\t".join(fields) + "\n").encode())
                    else:
                        fout.write((row_str + "\n").encode())
            except ValueError:
                pass

    tbx.close()
    return total_rows, total_modified


def _rewrite_baf_file(
    input_path: str,
    output_path: str,
    events: List[dict],
    label: str = "baf",
    n_workers: int = 1,
) -> int:
    """Rewrite a tabix-indexed BAF table, spiking per-site BAF values."""
    _ensure_baf_tabix_index(input_path)
    header_flag, column_names = _detect_baf_columns(input_path)
    pos_col = column_names.index("Pos")
    baf_col = column_names.index("BAF")
    sample_col = column_names.index("Sample")

    file_size = os.path.getsize(input_path)
    print(f"  [{label}] {input_path} ({file_size / 1e9:.2f} GB compressed)")

    tbx = pysam.TabixFile(input_path)
    contigs = list(tbx.contigs)
    tabix_header_lines = list(tbx.header)
    tbx.close()

    explicit_header_line = None
    with _open_textfile(input_path, "rt") as handle:
        first_line = handle.readline().rstrip("\n")
    if header_flag == 0:
        explicit_header_line = first_line

    resolved = _build_baf_interval_map(events)
    effective_workers = max(1, min(n_workers, len(contigs)))

    with tempfile.TemporaryDirectory() as tmpdir:
        contig_temp_paths: List[Tuple[str, str]] = [
            (ctg, os.path.join(tmpdir, f"contig_{i:04d}_{ctg}.bgzf"))
            for i, ctg in enumerate(contigs)
        ]

        groups: List[List[Tuple[str, str]]] = [[] for _ in range(effective_workers)]
        for i, ct_pair in enumerate(contig_temp_paths):
            groups[i % effective_workers].append(ct_pair)
        groups = [g for g in groups if g]
        n_groups = len(groups)

        print(f"  [{label}] {len(contigs)} contigs → {n_groups} worker groups ({effective_workers} workers)")

        header_path = os.path.join(tmpdir, "header.bgzf")
        with pysam.BGZFile(header_path, "w") as hf:
            for header_line in tabix_header_lines:
                hf.write((header_line + "\n").encode())
            if explicit_header_line:
                hf.write((explicit_header_line + "\n").encode())

        total_rows = 0
        total_modified = 0
        if effective_workers <= 1:
            for group in tqdm(groups, desc=f"  [{label}] Processing", unit=" groups", dynamic_ncols=True):
                n_rows, n_mod = _process_baf_contig_group(
                    input_path, group, resolved, pos_col, baf_col, sample_col,
                )
                total_rows += n_rows
                total_modified += n_mod
        else:
            with ProcessPoolExecutor(max_workers=effective_workers) as pool:
                futures = {
                    pool.submit(
                        _process_baf_contig_group,
                        input_path,
                        group,
                        resolved,
                        pos_col,
                        baf_col,
                        sample_col,
                    ): i
                    for i, group in enumerate(groups)
                }
                with tqdm(total=n_groups, desc=f"  [{label}] Processing", unit=" groups", dynamic_ncols=True) as pbar:
                    for future in as_completed(futures):
                        n_rows, n_mod = future.result()
                        total_rows += n_rows
                        total_modified += n_mod
                        pbar.update(1)

        ordered_parts = [header_path] + [tp for _, tp in contig_temp_paths]
        _concat_bgzf_parts(output_path, ordered_parts, label=label)

    print(f"  [{label}] Indexing {os.path.basename(output_path)} …")
    pysam.tabix_index(
        output_path,
        seq_col=0,
        start_col=1,
        end_col=1,
        meta_char="#",
        zerobased=False,
        line_skip=1 if explicit_header_line else 0,
        force=True,
    )

    print(f"  [{label}] Done: {total_rows:,} rows, {total_modified:,} modified")
    return total_modified


# ---------------------------------------------------------------------------
# Truth table writer
# ---------------------------------------------------------------------------

def _write_truth_table(
    assignments: Dict[str, Tuple[str, dict]],
    output_path: str,
) -> None:
    """Write a two-column truth table (sample_id, GD_ID).

    Args:
        assignments: Mapping from sample_id → (chrom, gd_entry).
        output_path: Destination TSV path.
    """
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample_id", "GD_ID"])
        for sid in sorted(assignments):
            _, entry = assignments[sid]
            writer.writerow([sid, entry["GD_ID"]])
    print(f"  Truth table: {output_path} ({len(assignments)} carriers)")


def _write_background_event_table(
    events: List[dict],
    output_path: str,
) -> None:
    """Write a manifest of nuisance/background events added during synthesis."""
    columns = [
        "sample_id",
        "event_id",
        "source",
        "chrom",
        "start",
        "end",
        "svtype",
        "multiplier",
        "baseline_ploidy",
        "cluster",
        "GD_ID",
        "body_start",
        "body_end",
        "flank_mode",
        "body_bp1",
        "body_bp2",
    ]
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for event in events:
            row = {column: event.get(column, "") for column in columns}
            writer.writerow(row)
    print(f"  Background events: {output_path} ({len(events)} events)")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    """Parse CLI arguments for the synthesize subcommand."""
    parser = argparse.ArgumentParser(
        description=(
            "Spike synthetic GD CNVs into optional read-depth and BAF tables. "
            "Writes modified output tables plus a truth table of spiked-in events."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--lo-res-counts", required=False,
        help="Low-resolution binned read-count table (.rd.txt.gz)",
    )
    parser.add_argument(
        "--hi-res-counts", required=False,
        help="High-resolution binned read-count table (.rd.txt.gz)",
    )
    parser.add_argument(
        "--baf-table", required=False,
        help="Optional bgzipped, tabix-indexed BAF table (.baf.txt.gz)",
    )
    parser.add_argument(
        "--ploidy-table", required=True,
        help="Preprocess ploidy estimates TSV (ploidy_estimates.tsv). "
             "Synthesize applies per-sample/per-contig ploidy-aware "
             "copy-step multipliers and haploid-aware BAF spiking.",
    )
    parser.add_argument(
        "-g", "--gd-table", required=True,
        help="GD locus definition table (TSV)",
    )
    parser.add_argument(
        "-o", "--output-dir", required=True,
        help="Output directory for modified tables and truth table",
    )
    parser.add_argument(
        "--gd-probability", type=float, default=0.5,
        help="Probability that a given sample will receive a GD event",
    )
    parser.add_argument(
        "--salted-event-probability", type=float, default=0.20,
        help="Fraction of samples that receive an additional nuisance DEL/DUP event "
             "that spans full or partial GD body intervals and bleeds into flanks.",
    )
    parser.add_argument(
        "--viable-trisomy-probability", type=float, default=0.10,
        help="Fraction of samples that receive a whole-chromosome viable trisomy/YY spike-in.",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducibility",
    )
    parser.add_argument(
        "--truth-table", required=False,
        help="Optional previously generated synthesize-format truth table. "
             "When provided, existing sample→GD assignments are reused instead "
             "of drawing new random events.",
    )
    parser.add_argument(
        "--region", dest="regions", action="append", default=None,
        metavar="REGION",
        help=(
            "Restrict eligible GD entries to those overlapping this region. "
            "Accepts a chromosome (e.g. chr22) or an interval "
            "(e.g. chr1:3000000-4000000).  May be specified multiple times."
        ),
    )
    parser.add_argument(
        "--del-multiplier", type=float, default=0.5,
        help="Count multiplier for deletions",
    )
    parser.add_argument(
        "--dup-multiplier", type=float, default=1.5,
        help="Count multiplier for duplications",
    )
    parser.add_argument(
        "--threads", type=int,
        default=max(1, multiprocessing.cpu_count() or 1),
        help="Number of parallel worker processes for matrix rewriting "
             "(default: number of CPUs)",
    )

    return parser.parse_args()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Entry-point for the ``synthesize`` subcommand."""
    args = parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if not any([args.lo_res_counts, args.hi_res_counts, args.baf_table]):
        print(
            "ERROR: Provide at least one input table via --lo-res-counts, "
            "--hi-res-counts, or --baf-table.",
            file=sys.stderr,
        )
        sys.exit(1)

    # ── Load GD table ────────────────────────────────────────────────
    print(f"Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"  {len(gd_table.loci)} loci loaded")
    gd_lookup = _build_gd_lookup(gd_table)

    # ── Collect eligible GD entries ──────────────────────────────────
    supported_gd_lookup = {
        gd_id: assignment
        for gd_id, assignment in gd_lookup.items()
        if assignment[1]["NAHR"] == "yes"
    }
    n_non_nahr = len(gd_lookup) - len(supported_gd_lookup)
    if n_non_nahr > 0:
        print(f"  Restricting synthesis to {len(supported_gd_lookup)} NAHR GD entries; "
              f"skipping {n_non_nahr} non-NAHR entries not modeled by preprocess")

    eligible_entries: List[Tuple[str, dict]] = []
    eligible_loci_by_cluster: Dict[str, GDLocus] = {}
    if args.truth_table:
        if args.regions:
            print("  NOTE: --region is ignored when --truth-table is provided")
        for cluster, locus in gd_table.get_all_loci().items():
            if locus.is_nahr:
                eligible_loci_by_cluster[cluster] = locus
    else:
        parsed_regions = None
        if args.regions:
            parsed_regions = [_parse_region(r) for r in args.regions]
            region_strs = []
            for chrom, start, end in parsed_regions:
                region_strs.append(
                    f"{chrom}" if start is None else f"{chrom}:{start:,}-{end:,}"
                )
            print(f"  Region filter: {', '.join(region_strs)}")

        for cluster, locus in gd_table.get_all_loci().items():
            if not locus.is_nahr:
                continue
            locus_has_eligible_entry = False
            for entry in locus.gd_entries:
                if parsed_regions and not _gd_entry_overlaps_regions(
                    entry, locus.chrom, parsed_regions
                ):
                    continue
                eligible_entries.append((locus.chrom, entry))
                locus_has_eligible_entry = True

            if locus_has_eligible_entry:
                eligible_loci_by_cluster[cluster] = locus

        if not eligible_entries:
            print("ERROR: No eligible GD entries after region filtering.", file=sys.stderr)
            sys.exit(1)
        print(f"  {len(eligible_entries)} eligible GD entries "
              f"({sum(1 for _, e in eligible_entries if e['svtype'] == 'DEL')} DEL, "
              f"{sum(1 for _, e in eligible_entries if e['svtype'] == 'DUP')} DUP)")

    # ── Discover sample universe / assignments ───────────────────────
    sample_ids = _discover_sample_ids(
        args.lo_res_counts,
        args.hi_res_counts,
        args.baf_table,
    )
    if sample_ids:
        print(f"\nDiscovered {len(sample_ids)} samples from input tables")
    elif args.truth_table:
        print("\nNo sample universe detected from inputs; proceeding with truth-table carriers only")
    else:
        print("ERROR: Unable to discover sample IDs from input tables", file=sys.stderr)
        sys.exit(1)

    rng = np.random.default_rng(args.seed)

    ploidy_lookup: Optional[Dict[Tuple[str, str], int]] = None
    if args.ploidy_table:
        print(f"\nLoading ploidy table: {args.ploidy_table}")
        ploidy_lookup = _load_ploidy_lookup(args.ploidy_table)
        print(f"  {len(ploidy_lookup)} sample/contig ploidy records")
    else:
        print("\nNo ploidy table provided; assuming diploid (ploidy=2) everywhere")

    if args.truth_table:
        print(f"\nLoading existing truth assignments: {args.truth_table}")
        assignments = _load_truth_assignments(args.truth_table, gd_lookup)
        dropped_non_nahr = sorted(
            sid for sid, (_, entry) in assignments.items()
            if entry["NAHR"] != "yes"
        )
        if dropped_non_nahr:
            print(f"  WARNING: dropping {len(dropped_non_nahr)} non-NAHR carrier assignments "
                  "because preprocess skips those loci")
            assignments = {
                sid: assignment
                for sid, assignment in assignments.items()
                if assignment[1]["NAHR"] == "yes"
            }
        if sample_ids:
            sample_set = set(sample_ids)
            dropped_samples = sorted(set(assignments) - sample_set)
            if dropped_samples:
                print(
                    f"  WARNING: dropping {len(dropped_samples)} truth-table carriers absent from input tables"
                )
                assignments = {
                    sid: assignment
                    for sid, assignment in assignments.items()
                    if sid in sample_set
                }
        print(f"  Reused {len(assignments)} carrier assignments from truth table")
    else:
        assignments = assign_gd_to_samples(
            sample_ids, eligible_entries, rng, args.gd_probability,
        )
        print(f"\nAssigned GD events to {len(assignments)}/{len(sample_ids)} samples")

    # Summary by GD_ID
    gd_counts: Dict[str, int] = {}
    for _, (_, entry) in assignments.items():
        gd_counts[entry["GD_ID"]] = gd_counts.get(entry["GD_ID"], 0) + 1
    for gd_id, count in sorted(gd_counts.items()):
        print(f"  {gd_id}: {count} carrier(s)")

    if not assignments:
        print("WARNING: No samples assigned a GD event (try increasing "
              "--gd-probability or the number of samples).", file=sys.stderr)

    primary_events = _build_assignment_events(
        assignments,
        del_multiplier=args.del_multiplier,
        dup_multiplier=args.dup_multiplier,
        ploidy_lookup=ploidy_lookup,
    )
    eligible_loci = list(eligible_loci_by_cluster.values())
    salted_events = generate_salted_flank_bleed_events(
        sample_ids,
        eligible_loci,
        rng,
        probability=args.salted_event_probability,
        del_multiplier=args.del_multiplier,
        dup_multiplier=args.dup_multiplier,
        ploidy_lookup=ploidy_lookup,
        assignments=assignments,
    )
    viable_trisomy_events = generate_viable_trisomy_events(
        sample_ids,
        [locus.chrom for locus in gd_table.get_all_loci().values()],
        rng,
        probability=args.viable_trisomy_probability,
        dup_multiplier=args.dup_multiplier,
        ploidy_lookup=ploidy_lookup,
        assignments=assignments,
    )
    background_events = salted_events + viable_trisomy_events
    all_events = primary_events + background_events

    if salted_events:
        print(f"  Added {len(salted_events)} salted flank-bleed event(s)")
    if viable_trisomy_events:
        print(f"  Added {len(viable_trisomy_events)} viable trisomy/YY event(s)")

    interval_map = _build_interval_map_from_events(all_events)

    lo_out = None
    if args.lo_res_counts:
        lo_out = os.path.join(args.output_dir, "lo_res_counts.synthesized.rd.txt.gz")
        print(f"\nRewriting low-res counts → {lo_out}")
        _rewrite_counts_file(
            args.lo_res_counts,
            lo_out,
            interval_map,
            label="lo-res",
            n_workers=args.threads,
        )

    hi_out = None
    if args.hi_res_counts:
        hi_out = os.path.join(args.output_dir, "hi_res_counts.synthesized.rd.txt.gz")
        print(f"\nRewriting high-res counts → {hi_out}")
        _rewrite_counts_file(
            args.hi_res_counts,
            hi_out,
            interval_map,
            label="hi-res",
            n_workers=args.threads,
        )

    baf_out = None
    if args.baf_table:
        baf_out = os.path.join(args.output_dir, "all_samples.synthesized.baf.txt.gz")
        print(f"\nRewriting BAF table → {baf_out}")
        _rewrite_baf_file(
            args.baf_table,
            baf_out,
            all_events,
            label="baf",
            n_workers=args.threads,
        )

    # ── Write truth table ────────────────────────────────────────────
    truth_path = os.path.join(args.output_dir, "truth_table.tsv")
    print(f"\nWriting truth table → {truth_path}")
    _write_truth_table(assignments, truth_path)

    background_truth_path = os.path.join(args.output_dir, "background_events.tsv")
    print(f"Writing background event manifest → {background_truth_path}")
    _write_background_event_table(background_events, background_truth_path)

    # ── Summary ──────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("SYNTHESIS COMPLETE")
    print(f"{'=' * 60}")
    if lo_out:
        print(f"  Low-res output:  {lo_out}")
    if hi_out:
        print(f"  High-res output: {hi_out}")
    if baf_out:
        print(f"  BAF output:      {baf_out}")
    print(f"  Truth table:     {truth_path}")
    print(f"  Background:      {background_truth_path}")
    if sample_ids:
        print(f"  Carriers:        {len(assignments)}/{len(sample_ids)} samples")
    else:
        print(f"  Carriers:        {len(assignments)} samples")
    print(f"  Salted events:   {len(salted_events)}")
    print(f"  Trisomy events:  {len(viable_trisomy_events)}")
    print(f"  Seed:            {args.seed}")
    print(f"  Threads:         {args.threads}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
