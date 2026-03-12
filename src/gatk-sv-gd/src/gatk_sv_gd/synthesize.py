"""
Synthesize GD CNVs into read-depth matrices.

Spikes simulated genomic-disorder CNV signal into high- and low-resolution
binned read-count matrices.  For each sample selected to carry a GD event
(controlled by ``--gd-probability``), a single randomly chosen GD entry is
assigned and the depth counts in bins overlapping the GD locus are scaled:

* **Deletions** – counts multiplied by 0.5
* **Duplications** – counts multiplied by 1.5

Counts are rounded to the nearest integer after scaling.

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
from typing import Dict, List, Optional, Tuple

import numpy as np
import pysam
from tqdm import tqdm

from gatk_sv_gd.models import GDTable


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
        ordered_parts = (
            [header_path]
            + [tp for _, tp in contig_temp_paths]
        )
        _concat_bgzf_parts(output_path, ordered_parts, label=label)

    # ── Create tabix index on the output ─────────────────────────────
    print(f"  [{label}] Indexing {os.path.basename(output_path)} …")
    pysam.tabix_index(
        output_path, seq_col=0, start_col=1, end_col=2,
        meta_char="#", zerobased=True, force=True,
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


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args():
    """Parse CLI arguments for the synthesize subcommand."""
    parser = argparse.ArgumentParser(
        description=(
            "Spike synthetic GD CNVs into read-depth count matrices. "
            "Writes modified high- and low-resolution depth tables and "
            "a truth table of spiked-in events."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--lo-res-counts", required=True,
        help="Low-resolution binned read-count table (.rd.txt.gz)",
    )
    parser.add_argument(
        "--hi-res-counts", required=True,
        help="High-resolution binned read-count table (.rd.txt.gz)",
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
        "--seed", type=int, default=42,
        help="Random seed for reproducibility",
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

    # ── Load GD table ────────────────────────────────────────────────
    print(f"Loading GD table: {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    print(f"  {len(gd_table.loci)} loci loaded")

    # ── Collect eligible GD entries ──────────────────────────────────
    parsed_regions = None
    if args.regions:
        parsed_regions = [_parse_region(r) for r in args.regions]
        region_strs = []
        for chrom, start, end in parsed_regions:
            region_strs.append(
                f"{chrom}" if start is None else f"{chrom}:{start:,}-{end:,}"
            )
        print(f"  Region filter: {', '.join(region_strs)}")

    eligible_entries: List[Tuple[str, dict]] = []
    for cluster, locus in gd_table.get_all_loci().items():
        for entry in locus.gd_entries:
            if parsed_regions and not _gd_entry_overlaps_regions(
                entry, locus.chrom, parsed_regions
            ):
                continue
            eligible_entries.append((locus.chrom, entry))

    if not eligible_entries:
        print("ERROR: No eligible GD entries after region filtering.", file=sys.stderr)
        sys.exit(1)
    print(f"  {len(eligible_entries)} eligible GD entries "
          f"({sum(1 for _, e in eligible_entries if e['svtype'] == 'DEL')} DEL, "
          f"{sum(1 for _, e in eligible_entries if e['svtype'] == 'DUP')} DUP)")

    # ── Discover sample columns from the lo-res header ───────────────
    print(f"\nReading sample columns from: {args.lo_res_counts}")
    with gzip.open(args.lo_res_counts, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
    meta_cols = {"#Chr", "Chr", "Start", "End", "source_file"}
    sample_ids = [c for c in header if c not in meta_cols]
    print(f"  {len(sample_ids)} samples")

    # ── Assign GDs to samples ────────────────────────────────────────
    rng = np.random.default_rng(args.seed)
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

    # ── Apply custom multipliers ─────────────────────────────────────
    # Patch entry-level multipliers if the user overrode defaults
    if args.del_multiplier != 0.5 or args.dup_multiplier != 1.5:
        patched: Dict[str, Tuple[str, dict]] = {}
        for sid, (chrom, entry) in assignments.items():
            patched[sid] = (chrom, entry)
        assignments = patched

    # Build the interval → [(sample, multiplier)] index
    # We override the multiplier based on svtype + CLI args
    interval_map: Dict[Tuple[str, int, int], List[Tuple[str, float]]] = {}
    for sid, (chrom, entry) in assignments.items():
        key = (chrom, int(entry["start_GRCh38"]), int(entry["end_GRCh38"]))
        mult = (args.del_multiplier if entry["svtype"] == "DEL"
                else args.dup_multiplier)
        interval_map.setdefault(key, []).append((sid, mult))

    # ── Rewrite lo-res counts ────────────────────────────────────────
    lo_out = os.path.join(args.output_dir, "lo_res_counts.synthesized.rd.txt.gz")
    print(f"\nRewriting low-res counts → {lo_out}")
    _rewrite_counts_file(args.lo_res_counts, lo_out, interval_map,
                         label="lo-res", n_workers=args.threads)

    # ── Rewrite hi-res counts ────────────────────────────────────────
    hi_out = os.path.join(args.output_dir, "hi_res_counts.synthesized.rd.txt.gz")
    print(f"\nRewriting high-res counts → {hi_out}")
    _rewrite_counts_file(args.hi_res_counts, hi_out, interval_map,
                         label="hi-res", n_workers=args.threads)

    # ── Write truth table ────────────────────────────────────────────
    truth_path = os.path.join(args.output_dir, "truth_table.tsv")
    print(f"\nWriting truth table → {truth_path}")
    _write_truth_table(assignments, truth_path)

    # ── Summary ──────────────────────────────────────────────────────
    print(f"\n{'=' * 60}")
    print("SYNTHESIS COMPLETE")
    print(f"{'=' * 60}")
    print(f"  Low-res output:  {lo_out}")
    print(f"  High-res output: {hi_out}")
    print(f"  Truth table:     {truth_path}")
    print(f"  Carriers:        {len(assignments)}/{len(sample_ids)} samples")
    print(f"  Seed:            {args.seed}")
    print(f"  Threads:         {args.threads}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
