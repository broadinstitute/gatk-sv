"""
Extract putative GD events from VCF files.

Reads a GD table to obtain genomic disorder coordinates, then queries
input VCFs (using pysam tabix fetch for efficiency) to pull out all
overlapping DEL and DUP records.  Each record is annotated with INFO
fields describing the type of GD overlap.

Outputs:
  - A VCF containing overlapping records with GD / NAHR_GD /
    NON_NAHR_GD / NAHR_GD_atypical INFO annotations.
  - A BED file summarising variants with carrier sample IDs.

Usage (via CLI)::

    gatk-sv-gd extract \\
        --vcf input.vcf.gz \\
        --gd-table gd_table.tsv \\
        --output-dir results/
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Set, Tuple

import pysam

from gatk_sv_gd.models import GDTable

# ── INFO header definitions ──────────────────────────────────────────

_INFO_HEADERS = [
    '##INFO=<ID=GD,Number=.,Type=String,Description="Genomic disorder region">',
    '##INFO=<ID=NAHR_GD,Number=0,Type=Flag,Description="Variant is a canonical NAHR-mediated genomic disorder CNV">',
    '##INFO=<ID=NON_NAHR_GD,Number=0,Type=Flag,Description="Variant is a non-NAHR-mediated genomic disorder CNV">',
    '##INFO=<ID=NAHR_GD_atypical,Number=0,Type=Flag,Description="Variant is an atypical NAHR-mediated genomic disorder CNV, i.e., it covers 50% of the canonical genomic disorder region but had other breakpoints">',
]

# INFO keys managed by this tool (stripped before re-annotation)
_GD_INFO_KEYS = frozenset({"GD", "NAHR_GD", "NON_NAHR_GD", "NAHR_GD_atypical"})

# ── Overlap helpers ──────────────────────────────────────────────────


def _overlap_bases(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Return the number of overlapping bases between two intervals."""
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _reciprocal_overlap(
    a_start: int, a_end: int, b_start: int, b_end: int
) -> float:
    """Return the reciprocal overlap: overlap / max(len_a, len_b).

    Both intervals must individually cover at least this fraction of the
    *larger* interval for the threshold to be met — equivalent to requiring
    both fractional overlaps to pass simultaneously.
    """
    ovl = _overlap_bases(a_start, a_end, b_start, b_end)
    if ovl == 0:
        return 0.0
    max_len = max(a_end - a_start, b_end - b_start)
    if max_len <= 0:
        return 0.0
    return ovl / max_len


def _fraction_covered(
    region_start: int, region_end: int, query_start: int, query_end: int
) -> float:
    """Fraction of *region* covered by *query*."""
    region_len = region_end - region_start
    if region_len <= 0:
        return 0.0
    return _overlap_bases(region_start, region_end, query_start, query_end) / region_len


# ── Annotation logic ────────────────────────────────────────────────


def _classify_variant(
    var_start: int,
    var_end: int,
    var_svtype: str,
    gd_entries: List[dict],
    ro_threshold: float,
    atypical_coverage: float,
    non_nahr_overlap: float,
) -> Tuple[Set[str], bool, bool, bool]:
    """Classify a variant against GD entries for a single locus.

    Returns
    -------
    gd_ids : set of str
        GD_IDs of matching entries.
    is_nahr : bool
        True if the variant is a canonical NAHR event.
    is_atypical : bool
        True if the variant is an atypical NAHR event.
    is_non_nahr : bool
        True if the variant overlaps sufficiently to be a non-NAHR event.
    """
    gd_ids: Set[str] = set()
    is_nahr = False
    is_atypical = False
    is_non_nahr = False

    # Track best RO match among NAHR='yes' entries only, to resolve
    # clustered GDs with multiple breakpoint combinations.
    best_ro = 0.0
    best_nahr_entry: Optional[dict] = None

    for entry in gd_entries:
        # Only match same SV type
        if entry["svtype"] != var_svtype:
            continue

        gd_start = entry["start_GRCh38"]
        gd_end = entry["end_GRCh38"]
        is_nahr_entry = entry["NAHR"] == "yes"

        # Non-NAHR entries: use the loose fraction-covered threshold.
        if not is_nahr_entry:
            frac = _fraction_covered(gd_start, gd_end, var_start, var_end)
            if frac >= non_nahr_overlap:
                gd_ids.add(entry["GD_ID"])
                is_non_nahr = True

        # NAHR entries: require reciprocal overlap.  The GD_ID is only
        # added below when the RO criterion is met (via best_nahr_entry).
        if is_nahr_entry:
            ro = _reciprocal_overlap(var_start, var_end, gd_start, gd_end)
            if ro >= ro_threshold and ro > best_ro:
                best_ro = ro
                best_nahr_entry = entry

    # Assign canonical / atypical NAHR flags from the best-matching
    # NAHR='yes' entry.
    if best_nahr_entry is not None:
        gd_ids.add(best_nahr_entry["GD_ID"])
        gd_start = best_nahr_entry["start_GRCh38"]
        gd_end = best_nahr_entry["end_GRCh38"]
        coverage = _fraction_covered(gd_start, gd_end, var_start, var_end)
        if coverage >= atypical_coverage:
            is_nahr = True
        else:
            is_atypical = True

    return gd_ids, is_nahr, is_atypical, is_non_nahr


# ── VCF reading helpers ─────────────────────────────────────────────


def _svtype_of(record: pysam.VariantRecord) -> Optional[str]:
    """Extract SVTYPE from a VCF record, or ``None`` if unavailable."""
    svtype = record.info.get("SVTYPE")
    if svtype is None:
        # Some callers encode type in ALT
        for alt in record.alts or []:
            alt_upper = alt.upper().strip("<>")
            if alt_upper in ("DEL", "DUP"):
                return alt_upper
        return None
    # pysam returns a tuple for Number=. INFO fields
    if isinstance(svtype, (tuple, list)):
        svtype = svtype[0] if svtype else None
    return str(svtype) if svtype is not None else None


def _ensure_info_headers(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """Add GD INFO fields to *header* if they are not already present."""
    existing_ids = set(header.info)
    for line in _INFO_HEADERS:
        info_id = line.split("ID=")[1].split(",")[0]
        if info_id not in existing_ids:
            header.add_line(line)
    return header


# ── Query region merging ─────────────────────────────────────────────


def _build_and_merge_query_regions(
    gd_table: GDTable,
) -> Dict[str, List[Tuple[int, int, List[str]]]]:
    """Build per-chrom query regions and merge overlapping ones.

    Returns a dict of ``chrom → [(start, end, [cluster, …]), …]`` where
    overlapping GD loci have been merged into a single fetch region that
    carries all contributing cluster names.  This avoids redundant tabix
    queries and duplicate record processing.
    """
    raw: Dict[str, List[Tuple[int, int, str]]] = {}
    for cluster, locus in gd_table.get_all_loci().items():
        raw.setdefault(locus.chrom, []).append((locus.start, locus.end, cluster))

    merged: Dict[str, List[Tuple[int, int, List[str]]]] = {}
    for chrom, regions in raw.items():
        regions.sort()
        result: List[Tuple[int, int, List[str]]] = []
        for start, end, cluster in regions:
            if result and start <= result[-1][1]:
                prev_start, prev_end, prev_clusters = result[-1]
                result[-1] = (prev_start, max(prev_end, end), prev_clusters + [cluster])
            else:
                result.append((start, end, [cluster]))
        merged[chrom] = result
    return merged


# ── Core extraction routine ──────────────────────────────────────────


def _process_vcfs(
    vcf_paths: List[str],
    gd_table: GDTable,
    ro_threshold: float,
    atypical_coverage: float,
    non_nahr_overlap: float,
) -> Tuple[List[tuple], Optional[pysam.VariantHeader], List[str]]:
    """Query all VCFs and collect annotated variant data.

    Merges overlapping GD query regions so each genomic interval is
    fetched at most once per VCF, and stores lightweight string
    representations of matched records to avoid holding pysam objects
    after the VCF handle is closed.

    Returns
    -------
    annotated : list of (chrom, start, end, var_id, svtype, vcf_line, ann)
        Sorted by contig order then position.  *vcf_line* is the raw
        VCF text line (no trailing newline).  *ann* is a dict with
        keys ``gd_ids`` (set), ``nahr``, ``atypical``, ``non_nahr`` (bool).
    out_header : pysam.VariantHeader
        Combined header (from first VCF) with GD INFO lines added.
    sample_names : list of str
        Sample IDs from the VCF header (used for carrier detection).
    """
    merged_regions = _build_and_merge_query_regions(gd_table)
    all_loci = gd_table.get_all_loci()

    seen_keys: Dict[str, int] = {}  # var_key → index into annotated
    annotated: List[tuple] = []
    out_header: Optional[pysam.VariantHeader] = None
    sample_names: List[str] = []

    for vcf_path in vcf_paths:
        vcf = pysam.VariantFile(vcf_path)
        if out_header is None:
            out_header = vcf.header.copy()
            _ensure_info_headers(out_header)
            sample_names = list(vcf.header.samples)

        for chrom, regions in merged_regions.items():
            if chrom not in vcf.header.contigs:
                continue

            for region_start, region_end, clusters in regions:
                try:
                    fetch_iter = vcf.fetch(chrom, region_start, region_end)
                except ValueError:
                    continue

                for record in fetch_iter:
                    svtype = _svtype_of(record)
                    if svtype not in ("DEL", "DUP"):
                        continue

                    var_start = record.start
                    var_end = record.stop

                    # Classify against every locus in this merged region
                    combined_gd_ids: Set[str] = set()
                    combined_nahr = False
                    combined_atypical = False
                    combined_non_nahr = False

                    for cluster in clusters:
                        locus = all_loci[cluster]
                        gd_ids, nahr, atypical, non_nahr = _classify_variant(
                            var_start, var_end, svtype,
                            locus.gd_entries,
                            ro_threshold, atypical_coverage, non_nahr_overlap,
                        )
                        combined_gd_ids |= gd_ids
                        combined_nahr = combined_nahr or nahr
                        combined_atypical = combined_atypical or atypical
                        combined_non_nahr = combined_non_nahr or non_nahr

                    if not combined_gd_ids:
                        continue

                    var_key = f"{record.chrom}:{record.pos}:{record.id}:{svtype}"
                    ann = {
                        "gd_ids": combined_gd_ids,
                        "nahr": combined_nahr,
                        "atypical": combined_atypical,
                        "non_nahr": combined_non_nahr,
                    }

                    if var_key in seen_keys:
                        idx = seen_keys[var_key]
                        prev_ann = annotated[idx][6]
                        prev_ann["gd_ids"] |= ann["gd_ids"]
                        prev_ann["nahr"] = prev_ann["nahr"] or ann["nahr"]
                        prev_ann["atypical"] = prev_ann["atypical"] or ann["atypical"]
                        prev_ann["non_nahr"] = prev_ann["non_nahr"] or ann["non_nahr"]
                    else:
                        vcf_line = str(record).rstrip("\n\r")
                        seen_keys[var_key] = len(annotated)
                        annotated.append((
                            record.chrom, var_start, var_end,
                            record.id or ".", svtype, vcf_line, ann,
                        ))

        vcf.close()

    # Sort by contig order (matching VCF header) then position
    if out_header is not None:
        contig_order = {c: i for i, c in enumerate(out_header.contigs)}
    else:
        contig_order = {}
    annotated.sort(key=lambda x: (contig_order.get(x[0], 999999), x[1]))

    return annotated, out_header, sample_names


# ── INFO field manipulation ──────────────────────────────────────────


def _modify_info_field(info_str: str, ann: dict) -> str:
    """Return *info_str* with GD annotations replaced / added."""
    parts = []
    for p in info_str.split(";"):
        if not p or p == ".":
            continue
        key = p.split("=", 1)[0]
        if key not in _GD_INFO_KEYS:
            parts.append(p)

    parts.append("GD=" + ",".join(sorted(ann["gd_ids"])))
    if ann["nahr"]:
        parts.append("NAHR_GD")
    if ann["atypical"]:
        parts.append("NAHR_GD_atypical")
    if ann["non_nahr"]:
        parts.append("NON_NAHR_GD")

    return ";".join(parts) if parts else "."


# ── Output writers ───────────────────────────────────────────────────


def _write_vcf(
    annotated: List[tuple],
    header: pysam.VariantHeader,
    output_path: str,
) -> None:
    """Write annotated records to a sorted, bgzipped & indexed VCF.

    Instead of creating new pysam records field-by-field (which is very
    slow for large cohorts), we modify the INFO column of the raw VCF
    text line directly and write the result.
    """
    header_text = str(header)

    # Write uncompressed VCF first, then bgzip + index.
    if output_path.endswith(".gz"):
        tmp_path = output_path[:-3]
    else:
        tmp_path = output_path

    with open(tmp_path, "w") as fh:
        fh.write(header_text)
        for _chrom, _start, _end, _vid, _svtype, vcf_line, ann in annotated:
            fields = vcf_line.split("\t")
            fields[7] = _modify_info_field(fields[7], ann)
            fh.write("\t".join(fields))
            fh.write("\n")

    if output_path.endswith(".gz"):
        pysam.tabix_compress(tmp_path, output_path, force=True)
        os.remove(tmp_path)
        pysam.tabix_index(output_path, preset="vcf", force=True)


def _carriers_from_fields(
    fields: List[str],
    sample_names: List[str],
) -> List[str]:
    """Extract carrier sample IDs from pre-split VCF line fields.

    Parses the GT sub-field directly from the tab-delimited text —
    significantly faster than round-tripping through pysam's per-sample
    dict interface for cohorts with many samples.
    """
    if len(fields) <= 9:
        return []

    # Locate GT in FORMAT (usually the first key)
    fmt_keys = fields[8].split(":")
    try:
        gt_idx = fmt_keys.index("GT")
    except ValueError:
        return []

    carriers: List[str] = []
    n_samples = min(len(fields) - 9, len(sample_names))

    if gt_idx == 0:
        # GT is the first FORMAT field (common case) — fast path:
        # read only up to the first ':' instead of splitting fully.
        for i in range(n_samples):
            sample_field = fields[9 + i]
            colon = sample_field.find(":")
            gt_str = sample_field[:colon] if colon >= 0 else sample_field
            for c in gt_str:
                if "1" <= c <= "9":
                    carriers.append(sample_names[i])
                    break
    else:
        for i in range(n_samples):
            parts = fields[9 + i].split(":", gt_idx + 1)
            if len(parts) <= gt_idx:
                continue
            gt_str = parts[gt_idx]
            for c in gt_str:
                if "1" <= c <= "9":
                    carriers.append(sample_names[i])
                    break

    return carriers


def _write_bed(
    annotated: List[tuple],
    sample_names: List[str],
    output_path: str,
) -> None:
    """Write a BED summary with carrier information."""
    with open(output_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tname\tsvtype\tgd_id\tsamples\tn_samples\t"
            "NAHR_GD\tNAHR_GD_atypical\tNON_NAHR_GD\n"
        )
        for chrom, start, end, var_id, svtype, vcf_line, ann in annotated:
            fields = vcf_line.split("\t")
            carriers = _carriers_from_fields(fields, sample_names)
            samples_str = ",".join(carriers) if carriers else "."
            n_carriers = len(carriers)
            gd_id_str = ",".join(sorted(ann["gd_ids"])) if ann["gd_ids"] else "."
            nahr_flag = "True" if ann["nahr"] else "False"
            atypical_flag = "True" if ann["atypical"] else "False"
            non_nahr_flag = "True" if ann["non_nahr"] else "False"
            fh.write(
                f"{chrom}\t{start}\t{end}\t{var_id}\t{svtype}\t"
                f"{gd_id_str}\t{samples_str}\t{n_carriers}\t{nahr_flag}\t{atypical_flag}\t{non_nahr_flag}\n"
            )


# ── CLI ──────────────────────────────────────────────────────────────


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="gatk-sv-gd extract",
        description=(
            "Extract putative GD events from VCF(s).  Queries GD locus "
            "regions using tabix and annotates overlapping DEL/DUP records."
        ),
    )

    # Input
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "--vcf",
        nargs="+",
        dest="vcfs",
        help="One or more input VCF/BCF files (must be indexed).",
    )
    input_group.add_argument(
        "--vcf-list",
        dest="vcf_list",
        help="Text file listing VCF paths, one per line.",
    )

    parser.add_argument(
        "--gd-table",
        required=True,
        help="Path to GD table TSV file.",
    )

    # Output
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory for output files.",
    )

    # Thresholds
    parser.add_argument(
        "--reciprocal-overlap",
        type=float,
        default=0.5,
        help="Minimum reciprocal overlap for NAHR classification (default: 0.5).",
    )
    parser.add_argument(
        "--atypical-coverage",
        type=float,
        default=0.7,
        help=(
            "Minimum fraction of the GD region that must be covered for a "
            "variant to be canonical NAHR rather than atypical (default: 0.7)."
        ),
    )
    parser.add_argument(
        "--non-nahr-overlap",
        type=float,
        default=0.01,
        help=(
            "Minimum fraction of GD region overlapped for non-NAHR "
            "classification (default: 0.01)."
        ),
    )

    return parser.parse_args()


def main() -> None:
    """Entry point for the *extract* subcommand."""
    args = _parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # Resolve VCF list
    if args.vcf_list:
        with open(args.vcf_list) as fh:
            vcf_paths = [line.strip() for line in fh if line.strip()]
    else:
        vcf_paths = args.vcfs

    # Validate inputs exist
    for path in vcf_paths:
        if not os.path.exists(path):
            print(f"Error: VCF not found: {path}", file=sys.stderr)
            sys.exit(1)
    if not os.path.exists(args.gd_table):
        print(f"Error: GD table not found: {args.gd_table}", file=sys.stderr)
        sys.exit(1)

    # Load GD table
    print(f"Loading GD table from {args.gd_table}")
    gd_table = GDTable(args.gd_table)
    n_loci = len(gd_table.get_all_loci())
    n_entries = len(gd_table.df)
    print(f"  {n_loci} loci, {n_entries} GD entries")

    # Process VCFs
    print(f"Processing {len(vcf_paths)} VCF(s) …")
    annotated, out_header, sample_names = _process_vcfs(
        vcf_paths,
        gd_table,
        ro_threshold=args.reciprocal_overlap,
        atypical_coverage=args.atypical_coverage,
        non_nahr_overlap=args.non_nahr_overlap,
    )

    if not annotated:
        print("No overlapping DEL/DUP records found.")
        sys.exit(0)

    print(f"  {len(annotated)} overlapping records found")
    n_nahr = sum(1 for *_, a in annotated if a["nahr"])
    n_atypical = sum(1 for *_, a in annotated if a["atypical"])
    n_non_nahr = sum(1 for *_, a in annotated if a["non_nahr"])
    print(f"    NAHR_GD:          {n_nahr}")
    print(f"    NAHR_GD_atypical: {n_atypical}")
    print(f"    NON_NAHR_GD:      {n_non_nahr}")

    # Write outputs
    vcf_out = os.path.join(args.output_dir, "gd_variants.vcf.gz")
    bed_out = os.path.join(args.output_dir, "gd_variants.bed")

    print(f"Writing VCF to {vcf_out}")
    _write_vcf(annotated, out_header, vcf_out)

    print(f"Writing BED to {bed_out}")
    _write_bed(annotated, sample_names, bed_out)

    print("Done.")
