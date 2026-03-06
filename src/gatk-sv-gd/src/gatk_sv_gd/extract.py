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

    gatk-sv-gd extract \
        --vcf input.vcf.gz \
        --gd-table gd_table.tsv \
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

# ── Overlap helpers ──────────────────────────────────────────────────


def _overlap_bases(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    """Return the number of overlapping bases between two intervals."""
    return max(0, min(a_end, b_end) - max(a_start, b_start))


def _reciprocal_overlap(
    a_start: int, a_end: int, b_start: int, b_end: int
) -> float:
    """Return the reciprocal overlap (minimum of the two fractional overlaps)."""
    ovl = _overlap_bases(a_start, a_end, b_start, b_end)
    if ovl == 0:
        return 0.0
    len_a = a_end - a_start
    len_b = b_end - b_start
    if len_a <= 0 or len_b <= 0:
        return 0.0
    return min(ovl / len_a, ovl / len_b)


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

    # Track best reciprocal overlap to resolve clustered GD loci
    best_ro = 0.0
    best_entry: Optional[dict] = None

    for entry in gd_entries:
        # Only match same SV type
        if entry["svtype"] != var_svtype:
            continue

        gd_start = entry["start_GRCh38"]
        gd_end = entry["end_GRCh38"]

        ro = _reciprocal_overlap(var_start, var_end, gd_start, gd_end)
        frac = _fraction_covered(gd_start, gd_end, var_start, var_end)

        # Check non-NAHR overlap (looser criterion)
        if frac >= non_nahr_overlap:
            gd_ids.add(entry["GD_ID"])

        # Track best reciprocal-overlap match for NAHR assignment
        if ro >= ro_threshold and ro > best_ro:
            best_ro = ro
            best_entry = entry

    # Assign NAHR flags based on the best-matching entry
    if best_entry is not None:
        gd_ids.add(best_entry["GD_ID"])
        gd_start = best_entry["start_GRCh38"]
        gd_end = best_entry["end_GRCh38"]
        coverage = _fraction_covered(gd_start, gd_end, var_start, var_end)
        if coverage >= atypical_coverage:
            is_nahr = True
        else:
            is_atypical = True

    # Non-NAHR: any GD_ID matched on the loose criterion but NOT already
    # flagged as (canonical or atypical) NAHR.
    if gd_ids and not is_nahr and not is_atypical:
        is_non_nahr = True

    return gd_ids, is_nahr, is_atypical, is_non_nahr


# ── Carrier detection ───────────────────────────────────────────────


def _get_carriers(record: pysam.VariantRecord) -> List[str]:
    """Return sample IDs carrying at least one non-reference allele.

    Accesses the raw GT tuple directly for speed — avoids string
    parsing or allele-object creation.
    """
    carriers: List[str] = []
    for sample in record.samples:
        gt = record.samples[sample].get("GT")
        if gt is None:
            continue
        # Any allele > 0 means non-ref; None means missing (treat as no-call)
        if any(a is not None and a > 0 for a in gt):
            carriers.append(sample)
    return carriers


# ── VCF reading / writing helpers ────────────────────────────────────


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
    return str(svtype)


def _sv_end(record: pysam.VariantRecord) -> int:
    """Return the end coordinate of an SV record."""
    # pysam exposes END via record.stop (0-based exclusive)
    return record.stop


def _ensure_info_headers(header: pysam.VariantHeader) -> pysam.VariantHeader:
    """Add GD INFO fields to *header* if they are not already present."""
    existing_ids = set(header.info)
    for line in _INFO_HEADERS:
        # Parse the ID from the header line
        info_id = line.split("ID=")[1].split(",")[0]
        if info_id not in existing_ids:
            header.add_line(line)
    return header


# ── Core extraction routine ──────────────────────────────────────────


def _build_query_regions(
    gd_table: GDTable,
) -> Dict[str, List[Tuple[int, int, str]]]:
    """Build a dict of chrom → [(start, end, cluster), …] from the GD table.

    Regions are the full span of each GDLocus (min-start to max-end of
    breakpoints), so that any overlapping DEL/DUP will be fetched.
    """
    regions: Dict[str, List[Tuple[int, int, str]]] = {}
    for cluster, locus in gd_table.get_all_loci().items():
        chrom = locus.chrom
        regions.setdefault(chrom, []).append((locus.start, locus.end, cluster))
    # Sort by start within each chrom for deterministic output
    for chrom in regions:
        regions[chrom].sort()
    return regions


def _process_vcfs(
    vcf_paths: List[str],
    gd_table: GDTable,
    ro_threshold: float,
    atypical_coverage: float,
    non_nahr_overlap: float,
) -> Tuple[List[Tuple[pysam.VariantRecord, dict]], Optional[pysam.VariantHeader]]:
    """Query all VCFs and collect annotated records.

    Returns
    -------
    annotated : list of (record, annotation_dict)
        Each annotation_dict has keys ``gd_ids``, ``nahr``, ``atypical``,
        ``non_nahr``.
    out_header : VariantHeader
        Combined header (from first VCF) with GD INFO lines added.
    """
    query_regions = _build_query_regions(gd_table)
    all_loci = gd_table.get_all_loci()

    # We'll collect results keyed by variant identity to deduplicate across
    # VCFs (same variant seen in multiple shards).
    seen_keys: Dict[str, int] = {}  # key → index into annotated list
    annotated: List[Tuple[pysam.VariantRecord, dict]] = []
    out_header: Optional[pysam.VariantHeader] = None

    for vcf_path in vcf_paths:
        vcf = pysam.VariantFile(vcf_path)
        if out_header is None:
            out_header = vcf.header.copy()
            _ensure_info_headers(out_header)

        for chrom, regions in query_regions.items():
            # Validate that the contig exists in this VCF
            if chrom not in vcf.header.contigs:
                continue

            for region_start, region_end, cluster in regions:
                locus = all_loci[cluster]
                try:
                    records = vcf.fetch(chrom, region_start, region_end)
                except ValueError:
                    # Region out of bounds for this VCF
                    continue

                for record in records:
                    svtype = _svtype_of(record)
                    if svtype not in ("DEL", "DUP"):
                        continue

                    var_start = record.start  # 0-based
                    var_end = _sv_end(record)

                    gd_ids, is_nahr, is_atypical, is_non_nahr = _classify_variant(
                        var_start,
                        var_end,
                        svtype,
                        locus.gd_entries,
                        ro_threshold,
                        atypical_coverage,
                        non_nahr_overlap,
                    )

                    if not gd_ids:
                        continue

                    # Deduplicate by (chrom, pos, id, svtype)
                    var_key = f"{record.chrom}:{record.pos}:{record.id}:{svtype}"

                    ann = {
                        "gd_ids": gd_ids,
                        "nahr": is_nahr,
                        "atypical": is_atypical,
                        "non_nahr": is_non_nahr,
                    }

                    if var_key in seen_keys:
                        # Merge annotations (union of GD_IDs, OR of flags)
                        idx = seen_keys[var_key]
                        prev_ann = annotated[idx][1]
                        prev_ann["gd_ids"] |= ann["gd_ids"]
                        prev_ann["nahr"] = prev_ann["nahr"] or ann["nahr"]
                        prev_ann["atypical"] = prev_ann["atypical"] or ann["atypical"]
                        prev_ann["non_nahr"] = prev_ann["non_nahr"] or ann["non_nahr"]
                    else:
                        seen_keys[var_key] = len(annotated)
                        annotated.append((record, ann))

        vcf.close()

    return annotated, out_header


# ── Output writers ───────────────────────────────────────────────────


def _write_vcf(
    annotated: List[Tuple[pysam.VariantRecord, dict]],
    header: pysam.VariantHeader,
    output_path: str,
) -> None:
    """Write annotated records to a VCF file."""
    with pysam.VariantFile(output_path, "w", header=header) as out:
        for record, ann in annotated:
            # Create a new record under the output header
            new_rec = out.new_record(
                contig=record.chrom,
                start=record.start,
                stop=record.stop,
                alleles=record.alleles,
                id=record.id,
                qual=record.qual,
                filter=record.filter,
            )

            # Copy existing INFO fields
            for key in record.info:
                if key in ("GD", "NAHR_GD", "NON_NAHR_GD", "NAHR_GD_atypical"):
                    continue  # will be set fresh below
                try:
                    new_rec.info[key] = record.info[key]
                except (KeyError, TypeError):
                    pass  # header mismatch – skip

            # Set GD annotations
            new_rec.info["GD"] = tuple(sorted(ann["gd_ids"]))
            if ann["nahr"]:
                new_rec.info["NAHR_GD"] = True
            if ann["atypical"]:
                new_rec.info["NAHR_GD_atypical"] = True
            if ann["non_nahr"]:
                new_rec.info["NON_NAHR_GD"] = True

            # Copy FORMAT / sample data
            for sample in record.samples:
                for key in record.samples[sample]:
                    try:
                        new_rec.samples[sample][key] = record.samples[sample][key]
                    except (KeyError, TypeError):
                        pass

            out.write(new_rec)

    # Index the output if bgzipped
    if output_path.endswith(".gz"):
        pysam.tabix_index(output_path, preset="vcf", force=True)


def _write_bed(
    annotated: List[Tuple[pysam.VariantRecord, dict]],
    output_path: str,
) -> None:
    """Write a BED summary with carrier information."""
    with open(output_path, "w") as fh:
        fh.write(
            "#chrom\tstart\tend\tname\tsvtype\tsamples\t"
            "NAHR_GD\tNAHR_GD_atypical\tNON_NAHR_GD\n"
        )
        for record, ann in annotated:
            svtype = _svtype_of(record) or "."
            carriers = _get_carriers(record)
            samples_str = ",".join(carriers) if carriers else "."
            name = record.id if record.id else "."
            nahr_flag = "True" if ann["nahr"] else "False"
            atypical_flag = "True" if ann["atypical"] else "False"
            non_nahr_flag = "True" if ann["non_nahr"] else "False"
            fh.write(
                f"{record.chrom}\t{record.start}\t{record.stop}\t"
                f"{name}\t{svtype}\t{samples_str}\t"
                f"{nahr_flag}\t{atypical_flag}\t{non_nahr_flag}\n"
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
    annotated, out_header = _process_vcfs(
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
    n_nahr = sum(1 for _, a in annotated if a["nahr"])
    n_atypical = sum(1 for _, a in annotated if a["atypical"])
    n_non_nahr = sum(1 for _, a in annotated if a["non_nahr"])
    print(f"    NAHR_GD:          {n_nahr}")
    print(f"    NAHR_GD_atypical: {n_atypical}")
    print(f"    NON_NAHR_GD:      {n_non_nahr}")

    # Write outputs
    vcf_out = os.path.join(args.output_dir, "gd_variants.vcf.gz")
    bed_out = os.path.join(args.output_dir, "gd_variants.bed")

    print(f"Writing VCF to {vcf_out}")
    _write_vcf(annotated, out_header, vcf_out)

    print(f"Writing BED to {bed_out}")
    _write_bed(annotated, bed_out)

    print("Done.")
