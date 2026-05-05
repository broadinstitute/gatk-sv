#!/usr/bin/env python3
"""Analyze VEP annotations and PREDICTED fields by gene.

Counts variants per gene grouped by consequence type (VEP-based and PREDICTED_*).
Uses vep_rank.tsv for consequence prioritization.
"""

import argparse
from collections import defaultdict
from typing import Dict, List, Set, Tuple


def load_vep_rank(vep_rank_path: str) -> List[str]:
    """Load VEP consequence ranking (top to bottom = highest to lowest priority)."""
    rank = []
    with open(vep_rank_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if cols:
                consequence = cols[0].strip()
                if consequence:
                    rank.append(consequence)
    return rank


def load_gene_list(gene_list_path: str) -> Set[str]:
    """Load gene list from file (1st column)."""
    genes = set()
    with open(gene_list_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if cols:
                gene = cols[0].strip()
                if gene:
                    genes.add(gene)
    return genes


def classify_consequence(consequences_str: str, vep_rank: List[str]) -> str:
    """
    Classify a consequence string (delimited by "&" or ",") using VEP ranking.
    Splits by both delimiters, then returns the highest-priority consequence.
    """
    if not consequences_str or consequences_str == ".":
        return "."

    # Split by both "&" and "," delimiters
    consequences = []
    for part in consequences_str.split("&"):
        for subpart in part.split(","):
            c = subpart.strip()
            if c:
                consequences.append(c)

    if not consequences:
        return "."

    # Find the first (highest priority) consequence in the rank
    for ranked_cons in vep_rank:
        for cons in consequences:
            if cons == ranked_cons:
                return cons

    # If no consequence matches the rank, return the first one
    return consequences[0] if consequences else "."


def parse_input_file(
    input_path: str,
    vep_consequence_col: int,
    vep_gene_col: int,
    vep_rank: List[str],
    predicted_cols: Dict[str, int],
) -> Dict[str, Dict[str, int]]:
    """
    Parse input file and build gene x consequence count matrix.
    Returns: {gene: {consequence_or_predicted_type: count}}
    """
    gene_counts: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))

    with open(input_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split("\t")
            # Only require the minimum columns we're going to access
            if len(cols) <= vep_consequence_col or len(cols) <= vep_gene_col:
                continue

            # Get gene
            gene = cols[vep_gene_col].strip()
            if not gene or gene == ".":
                continue

            # Get VEP consequence and classify
            vep_cons_raw = cols[vep_consequence_col].strip() if vep_consequence_col < len(cols) else "."
            vep_cons = classify_consequence(vep_cons_raw, vep_rank)
            if vep_cons and vep_cons != ".":
                gene_counts[gene][vep_cons] += 1

            # Count PREDICTED_* fields
            for pred_type, col_idx in predicted_cols.items():
                if col_idx < len(cols):
                    val = cols[col_idx].strip()
                    # Count if value is not "." or empty
                    if val and val not in (".", "0", ""):
                        gene_counts[gene][pred_type] += 1

    return gene_counts


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Analyze VEP annotations and PREDICTED fields by gene. "
            "Outputs gene x consequence count matrix."
        )
    )
    parser.add_argument("--input", required=True, help="Input VEP-annotated file (tab-delimited)")
    parser.add_argument("--vep-rank", required=True, help="VEP ranking file (consequences in priority order)")
    parser.add_argument("--gene-list", required=True, help="Gene list file (gene names in 1st column)")
    parser.add_argument(
        "--vep-consequence-col",
        type=int,
        default=10,
        help="0-based column index for VEP_Consequence (default: 10)",
    )
    parser.add_argument(
        "--vep-gene-col",
        type=int,
        default=9,
        help="0-based column index for VEP_Gene (default: 9)",
    )
    parser.add_argument(
        "--predicted-cols",
        type=str,
        help=(
            "Comma-separated list of PREDICTED_* column definitions "
            "(e.g., 'PREDICTED_LOF:11,PREDICTED_COPY_GAIN:12')"
        ),
    )
    parser.add_argument("--output", required=True, help="Output TSV file")

    args = parser.parse_args()

    # Load VEP ranking
    vep_rank = load_vep_rank(args.vep_rank)
    print(f"Loaded {len(vep_rank)} VEP consequence types from ranking", flush=True)

    # Load gene list
    gene_list = load_gene_list(args.gene_list)
    print(f"Loaded {len(gene_list)} genes from gene list", flush=True)

    # Parse PREDICTED_* columns
    predicted_cols: Dict[str, int] = {}
    if args.predicted_cols:
        for item in args.predicted_cols.split(","):
            parts = item.strip().split(":")
            if len(parts) == 2:
                pred_type = parts[0].strip()
                col_idx = int(parts[1].strip())
                predicted_cols[pred_type] = col_idx

    print(f"Tracking {len(predicted_cols)} PREDICTED_* fields", flush=True)

    # Parse input file
    gene_counts = parse_input_file(
        args.input,
        args.vep_consequence_col,
        args.vep_gene_col,
        vep_rank,
        predicted_cols,
    )

    print(f"Found {len(gene_counts)} genes with variants", flush=True)

    # Prepare output columns
    vep_consequence_types = set()
    for counts in gene_counts.values():
        vep_consequence_types.update(
            [k for k in counts.keys() if k not in predicted_cols]
        )
    vep_consequence_types = sorted(vep_consequence_types)
    predicted_types = sorted(predicted_cols.keys())

    header = ["gene"] + vep_consequence_types + predicted_types

    # Write output
    with open(args.output, "w") as fout:
        fout.write("\t".join(header) + "\n")

        for gene in sorted(gene_list):
            row = [gene]
            counts = gene_counts.get(gene, {})

            # VEP consequence columns
            for cons_type in vep_consequence_types:
                row.append(str(counts.get(cons_type, 0)))

            # PREDICTED_* columns
            for pred_type in predicted_types:
                row.append(str(counts.get(pred_type, 0)))

            # Only output if at least one data column (from index 1 onwards) is non-zero
            has_variant = any(int(val) > 0 for val in row[1:])
            if has_variant:
                fout.write("\t".join(row) + "\n")

    print({"output": args.output, "genes_written": sum(1 for _ in open(args.output)) - 1}, flush=True)


if __name__ == "__main__":
    main()
