"""
Shared data models for GD CNV analysis.

Contains the core data structures (GDLocus, GDTable) used across the
infer, call, plot, and eval modules.
"""

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import pandas as pd


@dataclass
class GDLocus:
    """Represents a genomic disorder locus with its breakpoints."""
    cluster: str
    chrom: str
    breakpoints: List[Tuple[int, int]]  # Sorted list of breakpoint ranges (start, end)
    breakpoint_names: List[str]  # Names of breakpoints (e.g., ['1', '2', '3'] or ['A', 'B', 'C'])
    gd_entries: List[dict]  # Original GD table entries for this locus
    is_nahr: bool
    is_terminal: bool

    @property
    def n_breakpoints(self) -> int:
        return len(self.breakpoints)

    @property
    def start(self) -> int:
        """Get overall start position (min of all breakpoint starts)."""
        return min(bp[0] for bp in self.breakpoints)

    @property
    def end(self) -> int:
        """Get overall end position (max of all breakpoint ends)."""
        return max(bp[1] for bp in self.breakpoints)

    def get_intervals(self) -> List[Tuple[int, int, str]]:
        """
        Get all intervals between adjacent breakpoints.
        Intervals span from the end of one breakpoint region to the start of the next.

        Returns:
            List of (start, end, interval_name) tuples
        """
        intervals = []
        for i in range(len(self.breakpoints) - 1):
            # Interval starts at end of BP_i and ends at start of BP_{i+1}
            start = self.breakpoints[i][1]  # End of current breakpoint
            end = self.breakpoints[i + 1][0]  # Start of next breakpoint
            # Use actual breakpoint names from the table
            interval_name = f"{self.breakpoint_names[i]}-{self.breakpoint_names[i + 1]}"
            intervals.append((start, end, interval_name))
        return intervals

    def get_intervals_between(self, bp1_name: str, bp2_name: str) -> List[Tuple[int, int, str]]:
        """
        Get the intervals covered between two named breakpoints.

        Uses the authoritative breakpoint ordering (``breakpoint_names``)
        rather than parsing interval name strings, so this is safe even if
        BP names contain hyphens or other special characters.

        Args:
            bp1_name: Name of the first breakpoint (from GD table ``BP1`` column).
            bp2_name: Name of the second breakpoint (from GD table ``BP2`` column).

        Returns:
            List of (start, end, interval_name) tuples for every interval
            whose left BP is at or after *bp1_name* and whose right BP is
            at or before *bp2_name* in the ordering.  Returns an empty
            list if either name is not found.
        """
        try:
            pos1 = self.breakpoint_names.index(bp1_name)
            pos2 = self.breakpoint_names.index(bp2_name)
        except ValueError:
            return []
        if pos1 > pos2:
            pos1, pos2 = pos2, pos1
        all_intervals = self.get_intervals()
        return [
            (start, end, name)
            for idx, (start, end, name) in enumerate(all_intervals)
            if idx >= pos1 and (idx + 1) <= pos2
        ]

    def get_gd_intervals(self) -> Dict[str, Tuple[int, int]]:
        """
        Get intervals for each GD entry (each DEL/DUP defined in the table).

        Returns:
            Dict mapping GD_ID to (start, end) tuple
        """
        return {
            entry["GD_ID"]: (entry["start_GRCh38"], entry["end_GRCh38"])
            for entry in self.gd_entries
        }

    def get_flanking_regions(self) -> List[Tuple[int, int, str]]:
        """
        Get flanking regions on either side of the locus (100% of locus size each side).

        These regions are used to detect large spanning CNVs that overlap the locus
        incidentally but are not true GD events. A true GD only affects the region
        between its breakpoints; a spanning variant would also show copy number change
        in these flanking regions.

        Note: The left flank is clipped to position 0 if the locus is near the
        chromosome start. Flanks at the chromosome end will simply have no bins.

        Returns:
            List of (start, end, name) tuples. Contains 0, 1, or 2 elements
            depending on chromosome proximity.
        """
        locus_size = self.end - self.start
        flanks = []

        # Left flank: extend left by 100% of locus size
        left_start = max(0, self.start - locus_size)
        left_end = self.start
        if left_end > left_start:
            flanks.append((left_start, left_end, "left_flank"))

        # Right flank: extend right by 100% of locus size
        right_start = self.end
        right_end = self.end + locus_size
        flanks.append((right_start, right_end, "right_flank"))

        return flanks

    @property
    def del_entries(self) -> List[dict]:
        """Get all deletion entries for this locus."""
        return [e for e in self.gd_entries if e["svtype"] == "DEL"]

    @property
    def dup_entries(self) -> List[dict]:
        """Get all duplication entries for this locus."""
        return [e for e in self.gd_entries if e["svtype"] == "DUP"]

    @property
    def svtypes(self) -> List[str]:
        """Get unique svtypes at this locus."""
        return list(set(e["svtype"] for e in self.gd_entries))


class GDTable:
    """Parser and container for GD locus definitions."""

    @staticmethod
    def _standalone_locus_key(row: pd.Series) -> str:
        """Return a stable locus key for GD rows without a cluster label."""
        chrom = str(row["chr"])
        start = int(row["start_GRCh38"])
        end = int(row["end_GRCh38"])
        return f"{chrom}:{start}-{end}"

    def __init__(self, filepath: str):
        """
        Load GD table from TSV file.

        Expected columns:
            chr, start_GRCh38, end_GRCh38, GD_ID, svtype, NAHR, terminal, cluster

        Args:
            filepath: Path to GD table TSV file
        """
        self.df = pd.read_csv(filepath, sep="\t")
        self._validate_columns()
        self.loci = self._parse_loci()

    def _validate_columns(self):
        """Validate required columns exist."""
        required_cols = ["chr", "start_GRCh38", "end_GRCh38", "GD_ID", "svtype",
                         "NAHR", "terminal", "cluster", "BP1", "BP2"]
        missing = set(required_cols) - set(self.df.columns)
        if missing:
            raise ValueError(f"Missing required columns in GD table: {missing}")

    def _parse_loci(self) -> Dict[str, GDLocus]:
        """
        Parse GD table into loci grouped by cluster.

        Rows without a cluster label are grouped by genomic coordinates so
        standalone DEL/DUP definitions for the same interval share one bin set.
        Breakpoints are extracted from BP1 and BP2 columns.

        Returns:
            Dict mapping cluster name to GDLocus object
        """
        loci = {}

        # Group by cluster, only processing entries with a defined cluster
        for _, row in self.df.iterrows():
            cluster = row["cluster"]
            
            # For entries without a cluster (standalone regions), use the
            # genomic interval as the locus key so DEL/DUP rows at the same
            # coordinates share one modeled locus and one bin set.
            if pd.isna(cluster) or cluster == "":
                cluster = self._standalone_locus_key(row)

            if cluster not in loci:
                loci[cluster] = {
                    "chrom": row["chr"],
                    "breakpoint_coords": {},  # Dict mapping BP name to list of coordinates
                    "entries": [],
                    "is_nahr": row["NAHR"] == "yes",
                    "is_terminal": row["terminal"] != "no" if pd.notna(row["terminal"]) else False,
                }

            # Get breakpoint names from BP1 and BP2 columns
            bp1 = row["BP1"]
            bp2 = row["BP2"]
            
            # If BP1 or BP2 is missing, use default values "1" and "2"
            if pd.isna(bp1) or pd.isna(bp2):
                bp1 = "1"
                bp2 = "2"
            else:
                # Convert to strings for consistent handling.
                # Pandas reads numeric columns with NaN gaps as float, so a
                # TSV value of "2" arrives as 2.0.  Normalise to clean integer
                # strings ("2", not "2.0") to avoid duplicate breakpoint names
                # and broken interval matching downstream.
                bp1 = str(int(bp1)) if isinstance(bp1, float) and bp1 == int(bp1) else str(bp1)
                bp2 = str(int(bp2)) if isinstance(bp2, float) and bp2 == int(bp2) else str(bp2)
            
            # Track coordinates for each breakpoint
            if bp1 not in loci[cluster]["breakpoint_coords"]:
                loci[cluster]["breakpoint_coords"][bp1] = []
            if bp2 not in loci[cluster]["breakpoint_coords"]:
                loci[cluster]["breakpoint_coords"][bp2] = []
            
            loci[cluster]["breakpoint_coords"][bp1].append(row["start_GRCh38"])
            loci[cluster]["breakpoint_coords"][bp2].append(row["end_GRCh38"])

            # Add entry
            loci[cluster]["entries"].append({
                "GD_ID": row["GD_ID"],
                "start_GRCh38": row["start_GRCh38"],
                "end_GRCh38": row["end_GRCh38"],
                "svtype": row["svtype"],
                "NAHR": row["NAHR"],
                "terminal": row["terminal"],
                "BP1": bp1,
                "BP2": bp2,
            })

        # Convert to GDLocus objects
        result = {}
        for cluster, data in loci.items():
            # Convert breakpoint coordinates to ranges
            bp_ranges = []
            
            if data["breakpoint_coords"]:
                # Sort breakpoints by their genomic position (min coordinate)
                # so that the ordering reflects physical location on the
                # chromosome.  Name-based sorting fails when numeric and
                # alphanumeric IDs are mixed (e.g. 1, 2, 3, 4, CHRNA7, 5, 6
                # at the 15q11-13 locus — CHRNA7 sits between BP4 and BP5
                # genomically but would sort after BP6 alphabetically).
                bp_names_sorted = sorted(
                    data["breakpoint_coords"].keys(),
                    key=lambda bp_name: min(data["breakpoint_coords"][bp_name]),
                )
                
                for bp_name in bp_names_sorted:
                    coords = data["breakpoint_coords"][bp_name]
                    # Breakpoint range spans from min to max observed coordinate
                    bp_range = (min(coords), max(coords))
                    bp_ranges.append(bp_range)
            else:
                # No breakpoints found at all - this shouldn't happen, but handle it
                print(f"  ERROR: No breakpoints found for {cluster}, skipping")
                continue
            
            result[cluster] = GDLocus(
                cluster=cluster,
                chrom=data["chrom"],
                breakpoints=bp_ranges,
                breakpoint_names=bp_names_sorted,
                gd_entries=data["entries"],
                is_nahr=data["is_nahr"],
                is_terminal=data["is_terminal"],
            )
            
            # Print breakpoint ranges for debugging
            print(f"  {cluster}: {len(bp_ranges)} breakpoints")
            for i, (start, end) in enumerate(bp_ranges):
                width = end - start
                bp_label = bp_names_sorted[i]
                print(f"    BP {bp_label}: {start:,} - {end:,} (width: {width:,} bp)")

        return result

    def get_locus(self, cluster: str) -> Optional[GDLocus]:
        """Get a specific locus by cluster name."""
        return self.loci.get(cluster)

    def get_all_loci(self) -> Dict[str, GDLocus]:
        """Get all loci."""
        return self.loci

    def get_loci_by_chrom(self, chrom: str) -> Dict[str, GDLocus]:
        """Get all loci on a specific chromosome."""
        return {k: v for k, v in self.loci.items() if v.chrom == chrom}


def validate_gd_table_for_preprocess(gd_table: GDTable) -> None:
    """Reject GD loci that preprocess does not support yet."""
    unsupported = [
        f"{cluster} ({locus.chrom}:{locus.start:,}-{locus.end:,})"
        for cluster, locus in gd_table.loci.items()
        if str(locus.chrom).strip().lower() in {"chry", "y"}
    ]
    if unsupported:
        joined = ", ".join(sorted(unsupported))
        raise ValueError(
            "gatk-sv-gd preprocess does not support GD loci on chrY. "
            "We do not model GDs on chrY yet, so please remove these entries "
            f"from the GD table before preprocessing: {joined}"
        )
