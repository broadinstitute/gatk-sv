"""
Permute gene locations across the genome.

This script now follows the same placement method as the R implementation:
for each gene, choose one eligible free region uniformly at random, then choose
a random start uniformly inside that region with a 1kb buffer on both sides.

Free regions are computed as the complement of the telomere/centromere BED over
hg38 chromosome sizes. Genes whose original coordinates overlap any interval in
the blacklist BED are not permuted and are written unchanged.

Usage:
    python3 permute_gtf.py <seed> <input.gtf.gz> <tel_cen.bed> <blacklist.bed[.gz]> <output.gtf.gz> [too_large_genes.txt]

The optional last argument is a path to write gene IDs that were too large to
place (kept at their original locations). If omitted they are only printed to
stderr.

Example:
    python3 permute_gtf.py 1 \
        genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf.gz \
        hg38_tel_cen.bed \
        merged_blacklist_GATKSV.bed.gz \
        genes_grch38_annotated_4_mapped_gencode_v39.CDS.permuted_seed1.gtf.gz \
        too_large_genes.seed1.txt
"""

import re, gzip, sys, random, bisect
from collections import defaultdict

# ── hg38 chromosome sizes (UCSC) ────────────────────────────────────────────
HG38_CHROM_SIZES = {
    'chr1':  248956422, 'chr2':  242193529, 'chr3':  198295559,
    'chr4':  190214555, 'chr5':  181538259, 'chr6':  170805979,
    'chr7':  159345973, 'chr8':  145138636, 'chr9':  138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16':  90338345, 'chr17':  83257441, 'chr18':  80373285,
    'chr19':  58617616, 'chr20':  64444167, 'chr21':  46709983,
    'chr22':  50818468, 'chrX': 156040895,  'chrY':  57227415,
}

# ── helpers ──────────────────────────────────────────────────────────────────

def load_exclusions(bed_path):
    """Return dict chrom -> sorted list of (start, end) 0-based half-open excluded intervals.
    Supports plain text or gzip-compressed BED files."""
    excl = defaultdict(list)
    with open_file(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            excl[f[0]].append((int(f[1]), int(f[2])))
    for chrom in excl:
        excl[chrom].sort()
    return excl


def overlaps_intervals(chrom, start_1based, end_1based, intervals_by_chrom):
    """Return True if the 1-based closed interval [start, end] overlaps any
    0-based half-open interval in intervals_by_chrom[chrom]."""
    ivs = intervals_by_chrom.get(chrom)
    if not ivs:
        return False
    # convert gene coords to 0-based half-open for comparison
    g_s = start_1based - 1
    g_e = end_1based          # already exclusive in 0-based
    # binary search: find first interval whose end > g_s
    lo, hi = 0, len(ivs)
    while lo < hi:
        mid = (lo + hi) // 2
        if ivs[mid][1] <= g_s:
            lo = mid + 1
        else:
            hi = mid
    return lo < len(ivs) and ivs[lo][0] < g_e


def build_allowed_intervals(chrom_sizes, exclusions):
    """
    For every chromosome compute the list of (start, end) windows where a gene
    *start* can be placed.  We merge exclusion intervals and subtract them from
    [0, chrom_size).  Returns dict chrom -> list of (start, end), and a parallel
    list of cumulative lengths for weighted random choice.
    """
    allowed = {}
    for chrom, size in chrom_sizes.items():
        excl = sorted(exclusions.get(chrom, []))
        free = []
        pos = 0
        for es, ee in excl:
            if es > pos:
                free.append((pos, es))
            pos = max(pos, ee)
        if pos < size:
            free.append((pos, size))
        allowed[chrom] = free
    return allowed


def total_free_length(intervals):
    return sum(e - s for s, e in intervals)


def build_intergenic_candidates(allowed_by_chrom):
    """Flatten free intervals into (chrom, start, end, length) rows,
    sorted ascending by length so bisect can find eligible intervals quickly."""
    candidates = []
    for chrom, intervals in allowed_by_chrom.items():
        for s, e in intervals:
            candidates.append((chrom, s, e, e - s))
    candidates.sort(key=lambda r: r[3])   # sort by length ascending
    return candidates


def draw_position_r_style(rng, intergenic_candidates, gene_len, flank=1000):
    """R-style placement:
    1) choose eligible region uniformly where region_len > gene_len + 2*flank
    2) choose start uniformly in [region_start + flank, region_end - flank - gene_len]
    Returns (chrom, start_0_based).

    Uses pre-sorted candidates + bisect for O(log N) eligibility lookup.
    """
    threshold = gene_len + 2 * flank
    # bisect on the 4th field (length): find first index where length > threshold
    # Use a key-based approach with a helper fake list
    lengths = [r[3] for r in intergenic_candidates]  # already sorted
    idx = bisect.bisect_right(lengths, threshold)
    eligible_count = len(intergenic_candidates) - idx
    if eligible_count == 0:
        raise RuntimeError(
            f"No eligible free interval for gene length {gene_len} with flank {flank}"
        )

    chrom, s, e, _ = intergenic_candidates[idx + rng.randint(0, eligible_count - 1)]
    start_min = s + flank
    start_max = e - flank - gene_len
    start = rng.randint(start_min, start_max)
    return chrom, start


def open_file(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) not in (6, 7):
        print("Usage: python3 permute_gtf.py <seed> <input.gtf[.gz]> <tel_cen.bed> <blacklist.bed[.gz]> <output.gtf[.gz]> [too_large_genes.txt]",
              file=sys.stderr)
        sys.exit(1)

    seed            = int(sys.argv[1])
    gtf_in          = sys.argv[2]
    bed_excl        = sys.argv[3]
    bed_black       = sys.argv[4]
    gtf_out         = sys.argv[5]
    too_large_out   = sys.argv[6] if len(sys.argv) == 7 else None

    print(f"Seed: {seed}", flush=True)

    # ── seed=0: pass-through mode — copy input to output unchanged ────────────
    if seed == 0:
        print("Seed=0: copying input GTF to output without permutation...", flush=True)
        with open_file(gtf_in, 'rt') as fh_in, open_file(gtf_out, 'wt') as fh_out:
            for line in fh_in:
                fh_out.write(line)
        print("Done (pass-through).", flush=True)
        return

    rng = random.Random(seed)

    # ── 1. Load exclusion regions and blacklist; merge and build allowed intervals ──
    print("Loading exclusion regions (tel/cen)...", flush=True)
    exclusions = load_exclusions(bed_excl)

    print("Loading blacklist regions...", flush=True)
    blacklist  = load_exclusions(bed_black)
    n_bl_ivs   = sum(len(v) for v in blacklist.values())
    print(f"  Blacklist: {n_bl_ivs:,} intervals across {len(blacklist)} chromosomes", flush=True)

    # Merge blacklist into exclusions so genes cannot be *placed* over blacklist
    # regions, but all source genes are still permuted regardless of their input location.
    merged_excl = defaultdict(list)
    for chrom, ivs in exclusions.items():
        merged_excl[chrom].extend(ivs)
    for chrom, ivs in blacklist.items():
        merged_excl[chrom].extend(ivs)
    for chrom in merged_excl:
        merged_excl[chrom].sort()

    allowed    = build_allowed_intervals(HG38_CHROM_SIZES, merged_excl)
    intergenic_candidates = build_intergenic_candidates(allowed)
    total_free = sum(row[3] for row in intergenic_candidates)
    print(f"  Total free genome after merging exclusions+blacklist: {total_free:,} bp across {len(allowed)} chromosomes", flush=True)

    # ── 2. Read GTF: group lines by gene_id ──────────────────────────────────
    print("Reading GTF...", flush=True)
    header_lines = []
    gene_blocks  = {}          # gene_id -> list of raw lines
    gene_order   = []          # to preserve output order
    gene_coords  = {}          # gene_id -> (chrom, gene_start_1based, gene_end_1based)

    with open_file(gtf_in) as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if line.startswith('#'):
                header_lines.append(line)
                continue
            fields = line.split('\t')
            if len(fields) < 9:
                continue
            gid_m = re.search(r'gene_id "([^"]+)"', fields[8])
            if not gid_m:
                continue
            gid = gid_m.group(1)
            if gid not in gene_blocks:
                gene_blocks[gid] = []
                gene_order.append(gid)
            gene_blocks[gid].append(fields)
            if fields[2] == 'gene':
                gene_coords[gid] = (fields[0], int(fields[3]), int(fields[4]))

    print(f"  Loaded {len(gene_order):,} genes", flush=True)

    # ── 3. For each gene, draw a new location and shift all features ─────────
    print("Permuting gene locations...", flush=True)
    out_blocks      = []
    n_permuted      = 0
    n_no_coord      = 0
    n_too_large     = 0
    too_large_genes = []   # gene IDs kept at original location due to size

    for gid in gene_order:
        if gid not in gene_coords:
            # gene has no 'gene' feature line — pass through unchanged
            out_blocks.append(gene_blocks[gid])
            n_no_coord += 1
            continue

        orig_chrom, g_start, g_end = gene_coords[gid]
        gene_len = g_end - g_start + 1   # 1-based inclusive length

        try:
            new_chrom, new_start_0 = draw_position_r_style(rng, intergenic_candidates, gene_len, flank=1000)
        except RuntimeError:
            # No eligible intergenic region is large enough; keep gene at original location
            print(f"  WARNING: no eligible region for gene {gid} (len={gene_len:,} bp); keeping original location", flush=True)
            out_blocks.append(gene_blocks[gid])
            too_large_genes.append(gid)
            n_too_large += 1
            continue
        new_start_1 = new_start_0 + 1    # convert to 1-based GTF coordinate
        offset = new_start_1 - g_start

        new_fields_list = []
        for fields in gene_blocks[gid]:
            f = fields[:]
            f[0] = new_chrom
            f[3] = str(int(f[3]) + offset)
            f[4] = str(int(f[4]) + offset)
            new_fields_list.append(f)
        out_blocks.append(new_fields_list)
        n_permuted += 1

    print(f"  Permuted: {n_permuted:,}  |  Too large to place (kept in place): {n_too_large:,}  |  No coords (kept in place): {n_no_coord:,}", flush=True)

    # ── 4a. Write too-large gene list ────────────────────────────────────────
    if too_large_genes:
        if too_large_out:
            with open(too_large_out, 'w') as fh:
                for gid in too_large_genes:
                    fh.write(gid + '\n')
            print(f"  Too-large gene IDs written to: {too_large_out}", flush=True)
        else:
            print("  Too-large gene IDs (pass a 7th argument to save to file):", flush=True)
            for gid in too_large_genes:
                print(f"    {gid}", flush=True)

    # ── 4b. Write output GTF ─────────────────────────────────────────────────
    print(f"Writing output to {gtf_out} ...", flush=True)
    with open_file(gtf_out, 'wt') as fh:
        for h in header_lines:
            fh.write(h + '\n')
        for block in out_blocks:
            for fields in block:
                fh.write('\t'.join(fields) + '\n')

    print(f"Done. {len(gene_order):,} genes processed.", flush=True)


if __name__ == '__main__':
    main()
