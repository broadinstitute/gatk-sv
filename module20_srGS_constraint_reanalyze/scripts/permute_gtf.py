"""
Permute gene locations across the genome.

For each gene in the input GTF, randomly reassign it to a new location on any
chromosome, avoiding telomere/centromere regions.  All sub-features of the gene
(transcript, exon, CDS, UTR, etc.) are shifted by the same offset so internal
gene structure is preserved exactly.

Usage:
    python3 permute_gtf.py <seed> <input.gtf.gz> <tel_cen.bed> <output.gtf.gz>

Example:
    python3 permute_gtf.py 1 \
        genes_grch38_annotated_4_mapped_gencode_v39.CDS.gtf.gz \
        hg38_tel_cen.bed \
        genes_grch38_annotated_4_mapped_gencode_v39.CDS.permuted_seed1.gtf.gz
"""

import re, gzip, sys, random
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
    """Return dict chrom -> list of (start, end) 0-based half-open excluded intervals."""
    excl = defaultdict(list)
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            f = line.split('\t')
            excl[f[0]].append((int(f[1]), int(f[2])))
    for chrom in excl:
        excl[chrom].sort()
    return excl


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


def draw_position(rng, allowed_by_chrom, chroms, cum_lengths, gene_len):
    """
    Pick a random (chrom, start) such that [start, start+gene_len) fits fully
    within a free interval.  Falls back to retry if no interval on the chosen
    chrom can fit the gene.
    """
    total = cum_lengths[-1]
    for _ in range(10000):           # safety limit
        r = rng.randint(0, total - 1)
        # binary search to find which chromosome
        lo, hi = 0, len(cum_lengths) - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if cum_lengths[mid] <= r:
                lo = mid + 1
            else:
                hi = mid
        chrom = chroms[lo]
        intervals = allowed_by_chrom[chrom]
        # find intervals large enough for this gene
        fitting = [(s, e) for s, e in intervals if e - s >= gene_len]
        if not fitting:
            continue
        # pick a random fitting interval weighted by length
        lengths = [e - s - gene_len + 1 for s, e in fitting]
        total_fit = sum(lengths)
        r2 = rng.randint(0, total_fit - 1)
        acc = 0
        for (s, e), l in zip(fitting, lengths):
            acc += l
            if r2 < acc:
                return chrom, s + (r2 - (acc - l))
    raise RuntimeError(f"Could not place gene of length {gene_len} after 10000 tries")


def open_file(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) != 5:
        print("Usage: python3 permute_gtf.py <seed> <input.gtf[.gz]> <tel_cen.bed> <output.gtf[.gz]>",
              file=sys.stderr)
        sys.exit(1)

    seed      = int(sys.argv[1])
    gtf_in    = sys.argv[2]
    bed_excl  = sys.argv[3]
    gtf_out   = sys.argv[4]

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

    # ── 1. Load exclusion regions and build allowed intervals ────────────────
    print("Loading exclusion regions...", flush=True)
    exclusions = load_exclusions(bed_excl)
    allowed    = build_allowed_intervals(HG38_CHROM_SIZES, exclusions)

    chroms      = sorted(allowed.keys(), key=lambda c: (0, int(c[3:])) if c[3:].isdigit() else (1, c[3:]))
    free_lens   = [total_free_length(allowed[c]) for c in chroms]
    cum_lengths = []
    acc = 0
    for fl in free_lens:
        acc += fl
        cum_lengths.append(acc)
    print(f"  Total free genome: {acc:,} bp across {len(chroms)} chromosomes", flush=True)

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
    out_blocks = []

    for gid in gene_order:
        if gid not in gene_coords:
            # gene has no 'gene' feature line — pass through unchanged
            out_blocks.append(gene_blocks[gid])
            continue

        orig_chrom, g_start, g_end = gene_coords[gid]
        gene_len = g_end - g_start + 1   # 1-based inclusive length

        new_chrom, new_start_0 = draw_position(rng, allowed, chroms, cum_lengths, gene_len)
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

    # ── 4. Write output GTF ──────────────────────────────────────────────────
    print(f"Writing output to {gtf_out} ...", flush=True)
    with open_file(gtf_out, 'wt') as fh:
        for h in header_lines:
            fh.write(h + '\n')
        for block in out_blocks:
            for fields in block:
                fh.write('\t'.join(fields) + '\n')

    print(f"Done. {len(gene_order):,} genes permuted.", flush=True)


if __name__ == '__main__':
    main()
