"""
Split a GTF file into 7 genomic annotation category BED files:
  CDS, coding_transcript, intron, promoter, transcript, utr_5, utr_3

Definitions (all merged per gene across all transcripts):
  CDS              : CDS features (merged across isoforms per gene)
  coding_transcript: merged exons of protein_coding genes
  intron           : gaps between merged exons per gene
  promoter         : 2 kb upstream of TSS (strand-aware)
                     + strand: [gene_start - 2000, gene_start)
                     - strand: (gene_end, gene_end + 2000]
  transcript       : full gene body (gene start to gene end)
  utr_5            : 5' UTR (strand-aware; UTR midpoint upstream of CDS)
  utr_3            : 3' UTR (strand-aware; UTR midpoint downstream of CDS)

Output BED columns: chr, start(0-based), end, strand, gene_id, gene_type, gene_name

Usage:
    python3 split_gtf_annotations.py <input.gtf[.gz]> <output_prefix>

Example:
    python3 split_gtf_annotations.py genes.gtf.gz ./annot/genes
    → writes genes.CDS.bed, genes.coding_transcript.bed, etc.
"""

import re, gzip, sys
from collections import defaultdict

# hg38 chromosome sizes (UCSC) — used to clamp promoter boundaries
HG38_CHROM_SIZES = {
    'chr1':  248956422, 'chr2':  242193529, 'chr3':  198295559,
    'chr4':  190214555, 'chr5':  181538259, 'chr6':  170805979,
    'chr7':  159345973, 'chr8':  145138636, 'chr9':  138394717,
    'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
    'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
    'chr16':  90338345, 'chr17':  83257441, 'chr18':  80373285,
    'chr19':  58617616, 'chr20':  64444167, 'chr21':  46709983,
    'chr22':  50818468, 'chrX':  156040895, 'chrY':   57227415,
}


def open_file(path, mode='rt'):
    return gzip.open(path, mode) if path.endswith('.gz') else open(path, mode)


def parse_attr(attrs, key):
    m = re.search(r'%s "([^"]+)"' % key, attrs)
    return m.group(1) if m else '.'


def merge_intervals(intervals):
    """Merge overlapping/adjacent 1-based inclusive intervals. Returns sorted merged list."""
    if not intervals:
        return []
    merged = sorted(intervals)
    result = [list(merged[0])]
    for s, e in merged[1:]:
        if s <= result[-1][1] + 1:
            result[-1][1] = max(result[-1][1], e)
        else:
            result.append([s, e])
    return result


def chrom_sort_key(row):
    c = row[0].replace('chr', '')
    chrom_key = (0, int(c)) if c.isdigit() else (1, c)
    return (chrom_key, int(row[1]), int(row[2]))


def write_bed(rows, path):
    rows_sorted = sorted(rows, key=chrom_sort_key)
    with open(path, 'w') as fh:
        fh.write('#chr\tstart\tend\tstrand\tgene_id\tgene_type\tgene_name\n')
        for r in rows_sorted:
            fh.write('\t'.join(str(x) for x in r) + '\n')
    return len(rows_sorted)


def main():
    if len(sys.argv) != 3:
        print("Usage: python3 split_gtf_annotations.py <input.gtf[.gz]> <output_prefix>",
              file=sys.stderr)
        sys.exit(1)

    gtf_in = sys.argv[1]
    prefix = sys.argv[2]

    # ── 1. Parse GTF ─────────────────────────────────────────────────────────
    print(f"Reading {gtf_in} ...", flush=True)

    gene_info = {}           # gene_id -> {chrom, strand, gene_type, gene_name, start, end}
    gene_exons = defaultdict(list)  # gene_id -> [(start, end)] 1-based inclusive
    gene_cds   = defaultdict(list)
    gene_utrs  = defaultdict(list)

    with open_file(gtf_in) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9:
                continue
            chrom, _, feat, start, end, _, strand, _, attrs = f
            start, end = int(start), int(end)

            gene_id   = parse_attr(attrs, 'gene_id')
            gene_type = parse_attr(attrs, 'gene_type')
            gene_name = parse_attr(attrs, 'gene_name')

            if feat == 'gene':
                gene_info[gene_id] = dict(
                    chrom=chrom, strand=strand,
                    gene_type=gene_type, gene_name=gene_name,
                    start=start, end=end
                )
            elif feat == 'exon':
                gene_exons[gene_id].append((start, end))
            elif feat == 'CDS':
                gene_cds[gene_id].append((start, end))
            elif feat == 'UTR':
                gene_utrs[gene_id].append((start, end))

    print(f"  Loaded {len(gene_info):,} genes", flush=True)

    # ── 2. Build output rows per category ────────────────────────────────────
    out = {k: [] for k in ('CDS', 'coding_transcript', 'intron',
                            'promoter', 'transcript', 'utr_5', 'utr_3')}

    for gene_id, info in gene_info.items():
        chrom     = info['chrom']
        strand    = info['strand']
        gene_type = info['gene_type']
        gene_name = info['gene_name']
        g_start   = info['start']   # 1-based
        g_end     = info['end']     # 1-based inclusive
        meta      = (strand, gene_id, gene_type, gene_name)

        # ── transcript: full gene body ────────────────────────────────────────
        out['transcript'].append((chrom, g_start - 1, g_end, *meta))

        # ── promoter: 2 kb upstream of TSS ───────────────────────────────────
        chrom_size = HG38_CHROM_SIZES.get(chrom, 999_999_999)
        if strand == '+':
            tss_0 = g_start - 1          # 0-based TSS position
            p_start = max(0, tss_0 - 2000)
            p_end   = tss_0
        else:                            # - strand: TSS is at the high-coord end
            tss_0 = g_end                # 0-based one-past-end = TSS for - strand
            p_start = tss_0
            p_end   = min(chrom_size, tss_0 + 2000)
        if p_start < p_end:
            out['promoter'].append((chrom, p_start, p_end, *meta))

        # ── CDS: merged CDS features ──────────────────────────────────────────
        if gene_cds[gene_id]:
            for s, e in merge_intervals(gene_cds[gene_id]):
                out['CDS'].append((chrom, s - 1, e, *meta))

        # ── coding_transcript: merged exons of protein_coding genes ──────────
        if gene_exons[gene_id] and gene_type == 'protein_coding':
            for s, e in merge_intervals(gene_exons[gene_id]):
                out['coding_transcript'].append((chrom, s - 1, e, *meta))

        # ── intron: gaps between merged exons ─────────────────────────────────
        if gene_exons[gene_id]:
            merged_ex = merge_intervals(gene_exons[gene_id])
            for i in range(len(merged_ex) - 1):
                i_start = merged_ex[i][1] + 1       # 1-based intron start
                i_end   = merged_ex[i + 1][0] - 1   # 1-based intron end
                if i_start <= i_end:
                    out['intron'].append((chrom, i_start - 1, i_end, *meta))

        # ── UTR: classify to 5' or 3' based on strand + CDS boundaries ───────
        if gene_utrs[gene_id] and gene_cds[gene_id]:
            cds_min = min(s for s, e in gene_cds[gene_id])
            cds_max = max(e for s, e in gene_cds[gene_id])

            for u_s, u_e in gene_utrs[gene_id]:
                mid = (u_s + u_e) / 2
                if strand == '+':
                    # 5' UTR: upstream of CDS (lower coords)
                    if mid < cds_min:
                        out['utr_5'].append((chrom, u_s - 1, u_e, *meta))
                    # 3' UTR: downstream of CDS (higher coords)
                    elif mid > cds_max:
                        out['utr_3'].append((chrom, u_s - 1, u_e, *meta))
                else:  # - strand: TSS is at high coords
                    # 5' UTR: upstream of CDS = higher coords
                    if mid > cds_max:
                        out['utr_5'].append((chrom, u_s - 1, u_e, *meta))
                    # 3' UTR: downstream of CDS = lower coords
                    elif mid < cds_min:
                        out['utr_3'].append((chrom, u_s - 1, u_e, *meta))

    # ── 3. Write output ───────────────────────────────────────────────────────
    print("\nWriting output BED files:", flush=True)
    for name in ('CDS', 'coding_transcript', 'intron', 'promoter',
                 'transcript', 'utr_5', 'utr_3'):
        path = f"{prefix}.{name}.bed"
        n = write_bed(out[name], path)
        print(f"  {name:<20s}: {n:>8,} intervals  →  {path}")

    print("\nDone.", flush=True)


if __name__ == '__main__':
    main()
