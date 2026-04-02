"""
Extract introns from a GTF file.
Introns are defined as gaps between CDS and UTR features within each gene.
CDS + UTR together cover all exonic sequence; gaps between them are intronic.

Output BED columns: chr, start(0-based), end(1-based), strand, gene_id, gene_type, gene_name

Usage: python3 gtf_to_introns.py <input.gtf[.gz]> <output.bed>
"""

import re
import gzip
import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print("Usage: python3 gtf_to_introns.py <input.gtf[.gz]> <output.bed>", file=sys.stderr)
    sys.exit(1)

GTF = sys.argv[1]
OUT = sys.argv[2]


def parse_attr(attrs, key):
    m = re.search(r'%s "([^"]+)"' % key, attrs)
    return m.group(1) if m else '.'


# Group by gene_id: collect all CDS + UTR intervals per gene
genes = defaultdict(lambda: {
    'chrom': None, 'strand': None,
    'gene_type': None, 'gene_name': None,
    'intervals': []
})

opener = gzip.open if GTF.endswith('.gz') else open
print("Reading GTF...", flush=True)
with opener(GTF, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if len(fields) < 9:
            continue
        chrom, source, feat, start, end, score, strand, frame, attrs = fields
        if feat not in ('CDS', 'UTR'):
            continue

        gene_id   = parse_attr(attrs, 'gene_id')
        gene_type = parse_attr(attrs, 'gene_type')
        gene_name = parse_attr(attrs, 'gene_name')

        g = genes[gene_id]
        g['chrom']     = chrom
        g['strand']    = strand
        g['gene_type'] = gene_type
        g['gene_name'] = gene_name
        g['intervals'].append((int(start), int(end)))  # 1-based, inclusive

print(f"Loaded {len(genes)} genes. Writing introns...", flush=True)

n_introns = 0
with open(OUT, 'w') as out:
    out.write('\t'.join(['#chr', 'start', 'end', 'strand', 'gene_id', 'gene_type', 'gene_name']) + '\n')

    # Sort genes by chrom then by first interval start for tidy output
    def gene_sort_key(item):
        g = item[1]
        chrom = g['chrom'] or ''
        # numeric sort for chr1..22, then chrX, chrY, chrM
        c = chrom.replace('chr', '')
        try:
            return (0, int(c), g['intervals'][0][0] if g['intervals'] else 0)
        except ValueError:
            return (1, c, g['intervals'][0][0] if g['intervals'] else 0)

    for gene_id, info in sorted(genes.items(), key=gene_sort_key):
        intervals = info['intervals']
        if len(intervals) < 2:
            continue  # single-exon gene — no introns

        # Merge overlapping or adjacent CDS/UTR intervals into exonic blocks
        merged = sorted(intervals)
        exons = [list(merged[0])]
        for s, e in merged[1:]:
            if s <= exons[-1][1] + 1:   # overlapping or immediately adjacent
                exons[-1][1] = max(exons[-1][1], e)
            else:
                exons.append([s, e])

        if len(exons) < 2:
            continue  # all intervals merged into one block — no introns

        # Gaps between consecutive exonic blocks = introns
        for i in range(len(exons) - 1):
            istart_1based = exons[i][1] + 1       # 1-based start of intron
            iend_1based   = exons[i + 1][0] - 1   # 1-based end   of intron
            if istart_1based > iend_1based:
                continue
            # Convert to 0-based half-open BED coordinates
            bed_start = istart_1based - 1
            bed_end   = iend_1based
            out.write('\t'.join([
                info['chrom'], str(bed_start), str(bed_end),
                info['strand'], gene_id, info['gene_type'], info['gene_name']
            ]) + '\n')
            n_introns += 1

print(f"Done. Wrote {n_introns} introns to:\n  {OUT}")
