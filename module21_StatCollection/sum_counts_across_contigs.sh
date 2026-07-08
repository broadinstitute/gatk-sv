#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $(basename "$0") OUTPUT_TSV PER_CONTIG_COUNT_1 [PER_CONTIG_COUNT_2 ...]" >&2
  exit 1
fi

OUT=$1
shift

cat "$@" \
  | awk '
      BEGIN { OFS="\t" }
      { total[$1] += $2; low[$1] += $3 }
      END {
        print "sample_id", "total_variant_count", "af_lt_cutoff_variant_count"
        for (sample in total) {
          print sample, total[sample], low[sample]
        }
      }
    ' \
  | (head -n 1; tail -n +2 | sort -k1,1) \
  > "$OUT"
