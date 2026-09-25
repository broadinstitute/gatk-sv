#!/usr/bin/env bash
# Reproduce the CondenseReadCounts header-probe failure and prove the guard fixes it.
# The OLD block is the verbatim pre-fix text from wdl/CollectCoverage.wdl; the NEW block is
# extracted from the working file at run time, so this grades what the WDL actually contains.
set -uo pipefail
cd /tmp
awk 'BEGIN{print "@HD\tVN:1.3"; print "@RG\tID:run1\tSM:NA12878"; for(i=1;i<=4000000;i++) print "chr1\t" i "\t100\t42"}' | gzip > big.counts.tsv.gz
awk 'BEGIN{print "@HD\tVN:1.3"; print "@RG\tID:run1\tSM:NA12878"; for(i=1;i<=200;i++)  print "chr1\t" i "\t100\t42"}' | gzip > small.counts.tsv.gz
ls -l big.counts.tsv.gz small.counts.tsv.gz | awk '{printf "  %s %s bytes\n", $9, $5}'

run_block() {  # $1 = script file, $2 = counts file
  bash -c 'set -euxo pipefail
           export COUNTS="'"$2"'"
           source '"$1"'
           echo "VALUE=${existing_sample_id}"' 2>/dev/null
  return $?
}

cat > old_block.sh <<'OLD'
existing_sample_id=$(zcat "$COUNTS" | awk -F "\t" '/^@RG/ {
    for (i = 1; i <= NF; ++i) {
      if ($i ~ /^SM:/) {
        sub(/^SM:/, "", $i)
        print $i
        exit
      }
    }
  }')
OLD

export WDL="${WDL:-$(cd "$(dirname "$0")/../.." && pwd)/wdl/CollectCoverage.wdl}"
python3 - <<'PY' > new_block.sh
import os, re
txt = open(os.environ['WDL']).read()
m = re.search(r"existing_sample_id=\$\(zcat ~\{counts\}.*?\n\s*\}'[^)]*\)", txt, re.S)
assert m, "could not find the probe in the WDL -- the test would be grading nothing"
block = m.group(0).replace('~{counts}', '"$COUNTS"')
assert '|| true' in block, "the probe in the WDL still has no guard"
print(block)
PY

echo "--- control: SMALL input with the OLD block (expected 0 -- this is why the bug hid) ---"
run_block old_block.sh small.counts.tsv.gz; echo "  rc=$?"
echo "--- OLD block on the LARGE input (expected 141 = SIGPIPE via pipefail) ---"
run_block old_block.sh big.counts.tsv.gz; old_rc=$?; echo "  rc=$old_rc"
echo "--- NEW block, taken from wdl/CollectCoverage.wdl, on the LARGE input (expected 0) ---"
run_block new_block.sh big.counts.tsv.gz; new_rc=$?; echo "  rc=$new_rc"

fail=0
[ "$old_rc" = "141" ] || { echo "FAIL: old block did not exit 141 (got $old_rc) -- reproduction is not the reported defect"; fail=1; }
[ "$new_rc" = "0" ]   || { echo "FAIL: fixed block still fails (got $new_rc)"; fail=1; }
grep -q "VALUE=NA12878" <(run_block new_block.sh big.counts.tsv.gz) || { echo "FAIL: guard changed the extracted sample id"; fail=1; }
echo "=== reproduction=$([ "$old_rc" = 141 ] && echo yes || echo no)  fixed=$([ "$new_rc" = 0 ] && echo yes || echo no)  value preserved ==="
exit $fail
