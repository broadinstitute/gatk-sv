#!/bin/bash

set -eu -o pipefail

# Filter input $VCF by passing its records through $RECORDS_FILTER_CMD, keeping the header intact
# Store the result in file named $FILTERED_VCF_NAME
# #### IMPORTANT NOTE ####
# If using a filter command that may return non-zero in normal operation (e.g. grep returning no records)
# you MUST ensure RECORDS_FILTER_CMD returns non-zero e.g
#   RECORDS_FILTER_CMD="grep $VARIANT_ID"
# should be replaced with
#   RECORDS_FILTER_CMD="(grep $VARIANT_ID || printf '')"
# (The alternative would be to turn off pipefail, but that would result in filters with genuine errors
#  always escaping detection)

VCF="$1"
RECORDS_FILTER_CMD="$2"
FILTERED_VCF_NAME="$3"

# make a temp directory to store uncompressed info without any possibility of clobbering calling scripts
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp/}$(basename $0).XXXXXXXXXXXX")
# clean up TMP_DIR on exit
function clean_up () {
  rm -rf "$TMP_DIR"
}
trap clean_up EXIT

# uncompress vcf
UNFILTERED_VCF="$TMP_DIR/uncompressed_unfiltered.vcf"
zcat "$VCF" > "$UNFILTERED_VCF"

# Extract vcf header:
ONLY_HEADER=false
HEADER="$TMP_DIR/uncompressed_vcf_header.txt"
#  search for first line not starting with '#', stop immediately,
#  take everything up to that point, then remove last line.
#  if there are no non-header lines, just cat the whole vcf
grep -B9999999999 -m1 -Ev "^#" "$UNFILTERED_VCF" | sed '$ d' > "$HEADER" \
  || ONLY_HEADER=true

if $ONLY_HEADER; then
  # no records, so filter is trivial, just copy original vcf
  cp "$VCF" "$FILTERED_VCF_NAME"
else
  N_HEADER=$(wc -l < "$HEADER")

  # Read all the records by skipping the header and using tail to extract the rest
  tail -n+$((N_HEADER+1)) "$UNFILTERED_VCF" \
    | eval "$RECORDS_FILTER_CMD" \
    | cat "$HEADER" - \
    | vcf-sort \
    | bgzip -c \
    > "$FILTERED_VCF_NAME"
fi
