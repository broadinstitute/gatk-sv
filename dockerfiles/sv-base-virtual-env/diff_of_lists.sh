#!/bin/bash
set -eu -o pipefail
# given two lists (input as strings of words, or paths to files),
# output all the words in the first file that are not in the second file

if [[ $# -ne 2 ]]; then
  echo "diff_of_lists.sh requires exactly two arguments (perhaps an un-quoted string of lists was passed?)" >&2
  exit 1
fi

TEMPD=$(mktemp -d)
trap "rm -rf $TEMPD" EXIT

# print out elements of first list not in second
{ if [[ -f "$1" ]]; then cat "$1"; else echo "$1"; fi; } | tr -s ' ' '\n' | sort > $TEMPD/list.a
{ if [[ -f "$2" ]]; then cat "$2"; else echo "$2"; fi; } | tr -s ' ' '\n' | sort > $TEMPD/list.b

comm -23 $TEMPD/list.a $TEMPD/list.b | tr -s '\n' ' '
