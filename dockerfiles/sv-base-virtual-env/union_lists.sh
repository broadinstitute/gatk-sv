#!/bin/bash
set -eu -o pipefail

for LIST in "$@"; do
    if [[ -f $LIST ]]; then
      sed 's/ /\n/g' $LIST
      printf "\n"
    else
      echo "$LIST" | sed 's/ /\n/g'
    fi
done | grep -v ^$ | sort -u | tr -s '\n' ' ' | sed 's/ *$//'
