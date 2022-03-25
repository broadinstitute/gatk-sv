#!/bin/bash
set -eu -o pipefail
MODULE_ARCHIVE_URL=$1
# override for debugging purposes
QUIET=${QUIET:-TRUE}

ARCHIVE_DIR=$(mktemp -d "${TMPDIR:-/tmp/}$(basename $0).XXXXXXXXXXXX")
trap "rm -rf $ARCHIVE_DIR" EXIT

MODULE_ARCHIVE_DEST="$ARCHIVE_DIR/$(basename "$MODULE_ARCHIVE_URL")"
curl "$MODULE_ARCHIVE_URL" --output "$MODULE_ARCHIVE_DEST"
Rscript -e "install.packages('$MODULE_ARCHIVE_DEST', repos = NULL, quiet = $QUIET)"
