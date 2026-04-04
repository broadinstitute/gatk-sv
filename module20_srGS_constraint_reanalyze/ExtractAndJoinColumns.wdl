version 1.0

import "Structs.wdl"

# For each column X in [37..101]:
#   - from each input tsv.gz (ordered by column_labels), extract cols 1-25 + col X
#   - join all files on cols 1-25, with one additional column per file (named by column_label)
#   - output: one joined tsv.gz per column X
workflow ExtractAndJoinColumns {
  input {
    Array[File]   tsv_gz_files    # same order as column_labels
    Array[String] column_labels   # one label per file; becomes the column header for col X
    String        docker          # docker with python3
    Int           first_col  = 37  # first SV feature column (1-based), inclusive
    Int           last_col   = 101 # last SV feature column (1-based), inclusive
    Int           key_cols   = 25  # number of leading key columns (1-based: 1..key_cols)
    RuntimeAttr?  runtime_attr_override
  }

  # scatter over column indices 37..101  (range gives 0..64, add offset)
  Int n_cols = last_col - first_col + 1   # = 65

  scatter (idx in range(n_cols)) {
    Int col_num = first_col + idx          # 37, 38, ..., 101

    call ExtractAndJoin {
      input:
        tsv_gz_files          = tsv_gz_files,
        column_labels         = column_labels,
        col_num               = col_num,
        key_cols              = key_cols,
        docker                = docker,
        runtime_attr_override = runtime_attr_override
    }
  }

  output {
    Array[File] joined_tables = ExtractAndJoin.joined_table
  }
}


# ---------------------------------------------------------------------------
# Task: extract cols 1..key_cols + col_num from every tsv.gz,
#       join all on cols 1..key_cols, output one wide table.
#
# Output column order:
#   col1 .. col<key_cols>   <column_labels[0]>  <column_labels[1]>  ...
# ---------------------------------------------------------------------------
task ExtractAndJoin {
  input {
    Array[File]   tsv_gz_files
    Array[String] column_labels
    Int           col_num       # 1-based column to extract from each file
    Int           key_cols      # number of leading key columns (1-based: 1..key_cols)
    String        docker
    RuntimeAttr?  runtime_attr_override
  }

  String out_file = "joined_col_" + col_num + ".tsv.gz"

  command <<<
    set -euo pipefail

    python3 <<'PYEOF'
import gzip
import sys
from collections import OrderedDict

tsv_files   = "~{sep=',' tsv_gz_files}".split(',')
col_labels  = "~{sep=',' column_labels}".split(',')
col_num     = ~{col_num}   # 1-based
key_cols    = ~{key_cols}  # number of key columns (1-based 1..key_cols)
out_file    = "~{out_file}"

col_idx     = col_num - 1             # 0-based index for the SV feature column
key_indices = list(range(key_cols))   # 0-based indices for key columns

if len(tsv_files) != len(col_labels):
    sys.exit(f"ERROR: {len(tsv_files)} files but {len(col_labels)} labels")

# Read each file: build dict keyed by tuple(key_cols) -> value at col_num
# First file also defines key column names
key_col_names = None
data = []   # list of dicts: key_tuple -> value

for fi, (fpath, label) in enumerate(zip(tsv_files, col_labels)):
    opener = gzip.open if fpath.endswith('.gz') else open
    d = OrderedDict()
    with opener(fpath, 'rt') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if len(header) <= col_idx:
            sys.exit(f"ERROR: file {fpath} has only {len(header)} columns; "
                     f"cannot extract column {col_num} (0-based {col_idx})")
        if key_col_names is None:
            key_col_names = [header[k] for k in key_indices]
        sv_col_name = header[col_idx]
        for line in fh:
            fields = line.rstrip('\n').split('\t')
            if len(fields) <= col_idx:
                continue
            key = tuple(fields[k] for k in key_indices)
            d[key] = fields[col_idx]
    data.append((label, d))

# Union of all keys, preserving first-file order, then adding extras
all_keys = list(data[0][1].keys())
seen = set(all_keys)
for label, d in data[1:]:
    for k in d:
        if k not in seen:
            all_keys.append(k)
            seen.add(k)

# Write output
with gzip.open(out_file, 'wt') as out:
    header_line = '\t'.join(key_col_names) + '\t' + '\t'.join(l for l, _ in data)
    out.write(header_line + '\n')
    for key in all_keys:
        row = list(key) + [d.get(key, 'NA') for _, d in data]
        out.write('\t'.join(row) + '\n')

print(f"Written {len(all_keys)} rows -> {out_file}", flush=True)
PYEOF
  >>>

  output {
    File joined_table = out_file
  }

  RuntimeAttr default_attr = object {
    cpu_cores:         2,
    mem_gb:            8,
    disk_gb:           20 + ceil(size(tsv_gz_files, "GiB") * 4),
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }

  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  runtime {
    cpu:            select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
    memory:         select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
    disks:          "local-disk " + select_first([runtime_attr.disk_gb, default_attr.disk_gb]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb,     default_attr.boot_disk_gb])
    docker:         docker
    preemptible:    select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
    maxRetries:     select_first([runtime_attr.max_retries,       default_attr.max_retries])
  }
}
