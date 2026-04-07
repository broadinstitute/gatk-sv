version 1.0

import "Structs.wdl"

# Workflow: for each SV feature column (first_col..last_col), scatter a task that:
#   1. Extracts that column + gene-key cols from training permutation files (seeds 1-N)
#      and from the prediction-seed file (seed 0).
#   2. Builds per-gene regression models (sv ~ SD + SR) across N training permutations.
#   3. Computes per-gene mean, median, and SD of the SV count across N training permutations.
#   4. Outputs four per-column tsv.gz files: lm_pred, mean, median, sd.
# A final gather task joins the per-column results into four genome-wide tables.
workflow RegressPermutationsPerColumn {
  input {
    # One tsv.gz per training permutation (seeds 1-N), in any order.
    # Must all have the same columns; col 1 = gene name.
    Array[File]  train_files           # e.g. seeds 1-100
    File         predict_file          # seed-0 tsv.gz (used for observed + SD/SR lookup)
    File         sd_table              # gene x seed columns: SD overlap proportions
    File         sr_table              # gene x seed columns: SR overlap proportions

    Int    first_col     = 37          # 1-based first SV feature column
    Int    last_col      = 426         # 1-based last SV feature column
    Int    key_cols      = 1           # number of leading key columns (gene id = col 1)
    String gene_sep      = ".permu"    # strip this suffix from gene names in col 1
    String predict_col   = "0"         # column name in sd_table / sr_table for predict seed
    Int    min_nonzero   = 5           # min non-zero obs per gene to attempt regression

    String r_docker                    # docker image with R (base R sufficient)

    RuntimeAttr? runtime_attr_regress
    RuntimeAttr? runtime_attr_gather
  }

  Int n_cols = last_col - first_col + 1

  # ── per-column scatter ────────────────────────────────────────────────────
  scatter (idx in range(n_cols)) {
    Int col_1based = first_col + idx   # 37, 38, ..., last_col

    call RegressOneColumn {
      input:
        train_files    = train_files,
        predict_file   = predict_file,
        sd_table       = sd_table,
        sr_table       = sr_table,
        col_1based     = col_1based,
        key_cols       = key_cols,
        gene_sep       = gene_sep,
        predict_col    = predict_col,
        min_nonzero    = min_nonzero,
        docker         = r_docker,
        runtime_attr_override = runtime_attr_regress
    }
  }

  # ── gather: join per-column outputs into four genome-wide tables ──────────
  call GatherColumnResults as GatherMean {
    input:
      per_col_files = RegressOneColumn.out_mean,
      out_name      = "permu_mean.tsv.gz",
      docker        = r_docker,
      runtime_attr_override = runtime_attr_gather
  }

  call GatherColumnResults as GatherMedian {
    input:
      per_col_files = RegressOneColumn.out_median,
      out_name      = "permu_median.tsv.gz",
      docker        = r_docker,
      runtime_attr_override = runtime_attr_gather
  }

  call GatherColumnResults as GatherSd {
    input:
      per_col_files = RegressOneColumn.out_sd,
      out_name      = "permu_sd.tsv.gz",
      docker        = r_docker,
      runtime_attr_override = runtime_attr_gather
  }

  call GatherColumnResults as GatherLmPred {
    input:
      per_col_files = RegressOneColumn.out_lm_pred,
      out_name      = "permu_lm_pred.tsv.gz",
      docker        = r_docker,
      runtime_attr_override = runtime_attr_gather
  }

  output {
    File mean_table    = GatherMean.merged_tsv_gz
    File median_table  = GatherMedian.merged_tsv_gz
    File sd_table_out  = GatherSd.merged_tsv_gz
    File lm_pred_table = GatherLmPred.merged_tsv_gz
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task 1: process a single SV-feature column
# ─────────────────────────────────────────────────────────────────────────────
task RegressOneColumn {
  input {
    Array[File] train_files
    File        predict_file
    File        sd_table
    File        sr_table
    Int         col_1based
    Int         key_cols    = 1
    String      gene_sep    = ".permu"
    String      predict_col = "0"
    Int         min_nonzero = 5
    String      docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    mem_gb:            4,
    cpu_cores:         1,
    disk_gb:           20,
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    Rscript - <<'REOF'
    col_idx   <- ~{col_1based}
    key_cols  <- ~{key_cols}
    gene_sep  <- "~{gene_sep}"
    pred_col  <- "~{predict_col}"
    min_nz    <- ~{min_nonzero}

    # ── helpers ──────────────────────────────────────────────────────────────
    strip_suffix <- function(x, sep) {
      idx <- regexpr(sep, x, fixed = TRUE)
      ifelse(idx > 0, substr(x, 1, idx - 1), x)
    }

    read_gz <- function(path) {
      read.table(gzfile(path), header = TRUE, sep = "\t",
                 comment.char = "", check.names = FALSE, quote = "")
    }

    # ── inputs ───────────────────────────────────────────────────────────────
    train_paths   <- readLines("~{write_lines(train_files)}")
    predict_path  <- "~{predict_file}"
    sd_path       <- "~{sd_table}"
    sr_path       <- "~{sr_table}"

    # ── load SD / SR tables ────────────────────────────────────────────────
    sd_tab <- read_gz(sd_path); rownames(sd_tab) <- sd_tab[[1]]
    sr_tab <- read_gz(sr_path); rownames(sr_tab) <- sr_tab[[1]]
    sd_tab = sd_tab[,-1]
    sr_tab = sr_tab[,-1]


    if (!pred_col %in% colnames(sd_tab)) stop("'", pred_col, "' not in SD table")
    if (!pred_col %in% colnames(sr_tab)) stop("'", pred_col, "' not in SR table")

    # ── load training files and extract target column + gene names ─────────
    # Result: mat_sv[gene, seed], mat_sd[gene, seed], mat_sr[gene, seed]
    cat(sprintf("Loading %d training files for col %d...\n", length(train_paths), col_idx))

    gene_names <- NULL
    seed_labels <- character(length(train_paths))

    sv_list <- vector("list", length(train_paths))
    for (i in seq_along(train_paths)) {
      dat <- read_gz(train_paths[i])
      if (ncol(dat) < col_idx) {
        warning("File ", basename(train_paths[i]), " has < ", col_idx, " cols — skipping")
        next
      }
      genes_raw <- strip_suffix(as.character(dat[[1]]), gene_sep)
      if (is.null(gene_names)) gene_names <- genes_raw
      sv_vec <- as.numeric(dat[[col_idx]])
      names(sv_vec) <- genes_raw
      sv_list[[i]] <- sv_vec
      seed_labels[i] <- as.character(i)   # track by index; seed label not needed
    }

    sv_list <- sv_list[!sapply(sv_list, is.null)]
    n_train <- length(sv_list)
    cat(sprintf("  %d training files loaded\n", n_train))

    if (is.null(gene_names)) stop("No training files successfully loaded")

    # Intersect genes across all sources
    common_genes <- Reduce(intersect, c(
      list(gene_names, rownames(sd_tab), rownames(sr_tab)),
      lapply(sv_list, names)
    ))
    cat(sprintf("  Common genes: %d\n", length(common_genes)))
    if (length(common_genes) == 0) stop("No common genes")

    # Build sv matrix: genes x train seeds
    sv_mat <- matrix(NA_real_, nrow = length(common_genes), ncol = n_train,
                     dimnames = list(common_genes, seq_len(n_train)))
    for (i in seq_along(sv_list)) {
      sv_mat[common_genes, i] <- sv_list[[i]][common_genes]
    }

    # SD/SR matrices for training seeds (use seed column names = 1..n_train)
    # Seed column labels in sd/sr table are the original seed numbers, so we use the
    # column names present in both train file list and sd/sr table.
    # We match by position: seed label from sv_list index vs. sd_tab column names.
    # For flexibility, match by string index of training-seed columns in sd/sr table.
    train_seed_cols <- colnames(sd_tab)[colnames(sd_tab) != pred_col]
    # keep only those present, up to n_train
    train_seed_cols <- head(train_seed_cols, n_train)

    sd_train <- as.matrix(sd_tab[common_genes, train_seed_cols, drop = FALSE])
    sr_train <- as.matrix(sr_tab[common_genes, train_seed_cols, drop = FALSE])

    # Prediction SD/SR (seed 0)
    sd0 <- as.numeric(sd_tab[common_genes, pred_col]); names(sd0) <- common_genes
    sr0 <- as.numeric(sr_tab[common_genes, pred_col]); names(sr0) <- common_genes

    # ── per-gene statistics ────────────────────────────────────────────────
    # Align sv_mat columns with sd_train/sr_train columns (both n_train wide)
    n_use <- min(ncol(sv_mat), ncol(sd_train))
    sv_use <- sv_mat[, seq_len(n_use), drop = FALSE]
    sd_use <- sd_train[, seq_len(n_use), drop = FALSE]
    sr_use <- sr_train[, seq_len(n_use), drop = FALSE]

    gene_mean   <- rowMeans(sv_use, na.rm = TRUE)
    gene_median <- apply(sv_use, 1, median, na.rm = TRUE)
    gene_sd     <- apply(sv_use, 1, sd,     na.rm = TRUE)

    # ── per-gene regression ────────────────────────────────────────────────
    gene_lm_pred <- rep(NA_real_, length(common_genes))
    names(gene_lm_pred) <- common_genes

    for (g in common_genes) {
      sv_vec <- sv_use[g, ]
      sd_vec <- sd_use[g, ]
      sr_vec <- sr_use[g, ]
      ok <- complete.cases(sv_vec, sd_vec, sr_vec)
      if (sum(ok) < 3) next
      if (sum(sv_vec[ok] != 0) < min_nz) next
      fit <- tryCatch(
        lm(sv ~ sd + sr, data = data.frame(sv = sv_vec[ok], sd = sd_vec[ok], sr = sr_vec[ok])),
        error = function(e) NULL
      )
      if (is.null(fit)) next
      gene_lm_pred[g] <- tryCatch(
        as.numeric(predict(fit, newdata = data.frame(sd = sd0[g], sr = sr0[g]))),
        error = function(e) NA_real_
      )
    }

    # ── load observed predict-seed column ─────────────────────────────────
    dat_pred <- read_gz(predict_path)
    if (ncol(dat_pred) < col_idx) stop("predict_file has < ", col_idx, " cols")
    genes_pred <- strip_suffix(as.character(dat_pred[[1]]), gene_sep)
    sv_pred <- as.numeric(dat_pred[[col_idx]]); names(sv_pred) <- genes_pred

    col_name <- colnames(dat_pred)[col_idx]
    cat(sprintf("  Column name: %s\n", col_name))

    # ── write four output files ────────────────────────────────────────────
    write_out <- function(vals, fname) {
      df <- data.frame(gene = names(vals), value = vals, stringsAsFactors = FALSE)
      colnames(df)[2] <- col_name
      write.table(df, gzfile(fname), sep = "\t", quote = FALSE, row.names = FALSE)
      cat("  Wrote:", fname, "\n")
    }

    write_out(gene_mean,    sprintf("col_%d_mean.tsv.gz",    col_idx))
    write_out(gene_median,  sprintf("col_%d_median.tsv.gz",  col_idx))
    write_out(gene_sd,      sprintf("col_%d_sd.tsv.gz",      col_idx))
    write_out(gene_lm_pred, sprintf("col_%d_lm_pred.tsv.gz", col_idx))

    cat("Done.\n")
    REOF
  >>>

  output {
    File out_mean    = "col_~{col_1based}_mean.tsv.gz"
    File out_median  = "col_~{col_1based}_median.tsv.gz"
    File out_sd      = "col_~{col_1based}_sd.tsv.gz"
    File out_lm_pred = "col_~{col_1based}_lm_pred.tsv.gz"
  }

  runtime {
    docker:       docker
    memory:       select_first([runtime_attr.mem_gb, 4]) + " GB"
    cpu:          select_first([runtime_attr.cpu_cores, 1])
    disks:        "local-disk " + select_first([runtime_attr.disk_gb, 20]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, 10])
    preemptible:  select_first([runtime_attr.preemptible_tries, 3])
    maxRetries:   select_first([runtime_attr.max_retries, 1])
  }
}

# ─────────────────────────────────────────────────────────────────────────────
# Task 2: join per-column two-column tsv.gz files → one wide gene x col table
# ─────────────────────────────────────────────────────────────────────────────
task GatherColumnResults {
  input {
    Array[File] per_col_files    # each: gene | colname  (two columns)
    String      out_name         # output filename
    String      docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr default_attr = object {
    mem_gb:            8,
    cpu_cores:         2,
    disk_gb:           30,
    boot_disk_gb:      10,
    preemptible_tries: 3,
    max_retries:       1
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])

  command <<<
    set -euo pipefail

    python3 - <<'PYEOF'
import gzip, sys, os

files = open("~{write_lines(per_col_files)}").read().splitlines()
files = [f for f in files if f]

# Read all files; each has header: gene <col_name>
tables = {}
col_order = []
all_genes = None

for f in files:
    with gzip.open(f, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        gene_col_name, val_col_name = header[0], header[1]
        col_order.append(val_col_name)
        col_data = {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            col_data[parts[0]] = parts[1]
        tables[val_col_name] = col_data
        if all_genes is None:
            all_genes = list(col_data.keys())

if all_genes is None:
    sys.exit("No input files loaded")

out_name = "~{out_name}"
with gzip.open(out_name, "wt") as out:
    out.write("gene\t" + "\t".join(col_order) + "\n")
    for gene in all_genes:
        vals = [tables[c].get(gene, "NA") for c in col_order]
        out.write(gene + "\t" + "\t".join(vals) + "\n")

print(f"Wrote {out_name}: {len(all_genes)} genes x {len(col_order)} columns")
PYEOF
  >>>

  output {
    File merged_tsv_gz = out_name
  }

  runtime {
    docker:       docker
    memory:       select_first([runtime_attr.mem_gb, 8]) + " GB"
    cpu:          select_first([runtime_attr.cpu_cores, 2])
    disks:        "local-disk " + select_first([runtime_attr.disk_gb, 30]) + " HDD"
    bootDiskSizeGb: select_first([runtime_attr.boot_disk_gb, 10])
    preemptible:  select_first([runtime_attr.preemptible_tries, 3])
    maxRetries:   select_first([runtime_attr.max_retries, 1])
  }
}
