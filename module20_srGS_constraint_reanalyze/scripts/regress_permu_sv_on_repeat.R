#!/usr/bin/env Rscript
# regress_permu_sv_on_repeat.R
#
# For each SV feature column (default: 37+):
#   For each gene:
#     - Collect N_train data points, one per training permutation:
#         sv_count[gene, seed],  SD_overlap[gene, seed],  SR_overlap[gene, seed]
#     - Fit: sv_count ~ SD + SR  (N_train observations = one model per gene)
#     - Predict: plug in SD[gene, "0"] and SR[gene, "0"] to get expected sv_count
#   Output tables of predicted vs observed for seed 0, plus per-gene model summaries.
#
# Outputs:
#   sv_predicted_seed0.tsv.gz  -- gene x SV_feature: model-predicted values for seed 0
#   sv_observed_seed0.tsv.gz   -- gene x SV_feature: observed values from seed 0
#   regression_summary.tsv     -- gene x SV_feature: intercept, beta_SD, beta_SR, R², N

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-d", "--dir"),        type = "character", default = NULL,
              help = "Directory containing permuted tsv.gz files and SD/SR overlap tables",
              metavar = "character"),
  make_option(c("--sd"),               type = "character", default = NULL,
              help = "gene_SD_overlap.tsv.gz path [default: <dir>/gene_SD_overlap.tsv.gz]",
              metavar = "character"),
  make_option(c("--sr"),               type = "character", default = NULL,
              help = "gene_SR_overlap.tsv.gz path [default: <dir>/gene_SR_overlap.tsv.gz]",
              metavar = "character"),
  make_option(c("--pattern"),          type = "character",
              default = "gnomAD_SV_v3\\.vs\\.hg38_gencode_v39\\.permuted_seed(\\d+)\\.tsv\\.gz$",
              help    = "Regex with one capture group for seed number [default: gnomAD naming]",
              metavar = "character"),
  make_option(c("--gene-sep"),         type = "character", default = ".permu",
              help    = "Separator to strip seed suffix from gene names in col 1 [default: .permu]",
              metavar = "character"),
  make_option(c("--first-sv-col"),     type = "integer",   default = 37L,
              help    = "1-based index of first SV feature column [default: 37]",
              metavar = "integer"),
  make_option(c("--train-seeds"),      type = "character", default = "1:100",
              help    = "Seed range for training, e.g. '1:100' or '1,2,3' [default: 1:100]",
              metavar = "character"),
  make_option(c("--predict-seed"),     type = "integer",   default = 0L,
              help    = "Seed to predict (SD/SR column matching this number) [default: 0]",
              metavar = "integer"),
  make_option(c("--min-nonzero"),      type = "integer",   default = 5L,
              help    = "Min non-zero SV counts across training seeds to attempt regression [default: 5]",
              metavar = "integer"),
  make_option(c("-o", "--output-dir"), type = "character", default = ".",
              help    = "Output directory [default: current dir]",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$dir)) stop("--dir is required")

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ── helpers ──────────────────────────────────────────────────────────────────
parse_seed_range <- function(s) {
  if (grepl("^\\d+:\\d+$", s)) { r <- as.integer(strsplit(s, ":")[[1]]); return(r[1]:r[2]) }
  as.integer(strsplit(s, ",")[[1]])
}

strip_suffix <- function(x, sep) {
  idx <- regexpr(sep, x, fixed = TRUE)
  ifelse(idx > 0, substr(x, 1, idx - 1), x)
}

read_tsv_gz <- function(path) {
  read.table(gzfile(path), header = TRUE, sep = "\t",
             comment.char = "", check.names = FALSE, quote = "")
}

train_seeds  <- parse_seed_range(opt$train_seeds)
predict_seed <- opt$predict_seed
first_sv_col <- opt$`first-sv-col`
gene_sep     <- opt$`gene-sep`
pattern      <- opt$pattern
min_nonzero  <- opt$`min-nonzero`

# ── load SD / SR tables ───────────────────────────────────────────────────────
sd_path <- if (!is.null(opt$sd)) opt$sd else file.path(opt$dir, "gene_SD_overlap.tsv.gz")
sr_path <- if (!is.null(opt$sr)) opt$sr else file.path(opt$dir, "gene_SR_overlap.tsv.gz")
for (p in c(sd_path, sr_path)) if (!file.exists(p)) stop("File not found: ", p)

cat("Loading SD overlap:", sd_path, "\n")
sd_tab <- read_tsv_gz(sd_path);  rownames(sd_tab) <- as.character(sd_tab[[1]])
cat("Loading SR overlap:", sr_path, "\n")
sr_tab <- read_tsv_gz(sr_path);  rownames(sr_tab) <- as.character(sr_tab[[1]])

predict_col <- as.character(predict_seed)
if (!predict_col %in% colnames(sd_tab)) stop("Column '", predict_col, "' not in SD table")
if (!predict_col %in% colnames(sr_tab)) stop("Column '", predict_col, "' not in SR table")

# ── discover permutation files ────────────────────────────────────────────────
all_files  <- list.files(opt$dir, full.names = TRUE)
seed_match <- regexec(pattern, basename(all_files), perl = TRUE)
file_seeds <- sapply(seed_match, function(m) if (length(m) >= 2) as.integer(m[2]) else NA_integer_)
valid      <- !is.na(file_seeds)
all_files  <- all_files[valid];  file_seeds <- file_seeds[valid]

needed <- c(train_seeds, predict_seed)
keep   <- file_seeds %in% needed
all_files  <- all_files[keep];  file_seeds <- file_seeds[keep]

if (length(all_files) == 0) stop("No matching tsv.gz files found in: ", opt$dir)
cat(sprintf("Found %d permutation files\n", length(all_files)))

# ── load all files ────────────────────────────────────────────────────────────
sv_list      <- list()
sv_col_names <- NULL

cat("Loading files...\n")
for (i in order(file_seeds)) {
  seed <- file_seeds[i]; f <- all_files[i]
  cat(sprintf("  seed %d: %s\n", seed, basename(f)))
  dat <- tryCatch(read_tsv_gz(f), error = function(e) { warning(e$message); NULL })
  if (is.null(dat)) next
  if (ncol(dat) < first_sv_col) {
    warning(sprintf("seed %d: only %d cols — skipping", seed, ncol(dat))); next
  }
  dat[[1]] <- strip_suffix(as.character(dat[[1]]), gene_sep)
  rownames(dat) <- dat[[1]]
  if (is.null(sv_col_names)) {
    sv_col_names <- colnames(dat)[first_sv_col:ncol(dat)]
    cat(sprintf("  -> %d SV feature columns (col %d-%d)\n",
                length(sv_col_names), first_sv_col, ncol(dat)))
  }
  sv_list[[as.character(seed)]] <- dat
}

if (is.null(sv_col_names)) stop("No files loaded.")
if (!as.character(predict_seed) %in% names(sv_list)) stop("Predict seed not loaded.")

train_seeds_avail <- intersect(as.character(train_seeds), names(sv_list))
cat(sprintf("Training seeds loaded: %d\n", length(train_seeds_avail)))

# ── common genes ──────────────────────────────────────────────────────────────
common_genes <- Reduce(intersect, c(
  list(rownames(sd_tab), rownames(sr_tab)),
  lapply(sv_list, rownames)
))
cat(sprintf("Genes common across all files: %d\n", length(common_genes)))
if (length(common_genes) == 0) stop("No common genes. Check --gene-sep.")

n_genes <- length(common_genes)
n_cols  <- length(sv_col_names)

# ── pre-build SD/SR matrices: genes x seeds (training + predict) ──────────────
# sd_mat[gene, seed]: SD overlap proportion
all_seeds_needed <- c(train_seeds_avail, predict_col)
sd_mat <- matrix(NA_real_, nrow = n_genes, ncol = length(all_seeds_needed),
                 dimnames = list(common_genes, all_seeds_needed))
sr_mat <- matrix(NA_real_, nrow = n_genes, ncol = length(all_seeds_needed),
                 dimnames = list(common_genes, all_seeds_needed))

for (s in all_seeds_needed) {
  if (s %in% colnames(sd_tab)) sd_mat[, s] <- as.numeric(sd_tab[common_genes, s])
  if (s %in% colnames(sr_tab)) sr_mat[, s] <- as.numeric(sr_tab[common_genes, s])
}

# SD/SR for prediction (seed 0)
sd0 <- sd_mat[, predict_col]
sr0 <- sr_mat[, predict_col]

# ── output matrices ───────────────────────────────────────────────────────────
pred_mat <- matrix(NA_real_, nrow = n_genes, ncol = n_cols,
                   dimnames = list(common_genes, sv_col_names))
obs_mat  <- matrix(NA_real_, nrow = n_genes, ncol = n_cols,
                   dimnames = list(common_genes, sv_col_names))

# regression summary: one row per (gene x SV_feature)
reg_list <- vector("list", n_cols)

dat0 <- sv_list[[predict_col]]

cat(sprintf("Running per-gene regressions for %d SV columns x %d genes...\n", n_cols, n_genes))

for (ci in seq_along(sv_col_names)) {
  sv_col <- sv_col_names[ci]
  if (ci %% 10 == 0) cat(sprintf("  col %d/%d: %s\n", ci, n_cols, sv_col))

  # Build sv_mat for this column: genes x train_seeds
  sv_mat_col <- matrix(NA_real_, nrow = n_genes, ncol = length(train_seeds_avail),
                       dimnames = list(common_genes, train_seeds_avail))
  for (s in train_seeds_avail) {
    dat_s <- sv_list[[s]]
    g_s   <- intersect(common_genes, rownames(dat_s))
    sv_mat_col[g_s, s] <- as.numeric(dat_s[g_s, sv_col])
  }

  # Observed seed 0
  g0 <- intersect(common_genes, rownames(dat0))
  obs_mat[g0, sv_col] <- as.numeric(dat0[g0, sv_col])

  # Per-gene regression
  gene_rows <- lapply(common_genes, function(gene) {
    sv_vec <- sv_mat_col[gene, ]   # length = n_train_seeds
    sd_vec <- sd_mat[gene, train_seeds_avail]
    sr_vec <- sr_mat[gene, train_seeds_avail]

    ok <- complete.cases(sv_vec, sd_vec, sr_vec)
    if (sum(ok) < 3) return(NULL)  # need at least 3 obs to fit intercept + 2 betas

    # Skip genes with near-zero variance in SV counts (all zeros / constant)
    if (sum(sv_vec[ok] != 0) < min_nonzero) return(NULL)

    fit <- tryCatch(lm(sv ~ sd + sr, data = data.frame(sv = sv_vec[ok], sd = sd_vec[ok], sr = sr_vec[ok])),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)

    # Predict seed 0 for this gene
    pred_val <- tryCatch(
      predict(fit, newdata = data.frame(sd = sd0[gene], sr = sr0[gene])),
      error = function(e) NA_real_
    )

    s  <- summary(fit)
    cf <- coef(s)
    list(
      pred    = as.numeric(pred_val),
      intercept = cf["(Intercept)", "Estimate"],
      beta_sd   = if ("sd" %in% rownames(cf)) cf["sd", "Estimate"]     else NA_real_,
      beta_sr   = if ("sr" %in% rownames(cf)) cf["sr", "Estimate"]     else NA_real_,
      pval_sd   = if ("sd" %in% rownames(cf)) cf["sd", "Pr(>|t|)"]     else NA_real_,
      pval_sr   = if ("sr" %in% rownames(cf)) cf["sr", "Pr(>|t|)"]     else NA_real_,
      r2        = s$r.squared,
      n         = sum(ok)
    )
  })
  names(gene_rows) <- common_genes

  # Fill prediction matrix and collect regression rows
  col_reg <- lapply(common_genes, function(gene) {
    r <- gene_rows[[gene]]
    if (is.null(r)) return(NULL)
    pred_mat[gene, sv_col] <<- r$pred
    data.frame(gene = gene, sv_feature = sv_col,
               intercept = r$intercept, beta_sd = r$beta_sd, beta_sr = r$beta_sr,
               pval_sd = r$pval_sd, pval_sr = r$pval_sr,
               r_squared = r$r2, n_obs = r$n,
               stringsAsFactors = FALSE)
  })
  reg_list[[ci]] <- do.call(rbind, col_reg[!sapply(col_reg, is.null)])
}

# ── write outputs ─────────────────────────────────────────────────────────────
to_df <- function(mat) {
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  cbind(gene = rownames(df), df)
}

pred_out <- file.path(opt$output_dir, "sv_predicted_seed0.tsv.gz")
obs_out  <- file.path(opt$output_dir, "sv_observed_seed0.tsv.gz")
reg_out  <- file.path(opt$output_dir, "regression_summary.tsv.gz")

cat("Writing", pred_out, "\n")
write.table(to_df(pred_mat), gzfile(pred_out), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Writing", obs_out, "\n")
write.table(to_df(obs_mat), gzfile(obs_out), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Writing", reg_out, "\n")
reg_tab <- do.call(rbind, reg_list[!sapply(reg_list, is.null)])
write.table(reg_tab, gzfile(reg_out), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")


suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-d", "--dir"),          type = "character", default = NULL,
              help = "Directory containing permuted tsv.gz files and SD/SR overlap tables",
              metavar = "character"),
  make_option(c("--sd"),                 type = "character", default = NULL,
              help = "gene_SD_overlap.tsv.gz path (overrides --dir default)",
              metavar = "character"),
  make_option(c("--sr"),                 type = "character", default = NULL,
              help = "gene_SR_overlap.tsv.gz path (overrides --dir default)",
              metavar = "character"),
  make_option(c("--pattern"),            type = "character",
              default = "gnomAD_SV_v3\\.vs\\.hg38_gencode_v39\\.permuted_seed(\\d+)\\.tsv\\.gz$",
              help    = "Regex with one capture group for seed number [default: gnomAD naming]",
              metavar = "character"),
  make_option(c("--gene-sep"),           type = "character", default = ".permu",
              help    = "Separator used to strip seed suffix from col 1 gene names [default: .permu]",
              metavar = "character"),
  make_option(c("--first-sv-col"),       type = "integer",   default = 37L,
              help    = "1-based index of first SV feature column [default: 37]",
              metavar = "integer"),
  make_option(c("--train-seeds"),        type = "character", default = "1:100",
              help    = "Seed range to use for training, e.g. '1:100' or '1,2,3' [default: 1:100]",
              metavar = "character"),
  make_option(c("--predict-seed"),       type = "integer",   default = 0L,
              help    = "Seed to predict (uses SD/SR col matching this seed) [default: 0]",
              metavar = "integer"),
  make_option(c("-o", "--output-dir"),   type = "character", default = ".",
              help    = "Output directory [default: current dir]",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$dir)) stop("--dir is required")

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- helper: parse seed range -----------------------------------------------
parse_seed_range <- function(s) {
  if (grepl("^\\d+:\\d+$", s)) {
    r <- as.integer(strsplit(s, ":")[[1]])
    return(r[1]:r[2])
  }
  as.integer(strsplit(s, ",")[[1]])
}

train_seeds   <- parse_seed_range(opt$train_seeds)
predict_seed  <- opt$predict_seed
first_sv_col  <- opt$`first-sv-col`
gene_sep      <- opt$`gene-sep`
pattern       <- opt$pattern

# ---- helper: strip gene suffix ----------------------------------------------
strip_suffix <- function(x, sep) {
  idx <- regexpr(sep, x, fixed = TRUE)
  ifelse(idx > 0, substr(x, 1, idx - 1), x)
}

# ---- helper: fast gz reader -------------------------------------------------
read_tsv_gz <- function(path) {
  read.table(gzfile(path), header = TRUE, sep = "\t",
             comment.char = "", check.names = FALSE, quote = "")
}

# ---- resolve SD / SR paths --------------------------------------------------
sd_path <- if (!is.null(opt$sd)) opt$sd else file.path(opt$dir, "gene_SD_overlap.tsv.gz")
sr_path <- if (!is.null(opt$sr)) opt$sr else file.path(opt$dir, "gene_SR_overlap.tsv.gz")

for (p in c(sd_path, sr_path)) {
  if (!file.exists(p)) stop("File not found: ", p)
}

# ---- load SD / SR tables ----------------------------------------------------
cat("Loading SD overlap table:", sd_path, "\n")
sd_tab <- read_tsv_gz(sd_path)
rownames(sd_tab) <- as.character(sd_tab[[1]])

cat("Loading SR overlap table:", sr_path, "\n")
sr_tab <- read_tsv_gz(sr_path)
rownames(sr_tab) <- as.character(sr_tab[[1]])

predict_col <- as.character(predict_seed)
if (!predict_col %in% colnames(sd_tab)) stop("Column '", predict_col, "' not found in SD table")
if (!predict_col %in% colnames(sr_tab)) stop("Column '", predict_col, "' not found in SR table")

# ---- discover permutation files ---------------------------------------------
all_files  <- list.files(opt$dir, full.names = TRUE)
matched    <- regmatches(basename(all_files),
                         regexpr(pattern, basename(all_files), perl = TRUE))
seed_match <- regmatches(basename(all_files),
                         regexec(pattern, basename(all_files), perl = TRUE))

file_seeds <- sapply(seed_match, function(m) if (length(m) >= 2) as.integer(m[2]) else NA_integer_)
valid      <- !is.na(file_seeds)
all_files  <- all_files[valid]
file_seeds <- file_seeds[valid]

needed_seeds <- c(train_seeds, predict_seed)
keep         <- file_seeds %in% needed_seeds
all_files    <- all_files[keep]
file_seeds   <- file_seeds[keep]

if (length(all_files) == 0) stop("No matching tsv.gz files found in: ", opt$dir)
cat(sprintf("Found %d files (seeds: %s ... %s)\n",
            length(all_files), min(file_seeds), max(file_seeds)))

# ---- load all files upfront -------------------------------------------------
# sv_list: named list by seed (character) -> data.frame with clean gene names as rownames

sv_list        <- list()
sv_col_names   <- NULL

cat("Loading permutation files...\n")

for (i in order(file_seeds)) {
  seed <- file_seeds[i]
  f    <- all_files[i]
  cat(sprintf("  seed %d: %s\n", seed, basename(f)))

  dat <- tryCatch(read_tsv_gz(f),
                  error = function(e) { warning("Failed: ", f, " -- ", e$message); NULL })
  if (is.null(dat)) next

  if (ncol(dat) < first_sv_col) {
    warning(sprintf("seed %d has only %d cols (< %d) — skipping", seed, ncol(dat), first_sv_col))
    next
  }

  # Strip seed suffix from col 1 gene names
  dat[[1]] <- strip_suffix(as.character(dat[[1]]), gene_sep)
  rownames(dat) <- dat[[1]]

  # Record SV column names from first file loaded
  if (is.null(sv_col_names)) {
    sv_col_names <- colnames(dat)[first_sv_col:ncol(dat)]
    cat(sprintf("  SV feature columns: %d (col %d to %d)\n",
                length(sv_col_names), first_sv_col, ncol(dat)))
  }

  sv_list[[as.character(seed)]] <- dat
}

if (length(sv_list) == 0) stop("No files loaded successfully.")
if (!as.character(predict_seed) %in% names(sv_list)) {
  stop("Prediction seed ", predict_seed, " not loaded.")
}

train_seeds_avail <- intersect(as.character(train_seeds), names(sv_list))
if (length(train_seeds_avail) == 0) stop("No training seeds loaded.")
cat(sprintf("Training seeds available: %d\n", length(train_seeds_avail)))

# ---- find common genes ------------------------------------------------------
common_genes <- Reduce(intersect, c(
  list(rownames(sd_tab), rownames(sr_tab)),
  lapply(sv_list, rownames)
))
cat(sprintf("Genes common across all files: %d\n", length(common_genes)))

if (length(common_genes) == 0) stop("No common genes found. Check gene name stripping (--gene-sep).")

# Pre-extract prediction SD/SR vectors (seed 0, common genes only)
sd0 <- as.numeric(sd_tab[common_genes, predict_col])
sr0 <- as.numeric(sr_tab[common_genes, predict_col])
names(sd0) <- common_genes
names(sr0) <- common_genes

# ---- regression loop --------------------------------------------------------
n_cols  <- length(sv_col_names)
pred_mat <- matrix(NA_real_, nrow = length(common_genes), ncol = n_cols,
                   dimnames = list(common_genes, sv_col_names))
obs_mat  <- matrix(NA_real_, nrow = length(common_genes), ncol = n_cols,
                   dimnames = list(common_genes, sv_col_names))

reg_rows <- vector("list", n_cols)

cat(sprintf("Running regressions for %d SV columns across %d training seeds...\n",
            n_cols, length(train_seeds_avail)))

for (ci in seq_along(sv_col_names)) {
  sv_col <- sv_col_names[ci]
  if (ci %% 10 == 0) cat(sprintf("  col %d/%d: %s\n", ci, n_cols, sv_col))

  # Build long training data.frame: gene x seed
  rows <- lapply(train_seeds_avail, function(s) {
    dat    <- sv_list[[s]]
    # only common genes that exist in this file
    g      <- intersect(common_genes, rownames(dat))
    sv_vec <- as.numeric(dat[g, sv_col])
    sd_vec <- as.numeric(sd_tab[g, s])
    sr_vec <- as.numeric(sr_tab[g, s])
    data.frame(sv = sv_vec, sd = sd_vec, sr = sr_vec, stringsAsFactors = FALSE)
  })
  train_df <- do.call(rbind, rows)
  train_df <- train_df[complete.cases(train_df), ]

  if (nrow(train_df) < 20) {
    warning(sprintf("Too few obs for '%s' (%d rows) — skipping", sv_col, nrow(train_df)))
    next
  }

  fit <- tryCatch(lm(sv ~ sd + sr, data = train_df),
                  error = function(e) { warning("lm failed for '", sv_col, "': ", e$message); NULL })
  if (is.null(fit)) next

  s <- summary(fit)
  coefs <- coef(s)

  reg_rows[[ci]] <- data.frame(
    sv_feature    = sv_col,
    intercept     = coefs["(Intercept)", "Estimate"],
    beta_sd       = if ("sd" %in% rownames(coefs)) coefs["sd", "Estimate"]         else NA_real_,
    beta_sr       = if ("sr" %in% rownames(coefs)) coefs["sr", "Estimate"]         else NA_real_,
    pval_sd       = if ("sd" %in% rownames(coefs)) coefs["sd", "Pr(>|t|)"]         else NA_real_,
    pval_sr       = if ("sr" %in% rownames(coefs)) coefs["sr", "Pr(>|t|)"]         else NA_real_,
    r_squared     = s$r.squared,
    adj_r_squared = s$adj.r.squared,
    n_obs         = nrow(train_df),
    stringsAsFactors = FALSE
  )

  # Predict for seed 0
  ok_genes   <- common_genes[complete.cases(data.frame(sd = sd0, sr = sr0))]
  newdata    <- data.frame(sd = sd0[ok_genes], sr = sr0[ok_genes])
  pred_vals  <- predict(fit, newdata = newdata)
  pred_mat[ok_genes, sv_col] <- pred_vals

  # Observed for seed 0
  dat0 <- sv_list[[as.character(predict_seed)]]
  g0   <- intersect(common_genes, rownames(dat0))
  obs_mat[g0, sv_col] <- as.numeric(dat0[g0, sv_col])
}

# ---- write outputs ----------------------------------------------------------
to_df <- function(mat) {
  df <- as.data.frame(mat, stringsAsFactors = FALSE)
  df <- cbind(gene = rownames(df), df)
  rownames(df) <- NULL
  df
}

pred_out <- file.path(opt$output_dir, "sv_predicted_seed0.tsv.gz")
obs_out  <- file.path(opt$output_dir, "sv_observed_seed0.tsv.gz")
reg_out  <- file.path(opt$output_dir, "regression_summary.tsv")

cat("Writing", pred_out, "\n")
write.table(to_df(pred_mat), gzfile(pred_out), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Writing", obs_out, "\n")
write.table(to_df(obs_mat), gzfile(obs_out), sep = "\t", quote = FALSE, row.names = FALSE)

cat("Writing", reg_out, "\n")
reg_tab <- do.call(rbind, reg_rows[!sapply(reg_rows, is.null)])
write.table(reg_tab, reg_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")
