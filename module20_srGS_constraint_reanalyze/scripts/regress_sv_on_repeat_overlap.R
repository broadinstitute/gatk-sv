#!/usr/bin/env Rscript
# regress_sv_on_repeat_overlap.R
#
# For each SV feature column (37+) across a series of permuted RData files:
#   - Extract SV counts per gene per seed
#   - Match each seed to the corresponding SD/SR overlap column (by seed number)
#   - Fit: SV_count ~ SD_overlap + SR_overlap  using all seeds EXCEPT seed 0
#   - Predict the expected SV_count for seed 0 using SD/SR from column "0"
#   - Output: predicted vs observed table + regression summary

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("-d", "--rdata-dir"),  type = "character", default = NULL,
              help = "Directory containing the permuted .rData files", metavar = "character"),
  make_option(c("--sd-overlap"),       type = "character", default = NULL,
              help = "gene_SD_overlap.tsv (gene x seed columns)", metavar = "character"),
  make_option(c("--sr-overlap"),       type = "character", default = NULL,
              help = "gene_SR_overlap.tsv (gene x seed columns)", metavar = "character"),
  make_option(c("-o", "--output-dir"), type = "character", default = ".",
              help = "Output directory [default: current dir]", metavar = "character"),
  make_option(c("--first-sv-col"),     type = "integer",   default = 37L,
              help = "1-based index of first SV feature column in RData [default: 37]", metavar = "integer"),
  make_option(c("--rdata-pattern"),    type = "character",
              default = "gnomAD_SV_v3\\.vs\\.hg38_gencode_v39\\.permuted_seed(\\d+)\\.rData",
              help = "Regex with one capture group for seed number [default matches gnomAD naming]",
              metavar = "character"),
  make_option(c("--rdata-obj"),        type = "character", default = "gene.data.reanno.permu",
              help = "Name of object inside each .rData file [default: gene.data.reanno.permu]",
              metavar = "character")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("rdata_dir", "sd_overlap", "sr_overlap")
missing  <- required[sapply(required, function(x) is.null(opt[[x]]) || opt[[x]] == "")]
if (length(missing) > 0) stop(paste("Missing required args:", paste(paste0("--", gsub("_", "-", missing)), collapse = ", ")))

dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Discover RData files and extract seed numbers ----------------------

rdata_files <- list.files(opt$rdata_dir, pattern = "\\.rData$", full.names = TRUE)
if (length(rdata_files) == 0) stop("No .rData files found in: ", opt$rdata_dir)

seed_pattern <- opt$rdata_pattern
seeds <- regmatches(basename(rdata_files), regexpr(seed_pattern, basename(rdata_files), perl = TRUE))
seed_nums <- as.integer(gsub(seed_pattern, "\\1", seeds, perl = TRUE))

valid <- !is.na(seed_nums)
rdata_files <- rdata_files[valid]
seed_nums   <- seed_nums[valid]

if (length(rdata_files) == 0) stop("No RData files matched the seed pattern: ", seed_pattern)

cat(sprintf("Found %d RData files (seeds: %s)\n", length(rdata_files),
            paste(sort(seed_nums), collapse = ", ")))

# ---- 2. Load SD / SR overlap tables ----------------------------------------

cat("Loading SD overlap table...\n")
sd_tab <- read.table(opt$sd_overlap, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
rownames(sd_tab) <- sd_tab$gene

cat("Loading SR overlap table...\n")
sr_tab <- read.table(opt$sr_overlap, header = TRUE, sep = "\t", check.names = FALSE, comment.char = "")
rownames(sr_tab) <- sr_tab$gene

# Verify seed 0 column exists
if (!"0" %in% colnames(sd_tab)) stop("Column '0' not found in gene_SD_overlap.tsv")
if (!"0" %in% colnames(sr_tab)) stop("Column '0' not found in gene_SR_overlap.tsv")

# Common genes across SD and SR
common_genes_repeat <- intersect(rownames(sd_tab), rownames(sr_tab))
cat(sprintf("Genes in SD: %d, SR: %d, intersection: %d\n",
            nrow(sd_tab), nrow(sr_tab), length(common_genes_repeat)))

# ---- 3. Load all RData files and collect SV feature matrices ---------------

first_sv_col <- opt$first_sv_col
obj_name     <- opt$rdata_obj

cat("Loading RData files...\n")

# sv_data: named list by seed -> data.frame (genes x SV features)
sv_data <- list()

sv_col_names <- NULL  # will be set on first load

for (i in seq_along(rdata_files)) {
  seed <- seed_nums[i]
  cat(sprintf("  seed %d: %s\n", seed, basename(rdata_files[i])))

  env <- new.env(parent = emptyenv())
  load(rdata_files[i], envir = env)

  if (!exists(obj_name, envir = env)) {
    warning(sprintf("Object '%s' not found in %s — skipping", obj_name, basename(rdata_files[i])))
    next
  }

  dat <- get(obj_name, envir = env)

  if (ncol(dat) < first_sv_col) {
    warning(sprintf("File %s has only %d columns (expected >= %d) — skipping",
                    basename(rdata_files[i]), ncol(dat), first_sv_col))
    next
  }

  # Store gene column name (col 1 expected to be 'gene')
  gene_col <- colnames(dat)[1]
  sv_cols  <- colnames(dat)[first_sv_col:ncol(dat)]

  if (is.null(sv_col_names)) {
    sv_col_names <- sv_cols
  } else if (!identical(sv_col_names, sv_cols)) {
    warning(sprintf("SV column names differ in seed %d — using intersection", seed))
    sv_col_names <- intersect(sv_col_names, sv_cols)
  }

  rownames(dat) <- as.character(dat[[gene_col]])
  sv_data[[as.character(seed)]] <- dat
}

if (length(sv_data) == 0) stop("No valid RData files loaded.")
if (!"0" %in% names(sv_data)) stop("Seed 0 RData not found or failed to load.")

cat(sprintf("Loaded %d seeds, %d SV feature columns\n", length(sv_data), length(sv_col_names)))

# ---- 4. Build combined dataset and run regressions -------------------------
# For each SV feature column:
#   Training rows: genes * seeds (seed != 0), columns: sv_count, sd, sr
#   Prediction: genes in seed 0 with sd[,0] and sr[,0]

# Common genes across all loaded seeds + repeat overlap tables
all_genes <- Reduce(intersect, c(
  list(common_genes_repeat),
  lapply(sv_data, function(d) rownames(d))
))

cat(sprintf("Genes available for regression: %d\n", length(all_genes)))

# Seeds for training (all except 0)
train_seeds <- setdiff(names(sv_data), "0")
if (length(train_seeds) == 0) stop("Need at least one seed other than 0 for regression training.")

# Pre-extract SD and SR for seed 0 (prediction)
sd0 <- as.numeric(sd_tab[all_genes, "0"])
sr0 <- as.numeric(sr_tab[all_genes, "0"])

# Output containers
pred_df    <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
obs_df     <- data.frame(gene = all_genes, stringsAsFactors = FALSE)
reg_summary <- list()

cat(sprintf("Running regressions for %d SV feature columns across %d training seeds...\n",
            length(sv_col_names), length(train_seeds)))

for (sv_col in sv_col_names) {
  # Build training data frame
  train_rows <- do.call(rbind, lapply(train_seeds, function(s) {
    sv_vec <- as.numeric(sv_data[[s]][all_genes, sv_col])
    sd_vec <- as.numeric(sd_tab[all_genes, s])
    sr_vec <- as.numeric(sr_tab[all_genes, s])
    data.frame(
      gene  = all_genes,
      seed  = as.integer(s),
      sv    = sv_vec,
      sd    = sd_vec,
      sr    = sr_vec,
      stringsAsFactors = FALSE
    )
  }))

  # Remove rows with NA
  train_rows <- train_rows[complete.cases(train_rows), ]

  if (nrow(train_rows) < 10) {
    warning(sprintf("Too few complete observations for column '%s' — skipping", sv_col))
    next
  }

  # Fit model
  fit <- tryCatch(
    lm(sv ~ sd + sr, data = train_rows),
    error = function(e) { warning(sprintf("lm failed for '%s': %s", sv_col, e$message)); NULL }
  )
  if (is.null(fit)) next

  reg_summary[[sv_col]] <- summary(fit)

  # Predict for seed 0
  newdata <- data.frame(sd = sd0, sr = sr0)
  newdata <- newdata[complete.cases(newdata), , drop = FALSE]
  pred_genes <- all_genes[complete.cases(data.frame(sd = sd0, sr = sr0))]

  pred_vals <- predict(fit, newdata = newdata)

  # Observed for seed 0
  obs_vals <- as.numeric(sv_data[["0"]][pred_genes, sv_col])

  # Store
  tmp_pred          <- setNames(rep(NA_real_, length(all_genes)), all_genes)
  tmp_pred[pred_genes] <- pred_vals
  pred_df[[sv_col]] <- tmp_pred

  tmp_obs          <- setNames(rep(NA_real_, length(all_genes)), all_genes)
  tmp_obs[pred_genes] <- obs_vals
  obs_df[[sv_col]] <- tmp_obs
}

# ---- 5. Write outputs ------------------------------------------------------

pred_out <- file.path(opt$output_dir, "sv_predicted_seed0.tsv")
obs_out  <- file.path(opt$output_dir, "sv_observed_seed0.tsv")
reg_out  <- file.path(opt$output_dir, "regression_summary.tsv")

cat(sprintf("Writing predicted values -> %s\n", pred_out))
write.table(pred_df, pred_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Writing observed values  -> %s\n", obs_out))
write.table(obs_df,  obs_out,  sep = "\t", quote = FALSE, row.names = FALSE)

# Summary table: one row per SV feature
cat(sprintf("Writing regression summary -> %s\n", reg_out))
reg_rows <- lapply(names(reg_summary), function(sv_col) {
  s     <- reg_summary[[sv_col]]
  coefs <- coef(s)
  data.frame(
    sv_feature      = sv_col,
    intercept       = coefs["(Intercept)", "Estimate"],
    beta_sd         = coefs["sd",          "Estimate"],
    beta_sr         = coefs["sr",          "Estimate"],
    pval_sd         = coefs["sd",          "Pr(>|t|)"],
    pval_sr         = coefs["sr",          "Pr(>|t|)"],
    r_squared       = s$r.squared,
    adj_r_squared   = s$adj.r.squared,
    n_obs           = nrow(s$residuals),
    stringsAsFactors = FALSE
  )
})
reg_tab <- do.call(rbind, reg_rows)
write.table(reg_tab, reg_out, sep = "\t", quote = FALSE, row.names = FALSE)

cat("Done.\n")
