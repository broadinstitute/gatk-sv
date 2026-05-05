#!/usr/bin/env Rscript
# Reorganize STR/VNTR disease table by Repeat Unit.
# Groups genes by Repeat Unit and lists disorders + gene count.

reorganize_str_by_repeat_unit <- function(
  input_file,
  output_file = NULL,
  sep = "\t"
) {
  # Read the file
  tbl <- read.table(input_file, header = TRUE, sep = sep, quote = "", stringsAsFactors = FALSE, check.names = FALSE)
  
  # Normalize and find columns
  col_names <- colnames(tbl)
  cat("Available columns:\n")
  print(col_names)
  
  # Find Repeat Unit column (case insensitive, whitespace-aware)
  repeat_unit_col <- grep("^Repeat\\s+Unit$", col_names, ignore.case = TRUE)
  if (length(repeat_unit_col) == 0) {
    repeat_unit_col <- grep("Repeat", col_names, ignore.case = TRUE)[1]
  }
  
  # Find Gene column
  gene_col <- grep("^Gene$", col_names, ignore.case = TRUE)[1]
  
  # Find Primary Disorder column (case insensitive)
  disorder_col <- grep("Primary.*Disorder|Disorder", col_names, ignore.case = TRUE)[1]
  
  if (is.na(repeat_unit_col) || is.na(gene_col) || is.na(disorder_col)) {
    cat("Repeat Unit col:", repeat_unit_col, "\n")
    cat("Gene col:", gene_col, "\n")
    cat("Disorder col:", disorder_col, "\n")
    stop("Could not identify required columns")
  }
  
  # Remove rows where Repeat Unit is empty
  tbl <- tbl[trimws(tbl[[repeat_unit_col]]) != "", ]
  
  # Group by Repeat Unit
  repeat_units <- unique(trimws(tbl[[repeat_unit_col]]))
  repeat_units <- repeat_units[nzchar(repeat_units)]
  repeat_units <- sort(repeat_units)
  
  result <- data.frame(
    Repeat_Unit = character(),
    Genes = character(),
    Primary_Disorders = character(),
    Gene_Count = integer(),
    stringsAsFactors = FALSE
  )
  
  for (ru in repeat_units) {
    subset_rows <- trimws(tbl[[repeat_unit_col]]) == ru
    subset_tbl <- tbl[subset_rows, ]
    
    # Aggregate genes
    genes <- trimws(subset_tbl[[gene_col]])
    genes <- genes[nzchar(genes)]
    genes_list <- paste(unique(genes), collapse = ",")
    gene_count <- length(unique(genes))
    
    # Aggregate disorders
    disorders <- trimws(subset_tbl[[disorder_col]])
    disorders <- disorders[nzchar(disorders)]
    disorders_list <- paste(unique(disorders), collapse = ";")
    
    result <- rbind(result, data.frame(
      Repeat_Unit = ru,
      Genes = genes_list,
      Primary_Disorders = disorders_list,
      Gene_Count = gene_count,
      stringsAsFactors = FALSE
    ))
  }
  
  rownames(result) <- NULL
  
  # Output to console
  print(result)
  
  # Save to file if specified
  if (!is.null(output_file)) {
    write.table(result, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("\nOutput saved to:", output_file, "\n")
  }
  
  invisible(result)
}

# Main execution
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript reorganize_str_by_repeat_unit.R <input_file> [output_file]\n")
    quit(status = 1)
  }
  
  input_file <- args[1]
  output_file <- if (length(args) > 1) args[2] else NULL
  
  if (!file.exists(input_file)) {
    stop(paste0("Input file not found: ", input_file))
  }
  
  reorganize_str_by_repeat_unit(input_file, output_file)
}
