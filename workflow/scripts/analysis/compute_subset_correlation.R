
get_script_dir <- function() {
  if (exists("snakemake") && !is.null(snakemake@script)) {
    return(dirname(snakemake@script))
  }
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(sub("^--file=", "", file_arg[1])))
  }
  possible_dirs <- c("workflow/scripts", "scripts", ".")
  for (d in possible_dirs) {
    if (dir.exists(d)) return(d)
  }
  return(".")
}

source_script <- function(script_name) {
  categories <- c(".", "plotting", "analysis", "normalization", "preprocessing", "quantification", "summary", "utils")
  s_dir <- get_script_dir()
  for (cat in categories) {
    cands <- c(
      file.path(s_dir, cat, script_name),
      file.path(s_dir, script_name),
      file.path("workflow/scripts", cat, script_name),
      file.path("scripts", cat, script_name),
      file.path(cat, script_name)
    )
    for (cand in cands) {
      if (file.exists(cand)) {
        source(cand)
        return(invisible(TRUE))
      }
    }
  }
  source(script_name)
}

source_script("corplot.R")

assign_length_group <- function(lengths) {
  # Define static breaks and labels
  # Breaks: 0 to 1500, 1501 to 3000, 3001 to 4500, and > 4500
  breaks <- c(0, 1500, 3000, 4500, Inf)
  labels <- c("0-1.5k", "1.5k-3k", "3k-4.5k", ">4.5k")
  
  as.character(cut(lengths, breaks = breaks, include.lowest = TRUE, labels = labels))
}

args <- commandArgs(trailingOnly = TRUE)

# Help/Usage check
if (length(args) < 4) {
  stop("Usage: Rscript compute_subset_correlation.R <eset_path> <short_obj> <tsv_path> <output_path>")
}

eset_path <- args[1]
short_obj <- args[2]
tsv_path <- args[3]
output_path <- args[4]

# Load required libraries
# using suppressPackageStartupMessages to keep logs clean
suppressPackageStartupMessages(library(Biobase))

# 1. Read the eset object and the TSV file
if (!file.exists(eset_path)) {
  stop(paste("Eset file not found:", eset_path))
}
if (!file.exists(tsv_path)) {
  stop(paste("TSV file not found:", tsv_path))
}

message(paste("Loading eset from:", eset_path))
eset <- readRDS(eset_path)

message(paste("Loading short reads object from:", short_obj))
short_obj <- readRDS(short_obj)

message(paste("Loading ID mapping from:", tsv_path))
mapping <- read.table(tsv_path, header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "")

# 2. Extract transcript names and clean version
# Get feature names (transcript names)
eset_ids <- featureNames(eset)
message(paste("Total features in eset:", length(eset_ids)))

# Clean the version of the transcript (remove anything after a point)
clean_ids <- sub("\\..*", "", eset_ids)

# 3. Identify which of the identifiers is matching the ids in the eset
match_col <- NULL
max_overlap <- 0
best_pct <- 0

message("Checking columns for matches...")
for (col in colnames(mapping)) {
  # Convert to character to ensure comparing strings
  col_vals <- as.character(mapping[[col]])
  col_vals <- col_vals[!is.na(col_vals) & col_vals != ""] # Ignore NAs and empty strings
  
  if (length(col_vals) == 0) next
  
  # Check intersection
  intersection <- intersect(clean_ids, col_vals)
  overlap_count <- length(intersection)
  
  # Calculate percentage based on the smaller of the two sets primarily, 
  # but here we care about finding the column that best explains the eset IDs.
  # Let's count how many eset IDs are found in this column.
  overlap_pct <- (overlap_count / length(clean_ids)) * 100
  
  message(paste("  Column:", col, "- Overlap:", overlap_count, "(", round(overlap_pct, 2), "%)"))
  
  if (overlap_count > max_overlap) {
    max_overlap <- overlap_count
    match_col <- col
    best_pct <- overlap_pct
  }
}

if (is.null(match_col) || max_overlap == 0) {
  warning("No matching column found with any overlap. Eset will be empty or script should fail.")
  # Depending on strictness, we might want to stop here.
  # But technically, the subset would just be empty.
} else {
  message(paste("Best matching column identified:", match_col))
}

# 4. Subset the eset to contain only the ids present in the tsv file
# The requirement is: "subset the eset to contain only the ids present in the tsv file"
# We match based on the cleaned ID
if (!is.null(match_col)) {
  valid_clean_ids <- unique(mapping[[match_col]])
  keep_mask <- clean_ids %in% valid_clean_ids
  
  eset_subset <- eset[keep_mask, ]
  
  message(paste("Original expression set dimensions:", paste(dim(eset), collapse = "x")))
  message(paste("Subset expression set dimensions:", paste(dim(eset_subset), collapse = "x")))
  # Clean the short reads ids and use the same column to subset the short reads object if needed
  clean_ids_short <- sub("\\..*", "", rownames(short_obj))
  keep_mask_short <- clean_ids_short %in% valid_clean_ids
  short_obj_subset <- short_obj[keep_mask_short, ]
 
  # 5 Compute correalation and save the results
  # Select main condition from pData
  cond = colnames(pData(eset_subset))[2]
  sample_col <- colnames(pData(eset_subset))[1]
  lr_vs_sr <- cor.dat(eset_subset, short_obj_subset, sample_col, full.join = T, norm = T)

  feature_data <- fData(eset_subset)
  has_length <- "Length" %in% colnames(feature_data)

  if (has_length) {
    transcript_lengths <- suppressWarnings(as.numeric(as.character(feature_data[, "Length"])))
    names(transcript_lengths) <- rownames(feature_data)
  }

  l_res <- list()
  for (sample in lr_vs_sr$Samples){
    sample_df <- lr_vs_sr$data2plot[[sample]]
    colnames(sample_df) <- c("Measured", "Expected")
    sample_df$Transcript <- rownames(sample_df)
    sample_df$Sample <- sample
    # Find the matching condition for this sample in the eset
    sample_ids <- pData(eset_subset)[, sample_col]
    idx <- match(sample, sample_ids)
    if (is.na(idx)) {
        stop("Sample '", sample, "' not found in pData(eset_subset)[, ", sample_col, "].")
    }
    if (sum(sample_ids == sample) > 1L) {
        stop("Sample '", sample, "' appears multiple times in pData(eset_subset)[, ", sample_col, "].")
    }
    condition <- as.character(pData(eset_subset)[idx, cond])
    sample_df$Condition <- condition

    if (has_length) {
      sample_df$Length <- transcript_lengths[sample_df$Transcript]
      sample_df$LengthQuantile <- assign_length_group(sample_df$Length)
    }

    l_res[[sample]] <- sample_df
}

summary_df <- do.call(rbind, l_res)
saveRDS(summary_df, output_path)

} else {
  stop("Could not identify a matching ID column to perform subsetting.")
}
