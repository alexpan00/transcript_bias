library(NOISeq)
library(ggplot2)

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

source_script("norm_methods.R")
source_script("corplot.R")

args <- commandArgs(trailingOnly = TRUE)
long_obj <- args[1]
short_obj <- args[2]
output_prefix <- args[3]

# Read long and short reads data
long_obj <- readRDS(long_obj)
short_obj <- readRDS(short_obj)

assign_length_quantile <- function(lengths, breaks) {
  if (length(breaks) == 0L) {
    return(rep(NA_character_, length(lengths)))
  }

  if (length(breaks) == 1L) {
    return(ifelse(is.na(lengths), NA_character_, "Q1"))
  }

  labels <- paste0("Q", seq_len(length(breaks) - 1L))
  as.character(cut(lengths, breaks = breaks, include.lowest = TRUE, labels = labels))
}

assign_length_group <- function(lengths) {
  # Define static breaks and labels
  # Breaks: 0 to 1500, 1501 to 3000, 3001 to 4500, and > 4500
  breaks <- c(0, 1500, 3000, 4500, Inf)
  labels <- c("0-1.5k", "1.5k-3k", "3k-4.5k", ">4.5k")
  
  as.character(cut(lengths, breaks = breaks, include.lowest = TRUE, labels = labels))
}

# Select main condition from pData
cond = colnames(pData(long_obj))[2]

## Compare long vs short
# Do it by sample rather than by condition
sample_col <- colnames(pData(long_obj))[1]
## Full Join -> True
p <- mycor.plot(cor.dat(long_obj, short_obj, sample_col, full.join = T, norm = T))
output_cor_full <- paste0(output_prefix, "_cor_all.png")
ggsave(filename = output_cor_full, plot = p, height = 6.5, width = 6.5)

## common transcripts
## Full Join -> False
p <- mycor.plot(cor.dat(long_obj, short_obj, sample_col, full.join = F, norm = T))
output_cor_common <- paste0(output_prefix, "_cor_common.png")
ggsave(filename = output_cor_common, plot = p, height = 6.5, width = 6.5)


lr_vs_sr <- cor.dat(long_obj, short_obj, sample_col, full.join = F, norm = T)
feature_data <- fData(long_obj)
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
  # Find the matching condition for this sample in the long_obj
  sample_ids <- pData(long_obj)[, sample_col]
  idx <- match(sample, sample_ids)
  if (is.na(idx)) {
    stop("Sample '", sample, "' not found in pData(long_obj)[, ", sample_col, "].")
  }
  if (sum(sample_ids == sample) > 1L) {
    stop("Sample '", sample, "' appears multiple times in pData(long_obj)[, ", sample_col, "].")
  }
  condition <- as.character(pData(long_obj)[idx, cond])
  sample_df$Condition <- condition

  if (has_length) {
    sample_df$Length <- transcript_lengths[sample_df$Transcript]
    sample_df$LengthQuantile <- assign_length_group(sample_df$Length)
  }

  l_res[[sample]] <- sample_df
}

summary_df <- do.call(rbind, l_res)
output_summary <- paste0(output_prefix, "_sr_summary.rds")
saveRDS(summary_df, output_summary)