library(NOISeq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

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

assign_length_quantile <- function(lengths, breaks) {
  if (length(breaks) == 0L) {
    return(rep(NA_character_, length(lengths)))
  }

  if (length(breaks) == 1L) {
    return(ifelse(is.na(lengths), NA_character_, as.character(breaks)))
  }

  # Create labels as ranges (e.g., "0-1500")
  labels <- paste(breaks[-length(breaks)], breaks[-1], sep = "-")
  as.character(cut(lengths, breaks = breaks, include.lowest = TRUE, labels = labels))
}

args <- commandArgs(trailingOnly = TRUE)
long_obj <- args[1]
ercc_obj <- args[2]
output_prefix <- args[3]

#show files paths
print("Input files:")
print(long_obj)
print(ercc_obj)


# Read long and short reads data
long_obj <- readRDS(long_obj)
ercc_obj <- readRDS(ercc_obj)


# Plot correlation
# Select main condition from pData
cond <- if (ncol(pData(long_obj)) >= 2) colnames(pData(long_obj))[2] else colnames(pData(long_obj))[1]
secondary_cond <- if (ncol(pData(long_obj)) >= 3) colnames(pData(long_obj))[3] else cond

## common transcripts
tryCatch({
  p <- mycor.plot(cor.dat(long_obj, ercc_obj, cond, full.join = FALSE))
  output_cor_common <- paste0(output_prefix, "_ercc_cor_common.png")
  ggsave(filename = output_cor_common, plot = p, height = 6.5, width = 6.5)
}, error = function(e) {
  warning("Error generating correlation plot: ", e$message, ". Creating empty plot.")
  p <- ggplot() + theme_void()
  output_cor_common <- paste0(output_prefix, "_ercc_cor_common.png")
  ggsave(filename = output_cor_common, plot = p, height = 6.5, width = 6.5)
})

# Save cor.data for summary plot, in this case the correlation is computed by sample
# First we filter the long object to keep only the ERCC transcripts, then we can
# do a full.join with the ercc object that includes the  erccs
sample_col <- colnames(pData(long_obj))[1]
ercc_feature_data <- fData(ercc_obj)
has_length <- "Length" %in% colnames(ercc_feature_data)

if (has_length) {
  transcript_lengths <- suppressWarnings(as.numeric(as.character(ercc_feature_data[, "Length"])))
  names(transcript_lengths) <- rownames(ercc_feature_data)
  valid_lengths <- transcript_lengths[!is.na(transcript_lengths)]

  if (length(valid_lengths) > 0L) {
    length_breaks <- unique(as.numeric(quantile(valid_lengths, probs = seq(0, 1, 0.25), na.rm = TRUE)))
  } else {
    length_breaks <- numeric()
  }
}

matching_ids <- intersect(rownames(fData(long_obj)), rownames(ercc_obj))

# in some cases no ercc transcripts are detected, so we need to check if matching_ids or ercc_obj is empty
if (length(matching_ids) == 0 || nrow(ercc_obj) == 0) {
  warning("No ERCC transcripts detected in the long reads data or empty ERCC object. Creating empty summary.")
  l_res <- list()
  if (nrow(ercc_obj) > 0) {
    log_cpm_ercc <- log(cpm(exprs(ercc_obj)) + 1)
    for (sample in colnames(exprs(ercc_obj))) {
      sample_ids <- pData(long_obj)[, sample_col]
      idx <- match(sample, sample_ids)
      condition <- if (!is.na(idx)) as.character(pData(long_obj)[idx, cond]) else "Unknown"
      l_res[[sample]] <- data.frame(Measured = 0, Expected = log_cpm_ercc[,sample], Sample = sample, Condition = condition, Transcript = rownames(exprs(ercc_obj)))
    }
  } else {
    l_res[["empty"]] <- data.frame(Measured = numeric(0), Expected = numeric(0), Sample = character(0), Condition = character(0), Transcript = character(0))
  }
} else {
  ercc_long <- long_obj[matching_ids,]
  ercc_vs_gt <- cor.dat(ercc_long, ercc_obj, full.join = T, norm = T) # if condition is not specified is computed by sample
  l_res <- list()
  for (sample in ercc_vs_gt$Samples){
    sample_df <- ercc_vs_gt$data2plot[[sample]]
    colnames(sample_df) <- c("Measured", "Expected")
    sample_df$Transcript <- rownames(sample_df)
    sample_df$Sample <- sample
    # Find the matching condition for this sample in the ercc_object
    sample_ids <- pData(long_obj)[, sample_col]
    idx <- match(sample, sample_ids)
    condition <- if (!is.na(idx)) as.character(pData(long_obj)[idx, cond]) else "Unknown"
    sample_df$Condition <- condition

    if (has_length) {
      sample_df$Length <- transcript_lengths[sample_df$Transcript]
      sample_df$LengthQuantile <- assign_length_quantile(sample_df$Length, length_breaks)
    }

    l_res[[sample]] <- sample_df
  }
}
summary_df <- do.call(rbind, l_res)
output_summary <- paste0(output_prefix, "_ercc_summary.rds")
saveRDS(summary_df, output_summary)

## Get count matrix for the all the ERCCs
long_exprs <- data.frame(exprs(long_obj))
common_erccs <- intersect(rownames(ercc_obj), rownames(long_exprs))
if (length(common_erccs) > 0) {
  ercc_exprs <- long_exprs[common_erccs, , drop = FALSE]
  ercc_exprs[is.na(ercc_exprs)] <- 0
  p <- ercc_exprs %>% 
    rownames_to_column("transcript_id") %>% 
    pivot_longer(cols = !transcript_id,
                 names_to = c("samples"),
                 values_to = "counts") %>% 
    mutate(counts = if_else(counts > 0, 1, 0)) %>% 
    ggplot(aes(x= transcript_id, y = samples, fill=counts)) +
    geom_tile(color="black") + 
    scale_fill_gradient(low = "#F44336", high = "#4CAF50") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "none")
} else {
  p <- ggplot() + theme_void() + labs(title = "No ERCC transcripts detected")
}
output_detection <- paste0(output_prefix, "_ercc_detection.png")
ggsave(filename = output_detection, plot = p, height = 5, width = 10)
