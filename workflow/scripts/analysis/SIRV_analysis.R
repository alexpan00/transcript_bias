library(NOISeq)
library(Biobase)
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
sirv_obj <- args[2]
sirv_info <- args[3]
output_prefix <- args[4]

#show files paths
print("Input files:")
print(long_obj)
print(sirv_obj)
print(sirv_info)


# Read long and short reads data
long_obj <- readRDS(long_obj)
sirv_obj <- readRDS(sirv_obj)

# Read the SIRV info
sirv_info <- read.csv(sirv_info)

# Plot correlation
# Select main condition from pData
cond <- if (ncol(pData(long_obj)) >= 2) colnames(pData(long_obj))[2] else colnames(pData(long_obj))[1]
secondary_cond <- if (ncol(pData(long_obj)) >= 3) colnames(pData(long_obj))[3] else "SIRV"
if (is.na(secondary_cond) || secondary_cond == "SIRV"){
  secondary_cond <- cond
}
## common transcripts
## Full Join
tryCatch({
  p <- mycor.plot(cor.dat(long_obj, sirv_obj, cond, full.join = F, norm = T))
  output_cor_common <- paste0(output_prefix, "_sirv_cor_common.png")
  ggsave(filename = output_cor_common, plot = p, height = 6.5, width = 6.5)
}, error = function(e) {
  warning("Error generating correlation plot: ", e$message, ". Creating empty plot.")
  p <- ggplot() + theme_void()
  output_cor_common <- paste0(output_prefix, "_sirv_cor_common.png")
  ggsave(filename = output_cor_common, plot = p, height = 6.5, width = 6.5)
})

# Save cor.data for summary plot, in this case the correlation is computed by sample
# First we filter the long object to keep only the SIRV transcripts, then we can
# do a full.join with the sirv object that includes the 69 sirvs
sample_col <- colnames(pData(long_obj))[1]
matching_ids <- intersect(rownames(fData(long_obj)), sirv_info$id)

sirv_feature_data <- fData(sirv_obj)
has_length <- "Length" %in% colnames(sirv_feature_data)

if (has_length) {
  transcript_lengths <- suppressWarnings(as.numeric(as.character(sirv_feature_data[, "Length"])))
  names(transcript_lengths) <- rownames(sirv_feature_data)
  valid_lengths <- transcript_lengths[!is.na(transcript_lengths)]

  if (length(valid_lengths) > 0L) {
    length_breaks <- unique(as.numeric(quantile(valid_lengths, probs = seq(0, 1, 0.25), na.rm = TRUE)))
  } else {
    length_breaks <- numeric()
  }
}

if (length(matching_ids) == 0 || nrow(sirv_obj) == 0) {
  warning("No SIRV transcripts detected in long reads data or empty SIRV object. Creating empty summary.")
  l_res <- list()
  if (nrow(sirv_obj) > 0) {
    log_cpm_sirv <- log(cpm(exprs(sirv_obj)) + 1)
    for (sample in colnames(exprs(sirv_obj))) {
      sample_ids <- pData(long_obj)[, sample_col]
      idx <- match(sample, sample_ids)
      condition <- if (!is.na(idx)) as.character(pData(long_obj)[idx, cond]) else "Unknown"
      l_res[[sample]] <- data.frame(Measured = 0, Expected = log_cpm_sirv[,sample], Sample = sample, Condition = condition, Transcript = rownames(exprs(sirv_obj)))
    }
  } else {
    l_res[["empty"]] <- data.frame(Measured = numeric(0), Expected = numeric(0), Sample = character(0), Condition = character(0), Transcript = character(0))
  }
  summary_df <- do.call(rbind, l_res)
  output_summary <- paste0(output_prefix, "_sirv_summary.rds")
  saveRDS(summary_df, output_summary)

  p_empty <- ggplot() + theme_void() + labs(title = "No SIRV transcripts detected")

  output_detection <- paste0(output_prefix, "_sirv_detection.png")
  ggsave(filename = output_detection, plot = p_empty, height = 5, width = 10)

  output_detection_10 <- paste0(output_prefix, "_sirv_detection_10_reads.png")
  ggsave(filename = output_detection_10, plot = p_empty, height = 5, width = 10)

  output_box <- paste0(output_prefix, "_sirv_boxplot.png")
  ggsave(filename = output_box, plot = p_empty, height = 4.5, width = 6.5)

  output_len <- paste0(output_prefix, "_sirv_vs_len.png")
  ggsave(filename = output_len, plot = p_empty, height = 5, width = 9)

  output_len_gt <- paste0(output_prefix, "_sirv_vs_len_gt.png")
  ggsave(filename = output_len_gt, plot = p_empty, height = 5, width = 9)

  quit(save = "no", status = 0)
}

sirv_long <- long_obj[matching_ids,]
sirv_vs_gt <- cor.dat(sirv_long, sirv_obj, full.join = T, norm = T) # if condition is not specified is computed by sample

l_res <- list()
for (sample in sirv_vs_gt$Samples){
  sample_df <- sirv_vs_gt$data2plot[[sample]]
  colnames(sample_df) <- c("Measured", "Expected")
  sample_df$Transcript <- rownames(sample_df)
  sample_df$Sample <- sample
  # Find the matching condition for this sample in the sirv_object
  sample_ids <- pData(sirv_long)[, sample_col]
  idx <- match(sample, sample_ids)
  if (is.na(idx)) {
    stop("Sample '", sample, "' not found in pData(sirv_long)[, ", sample_col, "].")
  }
  if (sum(sample_ids == sample) > 1L) {
    stop("Sample '", sample, "' appears multiple times in pData(sirv_long)[, ", sample_col, "].")
  }
  condition <- as.character(pData(sirv_long)[idx, cond])
  sample_df$Condition <- condition

  if (has_length) {
    sample_df$Length <- transcript_lengths[sample_df$Transcript]
    sample_df$LengthQuantile <- assign_length_quantile(sample_df$Length, length_breaks)
  }

  l_res[[sample]] <- sample_df
}

summary_df <- do.call(rbind, l_res)
output_summary <- paste0(output_prefix, "_sirv_summary.rds")
saveRDS(summary_df, output_summary)

## Get count matrix for the SIRVs
long_exprs <- data.frame(exprs(long_obj), check.names = FALSE)
sirv_exprs <- long_exprs[sirv_info$id,]

rownames(sirv_exprs) <- sirv_info$id
sirv_exprs[is.na(sirv_exprs)] <- 0

# Get a plot of detected SIRVs
p <- sirv_exprs %>% 
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
output_detection <- paste0(output_prefix, "_sirv_detection.png")
ggsave(filename = output_detection, plot = p, height = 5, width = 10)

p <- sirv_exprs %>% 
  rownames_to_column("transcript_id") %>% 
  pivot_longer(cols = !transcript_id,
               names_to = c("samples"),
               values_to = "counts") %>% 
  mutate(counts = if_else(counts > 10, 1, 0)) %>% 
  ggplot(aes(x= transcript_id, y = samples, fill=counts)) +
  geom_tile(color="black") + 
  scale_fill_gradient(low = "#F44336", high = "#4CAF50") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
output_detection <- paste0(output_prefix, "_sirv_detection_10_reads.png")
ggsave(filename = output_detection, plot = p, height = 5, width = 10)

# CPM-rescale within the SIRV subset, unless the values are already on a log
# scale (ratio_correction), where rescaling by a column sum is meaningless.
# CPM is idempotent, so no guard is needed for already-CPM input.
if (is_log_scale(sirv_exprs)) {
  sirv_exprs <- data.frame(sirv_exprs, check.names = FALSE)
} else {
  sums_exprs <- colSums(sirv_exprs)
  sums_exprs[sums_exprs == 0] <- 1
  sirv_exprs <- data.frame(t(10^6 * t(sirv_exprs) / sums_exprs), check.names = FALSE)
}
sirv_exprs$id <- rownames(sirv_exprs)

# Convert to long format and add metadata
sirv_exprs_long <- pivot_longer(sirv_exprs, cols = !id, names_to = "sample", values_to = "CPM")
sirv_exprs_long <- merge(sirv_exprs_long, pData(long_obj))
sirv_exprs_long$id_set <- paste(sirv_exprs_long$id, sirv_exprs_long$SIRV, sep = "_")

# Prepare SIRV info
sirv_info$E0 <- as.character(sirv_info$E0) # is all 1s so it is read as number
sirv_info_long <- pivot_longer(sirv_info, cols = starts_with("E"), names_to = "SIRV", values_to = "abundance")
sirv_info_long$id_set <- paste(sirv_info_long$id, sirv_info_long$SIRV, sep = "_")

# Combine SIRV info and expression
sirv_exprs_long <- merge(sirv_info_long, sirv_exprs_long)
sirv_exprs_long$Pool_set <- paste(sirv_exprs_long[,secondary_cond], sirv_exprs_long$SIRV, sep = "_")
abundance <- unique(sirv_exprs_long$abundance)
abundance <- abundance[order(sapply(abundance, function(x){eval(parse(text=x))}))]
sirv_exprs_long$abundance <- factor(sirv_exprs_long$abundance, levels=abundance)

sirv_exprs_long_mean <- sirv_exprs_long %>% 
  group_by(id, SIRV, Length, GC, abundance,.data[[secondary_cond]]) %>%
  summarise("CPM" = mean(CPM))

p <- ggplot(sirv_exprs_long_mean, aes(x=abundance, y=CPM)) +
  geom_boxplot() +
  geom_point(aes(size=Length), alpha = 0.25) +
  facet_grid(scales = "free_x", rows=vars(.data[[secondary_cond]]), cols = vars(SIRV), space = "free") +
  scale_radius(range = c(0.5,3)) +
  theme_light() +
  theme(legend.position = "bottom")

output_box <- paste0(output_prefix, "_sirv_boxplot.png")
ggsave(filename = output_box, plot = p, height = 4.5, width = 6.5)

my_colors <- RColorBrewer::brewer.pal(n = 9, name = "BuPu")[3:8]  # skip lightest
p <- ggplot(sirv_exprs_long_mean, aes(x=Length, y = log10(CPM+1))) +
  geom_smooth(aes(color = abundance), se = F) +
  geom_point(aes(color=abundance), alpha = .75) +
  facet_grid(rows=vars(.data[[secondary_cond]]), cols = vars(SIRV)) +
  theme_light() +
  scale_color_manual(values=my_colors) +
  ylim(c(0,log10(max(sirv_exprs_long_mean$CPM)) + 0.25))

output_len <- paste0(output_prefix, "_sirv_vs_len.png")
ggsave(filename = output_len, plot = p, height = 5, width = 9)

## SIRV ground truth as long
sirv_gt <- exprs(sirv_obj)
if (is_log_scale(sirv_gt)) {
  sirv_gt <- data.frame(sirv_gt, check.names = FALSE)
} else {
  sums_gt <- colSums(sirv_gt)
  sums_gt[sums_gt == 0] <- 1
  sirv_gt <- data.frame(t(10^6 * t(sirv_gt) / sums_gt), check.names = FALSE)
}
sirv_gt$id <- rownames(sirv_gt)

# Convert to long format and add metadata
sirv_gt_long <- pivot_longer(sirv_gt, cols = !id, names_to = "sample", values_to = "Expected")

# Merge gt and expression
sirv_gt_long$id_sample <- paste(sirv_gt_long$id, sirv_gt_long$sample, sep = "_")
sirv_exprs_long$id_sample <- paste(sirv_exprs_long$id, sirv_exprs_long$sample, sep = "_")
sirv_exprs_long <- merge(sirv_exprs_long, sirv_gt_long)

sirv_exprs_long_mean <- sirv_exprs_long %>% 
  group_by(id, SIRV, Length, GC, abundance,.data[[secondary_cond]]) %>%
  summarise("CPM_EXP_diff" = mean(CPM-Expected))
# Diference between expected and measured considering length
p <- ggplot(sirv_exprs_long_mean, aes(x=Length, y = CPM_EXP_diff)) +
  geom_smooth(aes(color = abundance), se = F) +
  geom_point(aes(color=abundance), alpha = .75) +
  facet_grid(rows=vars(.data[[secondary_cond]]), cols = vars(SIRV)) +
  ylab("CPM - Expected") +
  theme_light() +
  scale_color_manual(values=my_colors)

output_len_gt <- paste0(output_prefix, "_sirv_vs_len_gt.png")
ggsave(filename = output_len_gt, plot = p, height = 5, width = 9)
