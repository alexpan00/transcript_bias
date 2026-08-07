library(tidyverse)
library(NOISeq)
library(mgcv)

cpm <- function(df){
  df <- as.matrix(df)
  df_cpm <- t(10^6*t(df)/colSums(df))
  return(df_cpm)
}

# define functions to sample SIRVs, compute contrasts and plotting the results
evaluate_sampling_rmsd <- function(
  counts,
  metadata,
  factor_col,
  condition1,
  condition2,
  fc_reference,
  annotation,
  id_col = "id",
  n_replicates = 1e4,
  n_reads = 20000,
  seed = 1
) {
  # Check inputs
  if (!all(c(condition1, condition2) %in% unique(metadata[[factor_col]]))) {
    # return an empty list if conditions are not found
    return(list())
  }
  if (!fc_reference %in% colnames(annotation)) {
    stop("fc_reference must be a column in the annotation.")
  }

  # Apply multinomial sampling
  set.seed(seed)
  if (n_replicates < 1) {
    stop("n_replicates must be at least 1.")
  } else if (n_replicates == 1) {
    sampled_list <- list(
      apply(counts, 2, function(x){rmultinom(1, n_reads, prob = x)})
    )
  } else  {
    sampled_list <- replicate(n_replicates, {
      apply(counts, 2, function(x){rmultinom(1, n_reads, prob = x)})
    }, simplify = FALSE)
  }

  # Compute mean sampled counts per condition and replicate
  fc_list <- lapply(sampled_list, function(sampled_counts) {
    mean_counts <- data.frame(sapply(
      c(condition1, condition2), 
      function(cond) {
        rowMeans(sampled_counts[, metadata[[factor_col]] == cond, drop = FALSE])
      }
    ))
    colnames(mean_counts) <- c(condition1, condition2)
    mean_counts$FC <- log2((mean_counts[[condition1]] + 1) / (mean_counts[[condition2]] + 1))
    return(mean_counts$FC)
  })

  # Convert to matrix: rows = features, cols = replicates
  fc_matrix <- do.call(cbind, fc_list)
  mean_fc <- rowMeans(fc_matrix)
  sd_fc <- apply(fc_matrix, 1, sd)

  # Merge with annotation
  out_df <- data.frame(id = rownames(counts), mean_log2FC = mean_fc, sd_FC = sd_fc)
  final <- merge(annotation, out_df, by.x = id_col, by.y = "id")

  # Compute RMSD
  rmsd <- sqrt(mean((log2(final[[fc_reference]]) - (final$mean_log2FC))^2))

  return(list(
    results = final,
    rmsd = round(rmsd, 3)
  ))
}

plot_sampling_rmsd <- function(results_list, 
                               fc_col = "E1_E0", 
                               length_col = "Length", 
                               label_x = NULL,
                               label_y = NULL) {

  # if the list is empty, rturn empty plot
  if (length(results_list) == 0) {
    return(ggplot() + ggtitle("No data to plot"))
  }
  # Colors_constant
  COLORS <- c("#6036a4", "#9de26b", "#F5C290", "#B1B3B3")
  # Extract elements
  plot_df <- results_list$results
  rmsd <- results_list$rmsd
  
  # Generate unique log2FC levels
  unique_log2fc <- sort(unique(round(log2(plot_df[[fc_col]]), 2)))
  # Get the necessary colors
  colors <- COLORS[1:length(unique_log2fc)]

  # Default position for RMSD label
  if (is.null(label_x)) label_x <- max(plot_df[[length_col]]) * 0.8
  if (is.null(label_y)) label_y <- max(plot_df$mean_log2FC) * 0.9
  
  # Plot
  p <- ggplot(plot_df, aes_string(x = length_col, y = "mean_log2FC", 
                                  color = paste0("factor(round(log2(", fc_col, "), 2))"))) +
    geom_smooth(se =F, alpha = .5) +
    geom_point() +
    geom_hline(yintercept = unique_log2fc, color = colors, linetype = "dashed") +
    scale_color_manual("Expected log2FC", values = colors) +
    ggtitle(fc_col) +
    theme_light() +
    annotate("text", x = label_x, y = label_y, 
             label = paste("RMSD:", rmsd)) +
    labs(x = "Transcript Length", y = "Observed log2 Fold Change") +
    theme(legend.position = "bottom")
  
  return(p)
}

args <- commandArgs(trailingOnly = TRUE)
long_obj <- args[1]
sirv_info <- args[2]
output_prefix <- args[3]

# --- Prepare SIRV info to account for the contrasts --- #

# Read the SIRV info
sirv_info <- read.csv(sirv_info)


# conver franctions from str to number
sirv_info$E1 <- sapply(sirv_info$E1, function(x){eval(parse(text=x))})
sirv_info$E2 <- sapply(sirv_info$E2, function(x){eval(parse(text=x))})

# compute ratios
# sirv_info$E1_E0 <- sirv_info$E1/sirv_info$E0
# sirv_info$E2_E1 <- sirv_info$E2/sirv_info$E1
# sirv_info$E2_E0 <- sirv_info$E2/sirv_info$E0
# sirv_info$E0_E0 <- sirv_info$E0/sirv_info$E0


# Read long reads data and normalize using TMM
mydata <- readRDS(long_obj)
mydata_cpm <- exprs(mydata)

# keep only the true sirvs present in the info
mydata_cpm <- cpm(mydata_cpm[rownames(mydata_cpm) %in% sirv_info$id,])

# There are big gaps in the number of SIRVs reads per sample. In order to to the
# downstream analysis, I will sample the same number of reads per sample, that is
# the minimum number of reads in a sample
min_reads <- round(min(colSums(mydata_cpm)))

# define contrasts
main_factor <- colnames(pData(mydata))[2] 
conds <- levels(pData(mydata)[,2])
sirv_contrasts <- combn(conds, 2)

# reduce metadata
sirv_col <- if ("SIRV" %in% colnames(pData(mydata))) "SIRV" else (if (ncol(pData(mydata)) >= 3) colnames(pData(mydata))[3] else main_factor)
sirv_cond <- pData(mydata)[, unique(c(main_factor, sirv_col)), drop = FALSE] %>% distinct()
# --- Evaluate contrasts --- #
# E1 vs E0
for (i in 1:ncol(sirv_contrasts)){
  cond1 <- sirv_contrasts[1, i]
  cond2 <- sirv_contrasts[2, i]
  cond_contrast <- paste(cond1, cond2, sep = "_")
  
  sirv_cond1 <- if (sirv_col %in% colnames(sirv_cond)) as.character(sirv_cond[sirv_cond[,main_factor] == cond1, sirv_col]) else "Unknown"
  sirv_cond2 <- if (sirv_col %in% colnames(sirv_cond)) as.character(sirv_cond[sirv_cond[,main_factor] == cond2, sirv_col]) else "Unknown"
  
  sirv_contrast <- paste(sirv_cond1, sirv_cond2, sep = "_")
  if (!sirv_contrast %in% colnames(sirv_info)){
    sirv_info[sirv_contrast] <- sirv_info[, sirv_cond1]/sirv_info[,sirv_cond2]
  }
    p1 <- tryCatch({
      res_tmm <- evaluate_sampling_rmsd(mydata_cpm,
                             pData(mydata),
                             main_factor,
                             cond1,
                             cond2,
                             sirv_contrast,
                             sirv_info,
                             n_replicates = 1,
                             n_reads = min_reads)
      plot_sampling_rmsd(res_tmm, fc_col = sirv_contrast)
    }, error = function(e){
      message(paste("Error in contrast", cond_contrast, ":", e$message))
      ggplot() + ggtitle(paste("Error in contrast", cond_contrast))
    })
    
    ggsave(p1, filename = paste0(output_prefix, "_", cond_contrast ,".png"), height = 5, width = 7)
 
}
