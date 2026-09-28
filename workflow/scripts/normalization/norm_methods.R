library(NOISeq)
if (requireNamespace("matrixStats", quietly = TRUE)) {
  library(matrixStats)
}

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

source_script("lenplot.R")

cpm <- function(df){
  df <- as.matrix(df)
  sums <- colSums(df)
  sums[sums == 0] <- 1
  df_cpm <- t(10^6*t(df)/sums)
  return(df_cpm)
}

calculate_tpm <- function(counts, lengths) {
  # counts: a matrix or data frame with genes in rows and samples in columns
  # lengths: a numeric vector of transcript lengths (in base pairs), same order as rows of counts
  
  # Step 1: Convert lengths to kilobases
  lengths_kb <- lengths / 1000
  
  # Step 2: Calculate RPK (Reads Per Kilobase)
  rpk <- sweep(counts, 1, lengths_kb, FUN = "/")
  
  # Step 3: Calculate scaling factors (sum of RPK per sample)
  scaling_factors <- colSums(rpk)
  scaling_factors[scaling_factors == 0] <- 1
  
  # Step 4: Calculate TPM
  tpm <- sweep(rpk, 2, scaling_factors, FUN = "/") * 1e6
  
  return(tpm)
}

ratio_correction <- function(counts,
                             cpm_offset,
                             return.counts = FALSE) {
  # Compute cpm
  cpm_counts <- cpm(counts)

  # Compute the median per transcript across all samples
  trans_median <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMedians(as.matrix(cpm_counts))
  } else {
    apply(as.matrix(cpm_counts), 1, median)
  }

  # Convert to counts scaling by a constant factor
  if (return.counts) {
    sample_ratios <- (cpm_counts) / matrix(trans_median + cpm_offset, nrow = nrow(cpm_counts), ncol = ncol(cpm_counts))
    sample_counts <- cpm(sample_ratios)
    return(sample_counts)
  }

  # Compute log and apply offset to avoid log(0), then center by transcript median
  log_cpm_counts <- log2(cpm_counts + cpm_offset) - matrix(log2(trans_median + cpm_offset), nrow = length(trans_median), ncol = ncol(cpm_counts))

  return(log_cpm_counts)
}

gcLoess <- function(counts,gc, degree= 2, span=0.75) {
  ff <- function(y,x, degree, span) {
    xx <- x[(y>0)&(y<=quantile(y,probs=0.99))]
    yy <- log(y[(y>0)&(y<=quantile(y,probs=0.99))])
    
    l <- loess(yy~xx, degree = degree, span = span)
    y.fit <- predict(l,newdata=x)
    names(y.fit) <- names(y)
    y.fit[is.na(y.fit)] <- 0
    
    retval <- y/exp(y.fit-median(yy))
    return(retval)
  }
  apply(counts,2,ff,x=gc, degree, span)
}

# select top n percent expressed genes
select_top_n <- function(data, n = 3500){
  mean_exp <- rowMeans(data)
  order(mean_exp, decreasing = TRUE)[seq_len(min(n, length(mean_exp)))]
}
