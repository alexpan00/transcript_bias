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
                             return.counts = FALSE,
                             return.coutns = FALSE) {
  if (missing(return.counts) && !missing(return.coutns)) {
    return.counts <- return.coutns
  }
  # Compute cpm
  cpm_counts <- cpm(counts)

  # Compute the median per transcript across all samples
  trans_median <- if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::rowMedians(as.matrix(cpm_counts))
  } else {
    apply(as.matrix(cpm_counts), 1, median)
  }

  # Convert to counts scaling by a constant factor
  if (return.coutns) {
    sample_ratios <- (cpm_counts) / matrix(trans_median + cpm_offset, nrow = nrow(cpm_counts), ncol = ncol(cpm_counts))
    sample_counts <- cpm(sample_ratios)
    return(sample_counts)
  }

  # Compute log and apply offset to avoid log(0), then center by transcript median
  log_cpm_counts <- as.matrix(cpm_counts)

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
  order(mean_exp, decreasing = TRUE)[seq_len(n)]
}

# Helper function: correct counts with optimal epsilon
correct_counts <- function(seq_res_no_zero, reads_mu, reads_sig,
                           simulation_dist = "norm",
                           remove_outliers = TRUE,
                           w_power = 1,
                           ref_len = NULL,
                           epsilon_fraction = 0.01) {
  # Step 1: Initial relative frequencies
  seq_res_no_zero$freq0 <- seq_res_no_zero$counts / sum(seq_res_no_zero$counts)
  
  if (remove_outliers){
    outlier_range <- IQR(seq_res_no_zero$eff_length) * 1.5
    min_value <- quantile(seq_res_no_zero$eff_length, probs = 0.25) - outlier_range
    max_value <- quantile(seq_res_no_zero$eff_length, probs = 0.75) + outlier_range
    seq_res_no_zero$eff_length <- ifelse(seq_res_no_zero$eff_length < min_value, min_value, seq_res_no_zero$eff_length)
    seq_res_no_zero$eff_length <- ifelse(seq_res_no_zero$eff_length > max_value, max_value, seq_res_no_zero$eff_length)
  }
  
  # Step 2: Weight by length
  if (simulation_dist == "cauchy") {
    w <- dcauchy(seq_res_no_zero$eff_length, location = reads_mu, scale = reads_sig)
  } else if (simulation_dist == "norm") {
    w <- dnorm(seq_res_no_zero$eff_length, mean = reads_mu, sd = reads_sig)
  } else if (simulation_dist == "weibull"){
    w <- dweibull(seq_res_no_zero$eff_length, shape = reads_mu, scale = reads_sig)
  }
  
  # Step 3: Apply Epsilon Floor
  epsilon <- max(w, na.rm = TRUE) * epsilon_fraction
  w_safe <- w + epsilon 
  seq_res_no_zero$w <- w_safe
  
  # Step 4: Apply w_power using the safe denominator
  seq_res_no_zero$freqw <- (seq_res_no_zero$freq0 / (w_safe^w_power)) /
    sum(seq_res_no_zero$freq0 / (w_safe^w_power))
  seq_res_no_zero$counts_corrected <- seq_res_no_zero$freqw * sum(seq_res_no_zero$counts)
  
  return(seq_res_no_zero)
}

# Helper function: create NOISeq object from corrected data
create_noiseq_object <- function(seq_res) {
  # Get transcript length
  mylength <- seq_res$eff_length
  names(mylength) <- seq_res$gene_id
  
  # Get counts columns
  counts <- seq_res[, c("gene_id", "counts_corrected"), drop = FALSE]
  if (any(rownames(counts) != counts$gene_id)){
    counts <- tibble::column_to_rownames(counts, var = "gene_id")
  } else{
    counts <- counts[, -which(colnames(counts) == "gene_id"), drop = FALSE]
  }

  factors <- data.frame(sample = "counts_corrected")
  
  # Create NOISeq object
  noiseq_data <- NOISeq::readData(data = counts, length = mylength, 
                                  factors = factors)
  return(noiseq_data)
}

# Helper function: evaluate loss for epsilon optimization
evaluate_loss_epsilon <- function(epsilon_fraction, df1, reads_mu, reads_sig, 
                                 simulation_dist = "norm", remove_outliers = FALSE) {
  # Apply correction with the given epsilon_fraction
  corrected <- correct_counts(df1, 
                               reads_mu = reads_mu, 
                               reads_sig = reads_sig, 
                               simulation_dist = simulation_dist, 
                               remove_outliers = remove_outliers,
                               epsilon_fraction = epsilon_fraction)
  
  # Create NOISeq object and evaluate length bias
  suppressMessages({
    noi_obj <- create_noiseq_object(corrected)
    dat_length <- bias.dat(noi_obj, factor = NULL, norm = FALSE)
    # Get R² from the regression model for counts_corrected
    r2 <- summary(dat_length$RegressionModels$counts_corrected)$r.squared
  })
  
  return(r2)
}

# Main function: optimal epsilon length correction
optimal_epsilon_correction <- function(noiseq_obj, 
                                      sample_params,
                                      simulation_dist = "norm",
                                      remove_outliers = FALSE,
                                      epsilon_interval = c(0.0001, 1),
                                      min_count = 0.0) {
  # noiseq_obj: NOISeq object with counts and lengths
  # sample_params: data frame with rownames as sample names, columns: mean, sd
  # simulation_dist: distribution type ("norm", "cauchy", "weibull")
  # remove_outliers: whether to remove length outliers
  # epsilon_interval: optimization interval for epsilon_fraction
  # min_count: minimum count threshold to include genes
  
  # Extract data from NOISeq object
  counts_matrix <- exprs(noiseq_obj)
  lengths_vector <- noiseq_obj@featureData@data$Length
  gene_ids <- rownames(counts_matrix)
  sample_names <- colnames(counts_matrix)
  
  # Initialize output matrix
  corrected_counts <- matrix(0, nrow = nrow(counts_matrix), ncol = ncol(counts_matrix))
  rownames(corrected_counts) <- gene_ids
  colnames(corrected_counts) <- sample_names
  
  # Store optimal epsilon values
  optimal_epsilons <- numeric(length(sample_names))
  names(optimal_epsilons) <- sample_names
  
  # Process each sample
  for (i in seq_along(sample_names)) {
    sample_name <- sample_names[i]
    cat("Processing sample:", sample_name, "\n")
    
    # Get parameters for this sample
    if (!(sample_name %in% rownames(sample_params))) {
      stop(paste("Sample", sample_name, "not found in sample_params"))
    }
    reads_mu <- sample_params[sample_name, "mean"]
    reads_sig <- sample_params[sample_name, "sd"]
    
    # Prepare data frame for this sample (only non-zero genes)
    sample_data <- data.frame(
      gene_id = gene_ids,
      counts = counts_matrix[, i],
      eff_length = lengths_vector
    )
    sample_data_nonzero <- sample_data[sample_data$counts > min_count, ]
    
    if (nrow(sample_data_nonzero) == 0) {
      warning(paste("No genes above threshold for sample", sample_name))
      corrected_counts[, i] <- counts_matrix[, i]
      next
    }
    
    # Optimize epsilon_fraction
    cat("  Optimizing epsilon_fraction...\n")
    opt_result <- optimize(
      f = evaluate_loss_epsilon,
      interval = epsilon_interval,
      df1 = sample_data_nonzero,
      reads_mu = reads_mu,
      reads_sig = reads_sig,
      simulation_dist = simulation_dist,
      remove_outliers = remove_outliers
    )
    
    best_epsilon <- opt_result$minimum
    optimal_epsilons[i] <- best_epsilon
    cat("  Optimal epsilon_fraction:", best_epsilon, "(R²:", opt_result$objective, ")\n")
    
    # Apply correction with optimal epsilon
    corrected_data <- correct_counts(
      sample_data_nonzero,
      reads_mu = reads_mu,
      reads_sig = reads_sig,
      simulation_dist = simulation_dist,
      remove_outliers = remove_outliers,
      epsilon_fraction = best_epsilon
    )
    
    # Map corrected counts back to full matrix
    corrected_counts[corrected_data$gene_id, i] <- corrected_data$counts_corrected
  }
  
  # Return results
  return(list(
    corrected_counts = corrected_counts,
    optimal_epsilons = optimal_epsilons
  ))
}

# This is taken from edgeR
filterByExprs <- function (y, design = NULL, group = NULL, lib.size = NULL, min.count = 10, 
    min.total.count = 15, large.n = 10, min.prop = 0.7, ...) 
{
    y <- as.matrix(y)
    if (mode(y) != "numeric") 
        stop("y is not a numeric matrix")
    if (is.null(lib.size)) 
        lib.size <- colSums(y)
    if (is.null(group)) {
        if (is.null(design)) {
            message("No group or design set. Assuming all samples belong to one group.")
            MinSampleSize <- ncol(y)
        }
        else {
            h <- hat(design)
            MinSampleSize <- 1/max(h)
        }
    }
    else {
        group <- as.factor(group)
        n <- tabulate(group)
        MinSampleSize <- min(n[n > 0L])
    }
    if (MinSampleSize > large.n) 
        MinSampleSize <- large.n + (MinSampleSize - large.n) * 
            min.prop
    MedianLibSize <- median(lib.size)
    CPM.Cutoff <- min.count/MedianLibSize * 1e+06
    CPM <- cpm(y)
    tol <- 1e-14
    keep.CPM <- rowSums(CPM >= CPM.Cutoff) >= (MinSampleSize - 
        tol)
    keep.TotalCount <- (rowSums(y) >= min.total.count - tol)
    keep.CPM & keep.TotalCount
}
