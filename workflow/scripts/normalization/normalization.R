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
library(NOISeq)


args <- commandArgs(trailingOnly = TRUE)
long_obj <- args[1]
output_prefix <- args[2]
norm_method <- args[3]
post_filtering <- if (length(args) >= 4) as.numeric(args[4]) else 0
params_file <- if (length(args) >= 5) args[5] else NULL



# Read long reads data and normalize using TMM
mydata <- readRDS(long_obj)
# Compute library sizes from raw counts to use in post-filtering
lib_sizes <- colSums(exprs(mydata))

if (norm_method == "TMM"){
  norm_counts <- tmm(exprs(mydata))
} else if (norm_method == "TPM"){
  lengths <- fData(mydata)$Length
  names(lengths) <- rownames(fData(mydata))
  norm_counts <- calculate_tpm(exprs(mydata), lengths)
} else if (norm_method == "EDA"){
  lengths <- fData(mydata)$Length
  names(lengths) <- rownames(fData(mydata))
  norm_counts <- gcLoess(exprs(mydata), log10(lengths))
} else if (norm_method == "CPM"){
  norm_counts <- cpm(exprs(mydata))
} else if (norm_method == "ratio_correction"){
  cpm_offset <- 1
  norm_counts <- ratio_correction(exprs(mydata),
                                 cpm_offset)
} else if (norm_method == "cqn"){
  library(cqn)
  CQN_n_trans <- 3500 # Number of transcripts to use for CQN normalization

  lengths <- fData(mydata)$Length
  names(lengths) <- rownames(fData(mydata))
  gc_content <- fData(mydata)$GC
  names(gc_content) <- rownames(fData(mydata))
  sel_idx <- select_top_n(exprs(mydata), n = CQN_n_trans)
  cqn_obj <- cqn(exprs(mydata), lengths=lengths,
                x=gc_content, verbose=TRUE, 
                lengthMethod="smooth", subindex=sel_idx)
  # CQn returns log2 normalized values, so we need to transform back
  norm_counts <- 2**(cqn_obj$y + cqn_obj$offset)
} else if (norm_method == "read_density"){
  # Check that params_file is provided
  if (is.null(params_file)) {
    stop("For read_density normalization, you must provide a params_file as the 4th argument.\n",
         "The file should be a CSV/TSV with columns: sample, mean, sd")
  }
  
  # Read parameters file
  if (grepl("\\.csv$", params_file)) {
    sample_params <- read.csv(params_file, row.names = 1)
  } else {
    sample_params <- read.table(params_file, header = TRUE, row.names = 1, sep = "\t")
  }
  
  # Validate parameters file
  required_cols <- c("mean", "sd")
  if (!all(required_cols %in% colnames(sample_params))) {
    stop("params_file must contain columns: mean, sd")
  }
  
  # Validate that all samples in data have parameters
  missing_samples <- setdiff(colnames(exprs(mydata)), rownames(sample_params))
  if (length(missing_samples) > 0) {
    stop("Missing parameters for samples: ", paste(missing_samples, collapse = ", "))
  }
  
  # Apply optimal epsilon correction
  cat("Applying optimal epsilon correction...\n")
  result <- optimal_epsilon_correction(
    mydata,
    sample_params,
    simulation_dist = "norm",
    remove_outliers = FALSE,
    epsilon_interval = c(0.0001, 1),
    min_count = 0.0
  )
  
  norm_counts <- result$corrected_counts
  
  # Save optimal epsilon values for reference
  write.csv(data.frame(sample = names(result$optimal_epsilons), 
                      epsilon = result$optimal_epsilons),
           paste0(output_prefix, "_optimal_epsilons.csv"),
           row.names = FALSE)
  cat("Optimal epsilon values saved to:", paste0(output_prefix, "_optimal_epsilons.csv\n"))
} else if (norm_method == "ratio_counts") {
  cpm_offset <- 1
  norm_counts <- ratio_correction(exprs(mydata),
                                 cpm_offset,
                                 return.counts = TRUE)
} else {
  stop("Normalization method not recognized. Please use TMM, TPM, EDA, CPM, ratio_correction, ratio_counts, cqn, or read_density.")
}
# Compute ratio of lib sizes before and after normalization for reference
lib_sizes_after <- colSums(norm_counts)
lib_size_ratios <- lib_sizes_after / lib_sizes
print(lib_size_ratios)

# Adjust the post_filtering threshold based on the library size ratios
adjusted_post_filtering <- post_filtering * lib_size_ratios

# Save normalized counts
exprs(mydata) <- as.matrix(norm_counts)
print(adjusted_post_filtering)
if (post_filtering > 0) {
  myfactor <- if (ncol(pData(mydata)) >= 2) colnames(pData(mydata))[2] else colnames(pData(mydata))[1]
  # Filter using normalized expression but original library sizes for threshold calculation
  condition_counts <- sapply(unique(pData(mydata)[,myfactor]), function(cond){
    cond_mask <- pData(mydata)[,myfactor] == cond
    thresh <- adjusted_post_filtering[cond_mask]
    apply(norm_counts[, cond_mask, drop=FALSE], 1, 
    function(x){all(x >= thresh)})
  })
  condition_counts <- as.matrix(condition_counts)
  sel <- apply(condition_counts, 1, any)

  # 1. Extract the raw pieces
  exprs_matrix <- exprs(mydata)[sel, , drop=FALSE]
  pheno_data   <- pData(mydata)
  feat_data    <- fData(mydata)[sel, , drop=FALSE]

  # 2. Rebuild a fresh, clean ExpressionSet
  mydata_filtered <- ExpressionSet(
      assayData   = exprs_matrix,
      phenoData   = AnnotatedDataFrame(pheno_data),
      featureData = AnnotatedDataFrame(feat_data)
  )
  mydata <- mydata_filtered
  print(paste("Applied post-filtering with threshold:", adjusted_post_filtering))
  print(table(sel))
}
saveRDS(mydata, paste0(output_prefix, "_NOIseq.rds"))