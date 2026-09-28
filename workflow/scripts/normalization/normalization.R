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



# Read long reads data and apply the requested normalization
mydata <- readRDS(long_obj)
# Compute library sizes from raw counts to use in post-filtering
lib_sizes <- colSums(exprs(mydata))

if (norm_method == "TPM"){
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

  lengths <- fData(mydata)$Length
  names(lengths) <- rownames(fData(mydata))
  gc_content <- fData(mydata)$GC
  names(gc_content) <- rownames(fData(mydata))
  # cqn's own default subindex is which(rowMeans(counts) > 50). Long-read
  # counts often fall below that, and an empty or near-empty subindex makes
  # cqn fail with "'from' must be a finite number" or a singular design
  # matrix, so fall back to the most expressed transcripts instead.
  sel_idx <- which(rowMeans(exprs(mydata)) > 50)
  if (length(sel_idx) < 100) {
    sel_idx <- select_top_n(exprs(mydata), n = 3500)
    message(sprintf(
      "cqn: only %d transcripts have a mean count above 50; fitting on the %d most expressed instead.",
      sum(rowMeans(exprs(mydata)) > 50), length(sel_idx)))
  }
  cqn_obj <- cqn(exprs(mydata), lengths=lengths,
                x=gc_content, verbose=TRUE, 
                lengthMethod="smooth", subindex=sel_idx)
  # CQn returns log2 normalized values, so we need to transform back
  norm_counts <- 2**(cqn_obj$y + cqn_obj$offset)
} else if (norm_method == "ratio_counts") {
  cpm_offset <- 1
  norm_counts <- ratio_correction(exprs(mydata),
                                 cpm_offset,
                                 return.counts = TRUE)
} else {
  stop("Normalization method not recognized. Please use TPM, EDA, CPM, ratio_correction, ratio_counts, or cqn.")
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