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

create_noiseq_object <- function(counts, factors){ 
  # sort the colnames in counts
  counts <- counts[,sort(colnames(counts)), drop=FALSE]
  
  # order the samples in the way they are ordered in the counts matrix
  factors <- factors[match(colnames(counts), factors$sample),, drop=FALSE]
  factors[] <- lapply(factors, as.factor)
  # Create NOISeq object
  noiseq_data <- readData(data = counts, 
                          factors = factors)
  return(noiseq_data)
}


args <- commandArgs(trailingOnly = TRUE)
counts <- args[1]
factors <- args[2]
output <- args[3]

# Read factors
factors <- read.table(factors, header = TRUE, sep = ",")

# Read tsv counts file
counts <- read.table(counts, header = TRUE, row.names = 1, sep = ",", check.names = FALSE)

# counts only have one column, repeat as many times as rows in factors
if (ncol(counts) == 1){
  counts <- counts[, rep(1, nrow(factors)), drop = FALSE]
  colnames(counts) <- factors$sample
}
# Remove rows with all zeros and calculate CPM
counts <- counts[rowSums(counts) > 0, , drop=FALSE]
counts <- cpm(counts)


mydata <- create_noiseq_object(counts, factors)

saveRDS(mydata, file = output)