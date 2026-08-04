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
  s_dir <- get_script_dir()
  target <- file.path(s_dir, script_name)
  if (file.exists(target)) {
    source(target)
  } else if (file.exists(file.path("workflow/scripts", script_name))) {
    source(file.path("workflow/scripts", script_name))
  } else if (file.exists(file.path("scripts", script_name))) {
    source(file.path("scripts", script_name))
  } else {
    source(script_name)
  }
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
counts <- read.table(counts, header = TRUE, row.names = 1, sep = ",")

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