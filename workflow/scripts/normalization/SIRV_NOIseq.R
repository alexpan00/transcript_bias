library(NOISeq)


create_noiseq_object <- function(counts, factors, sirv_info){
  # Get transcript length
  mylength <- sirv_info$Length
  names(mylength) <- rownames(sirv_info)
  
  # Get GC% content
  mygc <- sirv_info$GC # %GC content
  names(mygc) <- rownames(sirv_info)
  # sort the colnames in counts
  counts <- counts[, sort(colnames(counts)), drop = FALSE]
  
  # order the samples in the way they are ordered in the counts matrix
  factors <- factors[match(colnames(counts), factors$sample), , drop = FALSE]
  factors[] <- lapply(factors, as.factor)
  # Create NOISeq object
  noiseq_data <- readData(data = counts,
                          factors = factors,
                          length = mylength,
                          gc = mygc)
  return(noiseq_data)
}


args <- commandArgs(trailingOnly = TRUE)
counts <- args[1]
factors <- args[2]
sirv_info <- args[3]
output <- args[4]

# Read tsv counts file
counts <- read.table(counts, header = TRUE, row.names = 1, sep = "\t")

# Remove rows with all zeros
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

# Read factors
factors <- read.table(factors, header = TRUE, sep = ",")
if ("Pool" %in% colnames(factors) && "SIRV" %in% colnames(factors)) {
  factors$Pool_Set <- paste(factors$Pool, factors$SIRV, sep = "_")
} else if ("condition" %in% colnames(factors) && "SIRV" %in% colnames(factors)) {
  factors$Pool_Set <- paste(factors$condition, factors$SIRV, sep = "_")
} else if (ncol(factors) >= 2) {
  factors$Pool_Set <- paste(factors[, 2], factors[, min(3, ncol(factors))], sep = "_")
} else {
  factors$Pool_Set <- factors[, 1]
}

sirv_info <- read.csv(sirv_info, row.names = 1)

mydata <- create_noiseq_object(counts, factors, sirv_info)

saveRDS(mydata, file = output)