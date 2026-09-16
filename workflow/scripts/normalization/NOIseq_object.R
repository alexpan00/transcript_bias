library(NOISeq)
library(tidyverse)

xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")

create_noiseq_object <- function(counts, factors, gclength, sq_categories, min_count_condition = 0){
  # Get transcript length
  mylength <- gclength$length
  names(mylength) <- rownames(gclength)
  
  # Get GC% content
  mygc <- gclength$GC # %GC content
  names(mygc) <- rownames(gclength)
  
  # Get SQ categories
  sq_categories$category <- factor(sq_categories$category, xaxislevelsF1, labels=xaxislabelsF1)
  
  # sort the colnames in counts
  counts <- counts[,sort(colnames(counts)), drop=FALSE]
  
  # order the samples in the way they are ordered in the counts matrix
  factors <- factors[match(colnames(counts), factors$sample),, drop=F]
  factors[] <- lapply(factors, as.factor)
  
  # grouping cond
  condition <- colnames(factors)[2]

  # Filter out transcripts that do not have at least min_count_condition counts
  # in at least all the samples of one condition
  if (min_count_condition > 0){
    condition_counts <- sapply(unique(factors[,condition]), function(cond){
      apply(counts[, factors[,condition] == cond, drop=FALSE], 1, 
      function(x){all(x >= min_count_condition)})
    })
    condition_counts <- as.matrix(condition_counts)
    keep_transcripts <- apply(condition_counts, 1, any)
    # removed vs retained transcripts
    print("Transcripts removed with counts lower than threshold -> False, transcripts kept -> True:")
    print(table(keep_transcripts))
    counts <- counts[keep_transcripts, , drop=F]
  }
  # Create NOISeq object
  noiseq_data <- readData(data = counts, length = mylength, 
                          gc = mygc, factors = factors,
                          biotype = sq_categories)
  fData(noiseq_data)$Biotype <- factor(fData(noiseq_data)$Biotype, levels = xaxislabelsF1)
  return(noiseq_data)
}



args <- commandArgs(trailingOnly = TRUE)
counts <- args[1]
factors <- args[2]
gclength <- args[3]
sq_categories <- args[4]
output_prefix <- args[5]
report_factors <- unlist(strsplit(args[7], ","))
sample_2_basenames <- args[8]
min_count_condition <- as.numeric(args[9])

# Read tsv counts file
counts <- read.table(counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)

# Remove rows with all zeros
counts <- counts[rowSums(counts) > 0,, drop=F]

# make sure that proper column names are used. If colnames not in sample_2_basenames
# sample column than try to find them in the basenames columns and assign the sample names
sample_2_basenames <- read.table(sample_2_basenames, header = TRUE, sep = "\t")
if (!all(colnames(counts) %in% sample_2_basenames$sample)){
  if (all(colnames(counts) %in% sample_2_basenames$bam_basename)){
    colnames(counts) <- sample_2_basenames$sample[match(colnames(counts), sample_2_basenames$bam_basename)]
  } else {
    stop("Column names in counts file do not match sample names or basenames in sample_2_basenames file")
  }
}

# Read factors
factors <- read.table(factors, header = TRUE, sep = ",")

for (factor in report_factors) {
  if (!factor %in% colnames(factors)) {
    constituent_factors <- unlist(strsplit(factor, "_"))
    if (length(constituent_factors) > 1 && all(constituent_factors %in% colnames(factors))) {
      factors[[factor]] <- apply(factors[, constituent_factors, drop = FALSE], 1, paste, collapse = "_")
    }
  }
}

# Read gclength
gclength <- read.table(gclength, header = F, row.names = 1, sep = "\t", col.names = c("transcript_id", "length", "GC"))

# Filter gclength
gclength <- gclength[rownames(counts),]

head(gclength)

sq_categories <- read.table(sq_categories, header = T, row.names = 1, sep = "\t", col.names = c("transcript_id", "category"))

if (nrow(sq_categories) == 0 || length(unique(sq_categories$category)) == 1 && unique(sq_categories$category) == "full-splice_match"){
  sq_categories <- data.frame(category = rep("full-splice_match", nrow(counts)), row.names = rownames(counts))
} else {
  sq_categories <- sq_categories[rownames(counts), , drop = FALSE]
  sq_categories$category[is.na(sq_categories$category)] <- "full-splice_match"
}
# Create NOISeq object
mydata <- create_noiseq_object(counts, factors, gclength, sq_categories, min_count_condition)

# Save the NOIseq object
saveRDS(mydata, file = paste0(output_prefix, "_NOIseq.rds"))

