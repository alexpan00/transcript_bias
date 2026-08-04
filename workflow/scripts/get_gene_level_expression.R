# From a Eset containing transcript-level expression, and a file with 2
# columns, and no header (transcript_id, gene_id) get gene-level expression
# by adding up the transcript-level expression for all transcripts belonging to the same gene.

library(NOISeq)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
eset_file <- args[1]
tx2gene_file <- args[2]
output_file <- args[3]

# Read the Eset and the tx2gene file
eset <- readRDS(eset_file)
tx2gene <- read.table(tx2gene_file,sep = "\t", header = FALSE, stringsAsFactors = FALSE, col.names = c("transcript_id", "gene_id"))

# Extract transcript lengths from fData of the input Eset
tx_lengths <- data.frame(
  transcript_id = rownames(fData(eset)),
  length = fData(eset)$Length
)

# Get the expression matrix from the Eset
exprs_matrix <- exprs(eset)
# Add gene_id to the expression matrix
exprs_df <- as.data.frame(exprs_matrix)
exprs_df$transcript_id <- rownames(exprs_df)
exprs_with_gene <- left_join(exprs_df, tx2gene, by = "transcript_id")

# Sum expression values for transcripts belonging to the same gene
gene_exprs <- exprs_with_gene %>%
  group_by(gene_id) %>%
  summarise(across(-transcript_id, sum, na.rm = TRUE))

# Convert back to a matrix and set gene_id as rownames
gene_exprs_matrix <- as.matrix(gene_exprs[,-1])
rownames(gene_exprs_matrix) <- gene_exprs$gene_id

# Compute median transcript length per gene
gene_lengths <- tx2gene %>%
  left_join(tx_lengths, by = "transcript_id") %>%
  group_by(gene_id) %>%
  summarise(median_length = median(length, na.rm = TRUE))

# Align gene_lengths to the row order of gene_exprs_matrix
gene_lengths_aligned <- gene_lengths[match(rownames(gene_exprs_matrix), gene_lengths$gene_id), ]

# Conver the gene-level expression matrix back to an Eset
gene_eset <- new("ExpressionSet", exprs = gene_exprs_matrix)

# add pData
pData(gene_eset) <- pData(eset)

# Add median transcript length as featureData
fdata <- data.frame(
  gene_id = gene_lengths_aligned$gene_id,
  Length = gene_lengths_aligned$median_length,
  row.names = rownames(gene_exprs_matrix)
)
featureData(gene_eset) <- new("AnnotatedDataFrame", data = fdata)

# Save the gene-level Eset
saveRDS(gene_eset, output_file)

# sanity check: sum of gene-level expression should be equal to sum of transcript-level expression for each sample
sum_tx <- colSums(exprs(eset))
sum_gene <- colSums(exprs(gene_eset))
if (!isTRUE(all.equal(sum_tx, sum_gene))) {
  warning("Sum of gene-level expression does not match sum of transcript-level expression")
} else {
  print("Sanity check passed: Sum of gene-level expression matches sum of transcript-level expression")
}