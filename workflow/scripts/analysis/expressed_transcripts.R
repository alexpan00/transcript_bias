# Arguments are:
# 1. Counts file path (TSV)
# 2. Output sqanti categories of expressed transcripts (TSV)
args <- commandArgs(trailingOnly = TRUE)
counts <- args[1]
output_categories <- args[2]

# read the counts
counts <- read.table(counts, header = TRUE, row.names = 1, sep = "\t")

# identify expressed transcripts (non-zero counts)
expressed_transcripts <- rownames(counts)[rowSums(counts) > 0]
output_df <- data.frame(isoform = expressed_transcripts, structural_category = "full-splice_match")

# write the output
write.table(output_df, output_categories, sep = "\t", row.names = FALSE, quote = FALSE)