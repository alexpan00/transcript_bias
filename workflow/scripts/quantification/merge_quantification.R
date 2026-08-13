
# Packages
library(tidyverse)
library(stringr)

# Input
args <- commandArgs(trailingOnly = TRUE)
quantification_fofn <- args[1]
tama_merge_file <- args[2]
out_file <- args[3]

quantification_fofn <- readLines(quantification_fofn)
tama_merge_table <- read.table(tama_merge_file, header = FALSE, sep = "\t")

# Format TAMA merge dataframe

# Get sample IDs and sort them by length (descending) to match the longest one first
all_sample_ids <- sapply(strsplit(quantification_fofn, "\t"), `[`, 1)
all_sample_ids <- all_sample_ids[order(nchar(all_sample_ids), decreasing = TRUE)]

# Parse TAMA merge table column 4 into new_id, sample_id, and old_id
tama_id_map <- bind_rows(lapply(as.character(tama_merge_table[, 4]), function(x) {
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  new_id <- parts[1]
  rest <- parts[2]
  
  matched_sample <- NA
  old_id <- NA
  
  for (sid in all_sample_ids) {
    prefix <- paste0(sid, "_")
    if (startsWith(rest, prefix)) {
      matched_sample <- sid
      old_id <- substr(rest, nchar(prefix) + 1, nchar(rest))
      break
    }
  }
  return(data.frame(new_id = new_id, sample_id = matched_sample, old_id = old_id, stringsAsFactors = FALSE))
}))

print("head of tama_id_map")
print(head(tama_id_map))

# Change old transcript id to merged transcript id
quant_mat_list <- list()
for (i in 1:length(quantification_fofn)) {
  quant_fofn_line <- strsplit(quantification_fofn[i], "\t")[[1]]
  sample_id <- quant_fofn_line[1]
  file <- quant_fofn_line[2]
  if (!file.exists(file) || file.info(file)$size == 0) {
    warning(paste("Quantification file missing or empty:", file, "- skipping."))
    next
  }
  quant_file <- read.table(file, header = TRUE, sep = "\t", comment.char = "")
  if (nrow(quant_file) == 0) {
    warning(paste("Quantification file has 0 rows:", file, "- skipping."))
    next
  }

  colnames(quant_file)[1] <- "pbid"
  print("colnames of quant file")
  print(colnames(quant_file))

  tama_name_by_sample <- tama_id_map[tama_id_map$sample_id == sample_id, ]
  match_ids <- match(quant_file$pbid, tama_name_by_sample$old_id)
  quant_file$transcript_id <- tama_name_by_sample$new_id[match_ids]
  if ("pbid" %in% colnames(quant_file)) {
    quant_file$pbid <- NULL
  }

  print("head of quant_file after transcript_id reassignment")
  print(head(quant_file))

  quant_file <- quant_file %>%
    pivot_longer(-transcript_id, names_to = "sample", values_to = "count") %>%
    group_by(transcript_id, sample) %>%
    summarise(tot_count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = "sample", values_from = "tot_count")

  quant_mat_list[[length(quant_mat_list) + 1]] <- as.data.frame(quant_file)
}

# Merge quantification matrices
quant_mat_list <- quant_mat_list[!sapply(quant_mat_list, is.null)]
if (length(quant_mat_list) == 0) {
  quant_merged_matrix <- data.frame(transcript_id = character(0))
} else {
  quant_merged_matrix <- quant_mat_list %>%
    reduce(full_join, by = "transcript_id")
  quant_merged_matrix[is.na(quant_merged_matrix)] <- 0
  quant_merged_matrix <- quant_merged_matrix %>% relocate(transcript_id)
}

# rename to fit sqanti FL format (will always be more than 1 sample)
names(quant_merged_matrix)[names(quant_merged_matrix) == 'transcript_id'] <- 'superPBID'

# Write output
write.table(quant_merged_matrix, file = out_file,
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
