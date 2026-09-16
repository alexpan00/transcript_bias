library(tidyverse)
library(NOISeq)

# Arguemnts are:
# 1. Comma separated list of long read RDS objects
# 2. Output PNG path
args <- commandArgs(trailingOnly = TRUE)
long_objs <- args[1]
output_summary_table <- args[2]

parse_path_metadata <- function(path_str) {
  tools_list <- c("bambu", "flair", "isoseq", "isoquant", "kallisto", "oarfish", "tama", "sqanti")
  norms_list <- c("raw", "cpm", "tpm", "ratio_correction", "ratio_counts", "cqn", "eda")
  parts <- unlist(strsplit(path_str, "/"))
  parts_lower <- tolower(parts)
  found_tool <- NA_character_
  for (t in tools_list) {
    idx <- which(parts_lower == t)
    if (length(idx) > 0) { found_tool <- parts[idx[1]]; break }
  }
  found_norm <- NA_character_
  for (n in norms_list) {
    idx <- which(parts_lower == n)
    if (length(idx) > 0) { found_norm <- parts[idx[1]]; break }
  }
  if (is.na(found_tool)) found_tool <- if (length(parts) >= 4) parts[length(parts) - 3] else if (length(parts) >= 2) parts[length(parts) - 1] else "unknown"
  if (is.na(found_norm)) found_norm <- if (length(parts) >= 2) parts[length(parts) - 1] else "raw"
  return(list(tool = found_tool, norm = found_norm))
}

long_objs <- unlist(strsplit(long_objs, ","))
l_res <- list()

# For each combination of tool + normalization get sirv summary and add tool and normalization info
for (long_obj in long_objs){
  # get tool and normalization method from the path
  meta <- parse_path_metadata(long_obj)
  tool <- meta$tool
  # Read sensitivity data
  sirv_summary <- read.csv(long_obj, header = TRUE, row.names = NULL)

  # add tool and normalization info
  sirv_summary$Tool <- tool
  l_res[[long_obj]] <- sirv_summary
}

# Combine results
final_df <- do.call(rbind, l_res)

# save summary table 
write.csv(final_df, output_summary_table, row.names = FALSE)
