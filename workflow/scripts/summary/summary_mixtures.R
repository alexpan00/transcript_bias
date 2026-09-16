library(tidyverse)
library(NOISeq)


cpm <- function(df){
  df <- as.matrix(df)
  df_cpm <- t(10^6*t(df)/colSums(df))
  return(df_cpm)
}

# Arguemnts are:
# 1. Comma separated list of long read RDS objects
# 2. Output PNG path
args <- commandArgs(trailingOnly = TRUE)
long_objs <- args[1]
output_png_correlation <- args[2]
output_summary_table <- args[3]

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

# For each combination of tool + normalization get mixture summary
# and add tool and normalization info
for (long_obj in long_objs){
  # get tool and normalization method from the path
  meta <- parse_path_metadata(long_obj)
  tool <- meta$tool
  normalization_method <- meta$norm
  # Read mixture data (columns are measured, expected and condtion)
  mixture_summary <- readRDS(long_obj)

  # add tool and normalization info
  mixture_summary$Tool <- tool
  mixture_summary$Normalization <- normalization_method
  l_res[[long_obj]] <- mixture_summary
}

# Combine results
final_df <- do.call(rbind, l_res)

# Compute correlation for each combination of tool and normalization method
final_df_summary <- final_df %>%
  group_by(Tool, Normalization, Condition, Sample) %>%
  summarise(
    Correlation = cor(Measured, Expected, method = "spearman")
  ) %>%
  ungroup()

# Create boxplot of correlation wrapped by Mixture
final_df_summary %>% 
  ggplot(aes(x = Tool, y = Correlation)) +
    geom_boxplot() +
    geom_jitter(aes(shape = Normalization), width = 0.2, size = 2, height = 0) +
    facet_grid(Condition ~ .) +
    theme_light() +
    labs(y = "Spearman Correlation", x = "Tool")

ggsave(output_png_correlation, height = 6, width = 8)

# write summary table to the same folder as the output PNG
write.csv(final_df_summary, output_summary_table, row.names = FALSE)