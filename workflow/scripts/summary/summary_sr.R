library(tidyverse)
library(NOISeq)


cpm <- function(df){
  df <- as.matrix(df)
  df_cpm <- t(10^6*t(df)/colSums(df))
  return(df_cpm)
}

# Arguemnts are:
# 1. Comma separated list of long read RDS objects
# 2. Output prefix for PNGs (e.g., results/summary_sr)
# 3. Output summary table path
args <- commandArgs(trailingOnly = TRUE)
long_objs <- args[1]
output_prefix <- args[2]
output_summary_table <- args[3]

parse_path_metadata <- function(path_str) {
  tools_list <- c("bambu", "flair", "isoseq", "isoquant", "kallisto", "oarfish", "tama", "sqanti")
  norms_list <- c("raw", "cpm", "tmm", "tpm", "ratio_correction", "ratio_counts", "cqn", "eda", "read_density", "optimal_epsilon")
  parts <- unlist(strsplit(path_str, "/"))
  parts_lower <- tolower(parts)
  found_tool <- NA_character_
  for (t in tools_list) { if (t %in% parts_lower) { found_tool <- t; break } }
  found_norm <- NA_character_
  for (n in norms_list) { if (n %in% parts_lower) { found_norm <- n; break } }
  if (is.na(found_tool)) found_tool <- if (length(parts) >= 4) parts[length(parts) - 3] else if (length(parts) >= 2) parts[length(parts) - 1] else "unknown"
  if (is.na(found_norm)) found_norm <- if (length(parts) >= 2) parts[length(parts) - 1] else "raw"
  return(list(tool = found_tool, norm = found_norm))
}

long_objs <- unlist(strsplit(long_objs, ","))
l_res <- list()

# For each combination of tool + normalization get short-reads summary and add tool and normalization info
for (long_obj in long_objs){
  # get tool and normalization method from the path
  meta <- parse_path_metadata(long_obj)
  tool <- meta$tool
  normalization_method <- meta$norm
  # Read short-reads data (columsn are measured, expected and condtion)
  sr_summary <- readRDS(long_obj)

  # add tool and normalization info
  sr_summary$Tool <- tool
  sr_summary$Normalization <- normalization_method
  l_res[[long_obj]] <- sr_summary
}

# Combine results
final_df <- do.call(rbind, l_res)

# Compute correlation for each combination of tool and normalization method
final_df_summary <- final_df %>%
  group_by(Tool, Normalization, Condition, Sample) %>%
  summarise(
    Correlation = cor(Measured, Expected, method = "spearman"),
    # include number of transcripts in the correlation
    N_Transcripts = n()
  ) %>%
  mutate(LengthQuantile = "All") %>%
  ungroup()

# Create boxplot of correlation wrapped by Condition
final_df_summary %>% 
  ggplot(aes(x = Tool, y = Correlation)) +
    geom_boxplot() +
    geom_jitter(aes(shape = Normalization, color = Sample), width = 0.2, size = 2, height = 0) +
    facet_grid(Condition ~ .) +
    theme_light() +
    labs(y = "Spearman Correlation", x = "Tool")

ggsave(paste0(output_prefix, "_normalization_sr_correlation.png"), height = 6, width = 8)

# we only need one normalization method to show N_Transcripts as they should be the same
final_df_transcripts <- final_df_summary %>%
  group_by(Tool, Condition, Sample) %>%
  filter(Normalization == first(Normalization)) %>%
  ungroup()

# Create barplot of N_Transcripts by sample and condition
final_df_transcripts %>%
  ggplot(aes(x = Sample, y = N_Transcripts)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_grid(Tool ~ .) +
    theme_light() +
    labs(y = "Number of Transcripts", x = "Sample", title = "Transcripts per Sample and Condition") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_prefix, "_barplot_ntranscripts.png"), height = 6, width = 10)

# Create boxplot with jitter for N_Transcripts, conditions in y-axis
final_df_transcripts %>%
  ggplot(aes(x = N_Transcripts, y = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Sample), height = 0.2, width = 0, size = 2) +
    facet_grid(Tool ~ .) +
    theme_light() +
    labs(x = "Number of Transcripts", y = "Condition", title = "Distribution of Transcripts by Condition and Tool")

ggsave(paste0(output_prefix, "_condition_boxplot_ntranscripts.png"), height = 6, width = 8)

if ("LengthQuantile" %in% colnames(final_df)) {
  final_df_summary_length_quantile <- final_df %>%
    filter(!is.na(LengthQuantile)) %>%
    group_by(Tool, Normalization, Condition, Sample, LengthQuantile) %>%
    summarise(
      Correlation = cor(Measured, Expected, method = "spearman"),
      N_Transcripts = n()
    ) %>%
    ungroup()
  final_df_summary <- bind_rows(final_df_summary, final_df_summary_length_quantile)

}

# write summary table to the same folder as the output PNG
write.csv(final_df_summary, output_summary_table, row.names = FALSE)

# sanity check: all the normalization methods should have the same number 
# of transcripts in the correlation, if not there is an issue with the data
final_df_summary %>%
  filter(LengthQuantile == "All") %>%
  group_by(Tool, Normalization) %>%
  summarise(mean_n_transcripts = mean(N_Transcripts),
            deviation_n_transcripts = sd(N_Transcripts)) %>%
  print()

