library(tidyverse)
library(NOISeq)

# Arguments are:
# 1. Comma separated list of SQANTI CSV summary objects
# 2. Output PNG path for stacked plot across tools
# 3. Output summary table path
args <- commandArgs(trailingOnly = TRUE)
sqanti_csvs <- args[1]
output_png_combined <- args[2]
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

sqanti_csvs <- unlist(strsplit(sqanti_csvs, ","))
l_res <- list()

# For each CSV file, get tool info from path and add it
for (csv_file in sqanti_csvs){
  # get tool from the path (e.g., results/bambu/NOIseq/raw/...)
  meta <- parse_path_metadata(csv_file)
  tool <- meta$tool
  
  # Read SQANTI data (Sample, Biotype, Count, Percentage)
  sq_data <- read.csv(csv_file)

  # add tool info
  sq_data$Tool <- tool
  l_res[[csv_file]] <- sq_data
}

# Combine results
final_df <- do.call(rbind, l_res)

# Source the palette if possible, otherwise define it (simplified for the summary)
# Using similar palette as sqanti_plots.R
cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")

# Create a combined plot comparing tools
# Instead of averaging, we plot each Sample and facet by Tool (in different rows)
# Switched from Percentage to Count as requested
p <- ggplot(final_df, aes(x = Sample, y = Count, fill = Biotype)) +
  geom_bar(stat = "identity", color = "black", size = 0.3, width = 0.7) +
  facet_grid(Tool ~ ., scales = "free_y") +
  xlab("Sample") +
  ylab("Transcripts, #") +
  theme_light() +
  scale_fill_manual(values = cat.palette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_png_combined, p, height = 8, width = 10)

# Create boxplot of total transcripts by condition
# Sum counts across all biotypes for each Tool, Condition, and Sample
final_df_aggregated <- final_df %>%
  group_by(Tool, Condition, Sample) %>%
  summarise(N_Transcripts = sum(Count), .groups = 'drop')

# Create boxplot with jitter for N_Transcripts, conditions in y-axis
final_df_aggregated %>%
  ggplot(aes(x = N_Transcripts, y = Condition)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = Sample), height = 0.2, width = 0, size = 2) +
    facet_grid(Tool ~ .) +
    theme_light() +
    labs(x = "Number of Transcripts", y = "Condition", title = "Distribution of Transcripts by Condition and Tool")

# Get the output path for the boxplot PNG
output_dir <- dirname(output_png_combined)
output_boxplot <- file.path(output_dir, "sqanti_condition_boxplot_ntranscripts.png")
ggsave(output_boxplot, height = 6, width = 8)

# Write full combined table
write.csv(final_df, output_summary_table, row.names = FALSE)
