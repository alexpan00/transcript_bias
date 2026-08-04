library(tidyverse)

# Arguments:
# 1. Path to input file (2 columns: Path to summary CSV, Dataset ID)
# 2. Output PNG path for Sum Counts
# 3. Output PNG path for Mean Expression
args <- commandArgs(trailingOnly = TRUE)
input_file_path <- args[1]
output_png_sum <- args[2]
output_png_mean <- args[3]

# Read input file
# Assuming no header, and space/tab separated
# Columns: 1 = File Path, 2 = Dataset ID
input_data <- read.table(input_file_path, header = FALSE, stringsAsFactors = FALSE)
colnames(input_data) <- c("Path", "DatasetID")

all_summaries <- list()

for (i in 1:nrow(input_data)) {
  file_path <- input_data$Path[i]
  dataset_id <- input_data$DatasetID[i]
  
  if (file.exists(file_path)) {
    # Read the summary CSV
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Add DatasetID column
    df$Dataset <- dataset_id
    
    all_summaries[[i]] <- df
  } else {
    warning(paste("File not found:", file_path))
  }
}

# Combine all dataframes
final_combined_df <- do.call(rbind, all_summaries)

if (is.null(final_combined_df) || nrow(final_combined_df) == 0) {
  stop("No valid data loaded. Check input file paths.")
}

# Ensure columns are factors for proper ordering/grouping
final_combined_df$Dataset <- as.factor(final_combined_df$Dataset)
final_combined_df$Tool <- as.factor(final_combined_df$Tool)
final_combined_df$Normalization <- as.factor(final_combined_df$Normalization)
final_combined_df$Length_bins <- factor(final_combined_df$Length_bins, levels = unique(final_combined_df$Length_bins ))



# Plot 1: Sum Counts
p1 <- final_combined_df %>% 
  ggplot(aes(x = Length_bins, y = sum_counts, color = Normalization, group = Normalization)) +
  geom_line(linewidth = 1) +
  facet_grid(Condition + Dataset ~ Tool) +
  theme_light() +
  labs(title = "Total Counts per Length Bin", x = "Length Bin", y = "Total Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        strip.text = element_text(size = 12))

ggsave(output_png_sum, plot = p1, height = max(8, 4 * length(levels(final_combined_df$Condition))), width = 16)


# Plot 2: Mean Expression
p2 <- final_combined_df %>% 
  ggplot(aes(x = Length_bins, y = mean_exp, color = Normalization, group = Normalization)) +
  geom_line(linewidth = 1) +
  facet_grid(Condition + Dataset ~ Tool, scales = "free_y") +
  theme_light() +
  labs(title = "Mean Counts per Length Bin", x = "Length Bin", y = "Mean Counts") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        strip.text = element_text(size = 12))

ggsave(output_png_mean, plot = p2, height = max(8, 4 * length(levels(final_combined_df$Condition))), width = 16)