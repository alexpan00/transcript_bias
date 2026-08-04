library(tidyverse)

# Arguments:
# 1. Path to input file (2 columns: Path to summary CSV, Dataset ID)
# 2. Output PNG path
args <- commandArgs(trailingOnly = TRUE)
input_file_path <- args[1]
output_png_path <- args[2]

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
mix_col <- which(colnames(final_combined_df) == "Mixture")
if (length(mix_col) != 0){
  colnames(final_combined_df)[mix_col] <- "Condition"
}
# Ensure Tool and Dataset are factors for proper ordering if needed
final_combined_df$Dataset <- as.factor(final_combined_df$Dataset)
final_combined_df$Tool <- as.factor(final_combined_df$Tool)

# Create the plot
# Using the style from summary_sr.R as a base
# Facet by Dataset and Condition to show comparisons
p <- final_combined_df %>% 
  filter(!is.na(Correlation)) %>% 
  ggplot(aes(x = Tool, y = Correlation)) +
  geom_boxplot() +
  # Use shape for Normalization, and maybe color for Sample or just keeping it simple
  geom_jitter(aes(color = Normalization), width = 0.2, size = 2, height = 0) +
  facet_grid(Dataset~Condition) +
  theme_light() +
  labs(y = "Spearman Correlation", x = "Tool") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        strip.text = element_text(size = 14))

# Save the plot
ggsave(output_png_path, plot = p, height = 8, width = 12) # Adjusted width/height for potential multiple datasets


if ("RMSE" %in% colnames(final_combined_df)){
  output_rmse_png_path <- if (length(args) >= 3) args[3] else sub("\\.png$", "_rmse.png", output_png_path)
  p_rmse <- final_combined_df %>% 
    ggplot(aes(x = Tool, y = RMSE)) +
    geom_boxplot() +
    # Use shape for Normalization, and maybe color for Sample or just keeping it simple
    geom_jitter(aes(color = Normalization), width = 0.2, size = 2, height = 0) +
    facet_grid(Dataset~Condition) +
    theme_light() +
    labs(y = "RMSE", x = "Tool") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
          strip.text = element_text(size = 14))
  ggsave(output_rmse_png_path, plot = p_rmse, height = 8, width = 12)
}