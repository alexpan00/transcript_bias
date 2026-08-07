library(tidyverse)
library(patchwork)

# Fisher z-transformation helpers
fisher_z <- function(r) atanh(r)
fisher_z_inv <- function(z) tanh(z)

# Compute mean correlation using Fisher z-transformation
mean_cor_fisher <- function(r_values) {
  r_values <- r_values[!is.na(r_values) & abs(r_values) < 1]
  if (length(r_values) == 0) return(NA_real_)
  fisher_z_inv(mean(fisher_z(r_values)))
}

# function to rank and convert rank to factor
rank_column <- function(x, n=3){
  ranking <- rank(x)
  ranking_factor <- factor(ifelse(ranking <= n, ranking, "Rest"))
  return(ranking_factor)
}

# Arguments:
# 1. Comma-separated list of sr_summary.rds files (transcript level)
# 2. Comma-separated list of gene_level_sr_summary.rds files
# 3. Output PNG path
# 4. Output CSV path
# 5. Comma-separated list of mixture_summary.rds files ("none" to skip)
# 6. Comma-separated list of sirv_summary.rds files ("none" to skip)
args <- commandArgs(trailingOnly = TRUE)
sr_transcript_objs  <- args[1]
sr_gene_objs        <- args[2]
output_png          <- args[3]
output_csv          <- args[4]
mixture_objs        <- if (length(args) >= 5 && args[5] != "none") args[5] else NULL
sirv_objs           <- if (length(args) >= 6 && args[6] != "none") args[6] else NULL
ercc_objs           <- if (length(args) >= 7 && args[7] != "none") args[7] else NULL
tusco_objs          <- if (length(args) >= 8 && args[8] != "none") args[8] else NULL

# Load summary CSV files, tagging each row with DataType
load_summary_files <- function(file_path, data_type) {
  df <- read.csv(file_path)
  # Ensure the data has the same columns
  # Some files use "Correlation" others use "MeanCorrelation"
  if ("Correlation" %in% colnames(df)) {
    df <- rename(df, MeanCorrelation = Correlation)
  }
  # Filter for 'All' length quantile if present
  if ("LengthQuantile" %in% colnames(df)) {
    df <- df[df$LengthQuantile == "All", ]
  }
  # If this is SIRV or ERCC data and a precomputed Similarity column exists,
  # use that as the metric to plot (store it in MeanCorrelation for downstream code)
  if (toupper(data_type) %in% c("SIRV", "ERCC") && "Similarity" %in% colnames(df)) {
    df$MeanCorrelation <- df$Similarity
  }
  
  # Step: average correlations with Fisher z-transformation if Sample is present
  # The input is now expected to be a pre-summarized CSV or a per-sample CSV.
  # If 'Sample' is present, we need to average it per Tool/Normalization/Condition.
  if ("Sample" %in% colnames(df)) {
    if (toupper(data_type) %in% c("SIRV", "ERCC")) {
      # For similarity-based metrics use arithmetic mean across samples
      df <- df %>%
        group_by(Tool, Normalization) %>%
        summarise(
          MeanCorrelation = mean(MeanCorrelation, na.rm = TRUE),
          N_Samples       = n(),
          .groups         = "drop"
        )
    } else {
      # For correlation-like metrics average using Fisher z-transform
      df <- df %>%
        group_by(Tool, Normalization) %>%
        summarise(
          MeanCorrelation = mean_cor_fisher(MeanCorrelation),
          N_Samples       = n(),
          .groups         = "drop"
        )
    }
  }

  df$DataType <- data_type
  return(df)
}

# Required data sources
df_sr_tx   <- load_summary_files(sr_transcript_objs, "SR (transcript)")
df_sr_tx$ranking <- rank_column(1-df_sr_tx$MeanCorrelation)
df_sr_gene <- load_summary_files(sr_gene_objs,       "SR (gene)")
df_sr_gene$ranking <- rank_column(1-df_sr_gene$MeanCorrelation)

all_dfs <- list(df_sr_tx, df_sr_gene)

if (!is.null(mixture_objs)) {
  df_mix <- load_summary_files(mixture_objs, "Mixture")
  df_mix$ranking <- rank_column(1-df_mix$MeanCorrelation)
  all_dfs <- c(all_dfs, list(df_mix))
}

if (!is.null(sirv_objs)) {
  df_sirv <- load_summary_files(sirv_objs, "SIRV")
  df_sirv$ranking <- rank_column(1-df_sirv$MeanCorrelation)
  all_dfs <- c(all_dfs, list(df_sirv))
}

if (!is.null(ercc_objs)) {
  df_ercc <- load_summary_files(ercc_objs, "ERCC")
  df_ercc$ranking <- rank_column(1-df_ercc$MeanCorrelation)
  all_dfs <- c(all_dfs, list(df_ercc))
}

if (!is.null(tusco_objs)) {
  df_tusco <- load_summary_files(tusco_objs, "TUSCO")
  df_tusco$ranking <- rank_column(1-df_tusco$MeanCorrelation)
  all_dfs <- c(all_dfs, list(df_tusco))
}

# Keep only the columns needed for visualization
keep_cols <- c("MeanCorrelation", "Tool", "Normalization", "Condition", "DataType", "ranking")
all_dfs <- lapply(all_dfs, function(df) {
  df[, intersect(colnames(df), keep_cols), drop = FALSE]
})

summary_df <- do.call(rbind, all_dfs)

# Ordered factor so DataType rows appear in a logical sequence
dt_levels <- intersect(
  c("SR (transcript)", "SR (gene)", "Mixture", "SIRV", "ERCC", "TUSCO"),
  unique(summary_df$DataType)
)
summary_df$DataType <- factor(summary_df$DataType, levels = dt_levels)

# Tool abbreviations helper
tool_abbrev <- c(
  "isoseq"   = "is",
  "kallisto" = "ka",
  "oarfish"  = "of",
  "isoquant" = "iq",
  "bambu"    = "ba",
  "flair"    = "fl"
)

summary_df <- summary_df %>%
  mutate(ToolAbbrev = ifelse(tolower(Tool) %in% names(tool_abbrev), 
                             tool_abbrev[tolower(Tool)], 
                             Tool))

p1 <- summary_df %>%
  ggplot(aes(x = Normalization, y = DataType, fill = MeanCorrelation, shape = ranking)) +
  geom_point(size = 3) +
  scale_fill_viridis_c(
    option   = "D",
    limits   = c(0, 1),
    na.value = "grey80",
    name     = "Mean Spearman\nCorrelation"
  ) +
  facet_grid(DataType ~ ToolAbbrev, scales = "free_x") +
  theme_light() +
  theme(
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust=1),
    axis.ticks.x     = element_blank(),
    strip.text       = element_text(size = 10),
    legend.position  = "right"
  ) +
  labs(
    x     = "",
    y     = ""
  ) + scale_shape_manual(values = c("1" = 24, "2" = 23, "3" = 22, "Rest" = 21))
n_data_types  <- length(unique(summary_df$DataType))
n_norms       <- length(unique(summary_df$Normalization))
n_tools       <- length(unique(summary_df$Tool))
ggsave(
  paste0(output_png, "_alternative.png"),
  plot   = p1,
  height = min(12, 1.5 * n_data_types + 1.5),
  width  = min(20, 1.2 * n_norms * n_tools + 2),
  dpi    = 650
)

annot <-  summary_df %>% 
  mutate(tool_norm = paste(Tool, Normalization, sep = "_")) %>% 
  select(tool_norm, Normalization, Tool) %>% 
  unique() %>% 
  column_to_rownames("tool_norm")

cormat <- as.matrix(summary_df %>% 
  mutate(tool_norm = paste(Tool, Normalization, sep = "_")) %>% 
  select(MeanCorrelation, DataType, tool_norm) %>% 
  pivot_wider(names_from = tool_norm, values_from = MeanCorrelation) %>% 
  column_to_rownames("DataType"))

# Generate dynamic color palettes for Tool and Normalization (non-overlapping hue ranges)
tool_unique <- unique(annot$Tool)
norm_unique <- unique(annot$Normalization)

# Tool: red to yellow range (0-120 hue)
tool_colors <- setNames(
  scales::hue_pal(h = c(0, 120), l = 60)(length(tool_unique)),
  tool_unique
)

# Normalization: Use a qualitative palette with high contrast (e.g., ColorBrewer Set1 or Dark2)
norm_colors <- setNames(
  scales::brewer_pal(palette = "Set1")(length(norm_unique)),
  norm_unique
)

annotation_colors <- list(
  Tool = tool_colors,
  Normalization = norm_colors
)

# ============================================================================
# ALTERNATIVE: ggplot2 version with geom_tile and shapes (ranking overlay)
# ============================================================================
# Prepare data in long format for ggplot
heatmap_data <- summary_df %>%
  mutate(tool_norm = paste(Tool, Normalization, sep = "_")) %>%
  select(MeanCorrelation, DataType, tool_norm, ranking, Tool, Normalization) %>%
  unique()

# Shape mapping (based on ranking)
shape_map <- c("1" = 24, "2" = 23, "3" = 22, "Rest" = NA)

# Build annotation data matching column order from heatmap
annot_long <- heatmap_data %>%
  select(tool_norm, Tool, Normalization) %>%
  distinct() %>%
  arrange(match(tool_norm, colnames(cormat)))

# Convert to factor to preserve order
annot_long$tool_norm <- factor(annot_long$tool_norm, levels = annot_long$tool_norm)

# Main heatmap plot
p_tile <- heatmap_data %>%
  mutate(tool_norm = factor(tool_norm, levels = levels(annot_long$tool_norm))) %>%
  ggplot(aes(x = tool_norm, y = DataType, fill = MeanCorrelation, shape = ranking)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_point(size = 4) +
  scale_fill_viridis_c(
    option = "D",
    limits = c(0, 1),
    na.value = "grey80",
    name = "Mean Spearman\nCorrelation"
  ) +
  scale_shape_manual(
    values = shape_map,
    name = "Ranking"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(0, 5.5, 0, 5.5)
  )
p_tile_annotated <- heatmap_data %>%
  mutate(tool_norm = factor(tool_norm, levels = levels(annot_long$tool_norm))) %>%
  ggplot(aes(x = Normalization, y = DataType, fill = MeanCorrelation, shape = ranking)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_point(size = 4) +
  scale_fill_viridis_c(
    option = "D",
    limits = c(0, 1),
    na.value = "grey80",
    name = "Mean Spearman\nCorrelation"
  ) +
  scale_shape_manual(
    values = shape_map,
    name = "Ranking"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_light() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
    axis.ticks.x = element_blank(),
    axis.title = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(0, 5.5, 0, 5.5)
  ) + facet_grid(DataType ~ Tool)


# Uncomment below to save the ggplot version with annotations instead:
ggsave(
   output_png,
   plot = p_tile_annotated,
   height =  min(10, 0.2 * n_data_types + 1.5),
   width = min(20, 0.3 * n_norms * n_tools + 2),
   dpi = 650
)

write.csv(heatmap_data, output_csv, row.names = FALSE)
