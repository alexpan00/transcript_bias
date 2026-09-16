library(tidyverse)
library(NOISeq)
library(mgcv)

get_script_dir <- function() {
  if (exists("snakemake") && !is.null(snakemake@script)) {
    return(dirname(snakemake@script))
  }
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(sub("^--file=", "", file_arg[1])))
  }
  possible_dirs <- c("workflow/scripts", "scripts", ".")
  for (d in possible_dirs) {
    if (dir.exists(d)) return(d)
  }
  return(".")
}

source_script <- function(script_name) {
  categories <- c(".", "plotting", "analysis", "normalization", "preprocessing", "quantification", "summary", "utils")
  s_dir <- get_script_dir()
  for (cat in categories) {
    cands <- c(
      file.path(s_dir, cat, script_name),
      file.path(s_dir, script_name),
      file.path("workflow/scripts", cat, script_name),
      file.path("scripts", cat, script_name),
      file.path(cat, script_name)
    )
    for (cand in cands) {
      if (file.exists(cand)) {
        source(cand)
        return(invisible(TRUE))
      }
    }
  }
  source(script_name)
}

source_script("scale_utils.R")

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
output_png_sum <- args[2]
output_png_mean <- args[3]
output_csv <- args[4]

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

# For each combination of tool + normalization get counts per length bin
for (long_obj in long_objs){
  # get tool and normalization method from the path
  meta <- parse_path_metadata(long_obj)
  tool <- meta$tool
  normalization_method <- meta$norm
  # Read long reads data and normalize using CPM
  mydata <- readRDS(long_obj)
  mydata_cpm <- exprs(mydata)

  # Put every method on a common CPM scale, so that what is compared between
  # methods is the part of the normalization that is not sequencing depth.
  # Values already on a log scale (ratio_correction) are left untouched.
  if (!is_log_scale(mydata_cpm)) {
    mydata_cpm <- cpm(mydata_cpm)
  }
  # select the main factor
  main_factor <- if (ncol(pData(mydata)) >= 2) pData(mydata)[,2] else pData(mydata)[,1]
  main_factor <- as.factor(main_factor)
  conds <- levels(main_factor)

  # get the mean counts per condition
  datos = data.frame(sapply(conds, 
                function (k) {
                  rowMeans(as.matrix(mydata_cpm[, main_factor == k]))
                }), check.names = FALSE)
  datos$Length <- fData(mydata)$Length
  # bin the transcript length every 500 and make everything above 10000 to be in the same bin
  datos$Length_bins <- cut(
                          datos$Length,
                          breaks = c(seq(0, 5000, 500), Inf),
                          labels = c(paste0(seq(0, 4.5, 0.5),"k-", seq(0.5, 5, 0.5),"k"), "5k+"),
                          include.lowest = TRUE)
  datos_melt <- pivot_longer(datos, cols = conds, names_to = "Condition", values_to = "Counts")
  datos_summary <- datos_melt %>% 
    group_by(Condition, Length_bins) %>% 
    summarise(sum_counts = sum(Counts),
              mean_exp = ifelse(all(Counts == 0), 0, mean(Counts[Counts != 0], trim = 0.025, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate("Tool" = tool, "Normalization" = normalization_method)
  # Inlcude global summary for all conditions
  datos_summary_global <- datos_melt %>% 
    group_by(Length_bins) %>% 
    summarise(sum_counts = sum(Counts),
              mean_exp = ifelse(all(Counts == 0), 0, mean(Counts[Counts != 0], trim = 0.025, na.rm = TRUE))) %>% 
    ungroup() %>% 
    mutate("Tool" = tool, "Normalization" = normalization_method, "Condition" = "All")
  datos_summary <- rbind(datos_summary, datos_summary_global)
  l_res[[long_obj]] <- datos_summary
}

# Combine results
final_df <- do.call(rbind, l_res)

final_df %>% 
  ggplot(aes(x = Length_bins, y = sum_counts, color = Normalization, group = Normalization)) +
    geom_line(position = position_dodge(width = 0.2), linewidth = 1) +
    facet_grid(Condition ~ Tool) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_png_sum, height = 8, width = 24)

final_df %>% 
  ggplot(aes(x = Length_bins, y = mean_exp, color = Normalization, group = Normalization)) +
    geom_line(position = position_dodge(width = 0.2), linewidth = 1) +
    facet_grid(Condition ~ Tool) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(output_png_mean, height = 8, width = 24)

# Sanity check: print total counts per tool + normalization
print(final_df %>% 
  group_by(Tool, Normalization) %>% 
  summarise(total_counts = sum(sum_counts)))

# Write summary to CSV
write.csv(final_df, output_csv, row.names = FALSE)