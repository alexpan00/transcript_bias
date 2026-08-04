library(NOISeq)
library(ggplot2)
library(dplyr)
library(tidyr)

remove_zeros <- function(datos, condition = NULL){
  ceros = which(rowSums(datos) == 0)
  if (length(ceros) > 0) {
    if (!is.null(condition)){
      print(paste("Warning:", length(ceros), 
                  "features with 0 counts in", condition, "are to be removed for this analysis."))
    } else{
      print(paste("Warning:", length(ceros), 
                "features with 0 counts in all samples are to be removed for this analysis."))
    }
    datos = datos[-ceros,]
  }
  return(datos)
}


cpm <- function(df){
  df <- as.matrix(df)
  df_cpm <- t(10^6*t(df)/colSums(df))
  return(df_cpm)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript mixture_correlation.R <observed_norm.rds> <synthetic_norm.rds> <mixture_definition.csv> <output_prefix>")
}

obs_norm_path <- args[1]
syn_norm_path <- args[2]
mix_def_path <- args[3]
output_prefix <- args[4]

# Read data
obs_obj <- readRDS(obs_norm_path)
syn_obj <- readRDS(syn_norm_path)
mix_def <- read.csv(mix_def_path)

required_cols <- c("Mixture", "Component", "Proportion")
if (!all(required_cols %in% colnames(mix_def))) {
  stop("Mixture definition file must contain columns: ", paste(required_cols, collapse = ", "))
}

# Get expressions
# If normalization was performed, it should be in exprs
get_exprs <- function(obj) {
  if (!is.null(Biobase::assayData(obj)$exprs)) return(Biobase::assayData(obj)$exprs)
  return(Biobase::assayData(obj)$counts)
}

data_obs <- get_exprs(obs_obj)
data_syn <- get_exprs(syn_obj)

# Get metadata for labels
metadata_obs <- pData(obs_obj)
cond_col <- colnames(metadata_obs)[2]
sample_col <- colnames(metadata_obs)[1]

mixture_conditions <- unique(mix_def$Mixture)
mixture_samples <- metadata_obs[[sample_col]][metadata_obs[[cond_col]] %in% mixture_conditions]

# Results storage
correlation_results <- data.frame()
l_res <- list()

# Match samples
common_samples <- intersect(colnames(data_obs), colnames(data_syn))
common_samples <- intersect(common_samples, mixture_samples)

if (length(common_samples) == 0) {
  stop("No mixture samples found in common between observed and synthetic objects.")
}

for (sname in common_samples) {
    df_comp <- data.frame(Measured = data_obs[, sname], Expected = data_syn[, sname])

    df_comp <- remove_zeros(df_comp, sname)
    
    # Log transformation (only if the method is not ratio correction, which can have negative values)
    if (all(df_comp >= 0)){
      df_comp_log <- log2(cpm(df_comp) + 1)
    } else {
      df_comp_log <- df_comp
    }
    
    mix_correlation <- cor(df_comp_log[, "Measured"], df_comp_log[, "Expected"], method = "spearman")
    
    cond <- metadata_obs[[cond_col]][match(sname, metadata_obs[[sample_col]])]
    
    correlation_results <- rbind(correlation_results, data.frame(
      Sample = sname,
      Condition = cond,
      Correlation = mix_correlation
    ))
    df_comp_log <- data.frame(df_comp_log)
    df_comp_log$Condition <- cond
    df_comp_log$Sample = sname
    l_res[[sname]] <- df_comp_log
}

res_df <- do.call(rbind, l_res)

# PLOTTING
if (nrow(correlation_results) > 0) {
  correlation_results$Label <- paste0("Rho: ", round(correlation_results$Correlation, 3))
  
  p <- ggplot(res_df, aes(x = Measured, y = Expected)) +
    geom_point(alpha = 0.5, size = 1) +
    theme_light() +
    labs(title = "Mixed Samples: Observed vs Synthetic Expected",
         y = "Synthetic Expected log2(Counts/CPM)",
         x = "Measured log2(Counts/CPM)") +
    theme(aspect.ratio = 1) +
    facet_wrap(~ Sample) +
    geom_text(data = correlation_results,
              aes(x = -Inf, y = Inf, label = Label),
              inherit.aes = FALSE,
              hjust = -0.1,
              vjust = 1.5,
              size = 4)
  
  ggsave(paste0(output_prefix, "_mixture_correlation.png"), p, width = 10, height = 10)
}

# 2. Summary file (RDS)
output_summary <- paste0(output_prefix, "_mixture_summary.rds")
saveRDS(res_df, output_summary)
