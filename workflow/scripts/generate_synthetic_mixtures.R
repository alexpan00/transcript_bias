library(Biobase)
library(NOISeq)
library(dplyr)

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
  s_dir <- get_script_dir()
  target <- file.path(s_dir, script_name)
  if (file.exists(target)) {
    source(target)
  } else if (file.exists(file.path("workflow/scripts", script_name))) {
    source(file.path("workflow/scripts", script_name))
  } else if (file.exists(file.path("scripts", script_name))) {
    source(file.path("scripts", script_name))
  } else {
    source(script_name)
  }
}

source_script("norm_methods.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript generate_synthetic_mixtures.R <raw_noiseq_obj.rds> <mixture_definition.csv> <output_rds>")
}

raw_obj_path <- args[1]
mix_def_path <- args[2]
output_path <- args[3]

# Read data
raw_obj <- readRDS(raw_obj_path)
mix_def <- read.csv(mix_def_path)

# Ensure required columns in mix_def
required_cols <- c("Mixture", "Component", "Proportion")
if (!all(required_cols %in% colnames(mix_def))) {
  stop("Mixture definition file must contain columns: ", paste(required_cols, collapse = ", "))
}

# Get raw count data
data_raw <- exprs(raw_obj)
data_cpm <- cpm(data_raw)

# Get metadata
metadata <- pData(raw_obj)
cond_col <- colnames(metadata)[2]
sample_col <- colnames(metadata)[1]

# Retrieve sample names per condition
sample_list <- split(metadata[[sample_col]], metadata[[cond_col]])

# Keep the pure component samples in the final matrix so normalization methods
# see the same full sample space as the original object.
pure_profiles <- list()
unique_pure_conditions <- unique(mix_def$Component)

for (pure_cond in unique_pure_conditions) {
  pure_samples <- sample_list[[pure_cond]]

  if (is.null(pure_samples) || length(pure_samples) == 0) {
    warning("No samples found for pure component condition in metadata: ", pure_cond)
    next
  }

  for (pure_sample in pure_samples) {
    pure_profiles[[as.character(pure_sample)]] <- data_raw[, pure_sample]
  }
}

# 1. GENERATE SYNTHETIC MIXTURES
synthetic_profiles <- list()
unique_mixtures <- unique(mix_def$Mixture)

for (mix_cond in unique_mixtures) {
  components <- mix_def %>% filter(Mixture == mix_cond)
  mix_samples <- sample_list[[mix_cond]]
  
  if (is.null(mix_samples) || length(mix_samples) == 0) {
    warning("No samples found for mixture condition in metadata: ", mix_cond)
    next
  }

  for (r_idx in seq_along(mix_samples)) {
    curr_mix_sample <- mix_samples[r_idx]
    
    expected_profile_cpm <- numeric(nrow(data_raw))
    names(expected_profile_cpm) <- rownames(data_raw)
    valid_replicate <- TRUE
    
    for (i in 1:nrow(components)) {
      comp_cond <- components$Component[i]
      prop <- components$Proportion[i]
      comp_samples <- sample_list[[comp_cond]]
      
      if (is.null(comp_samples) || r_idx > length(comp_samples)) {
        warning(paste("Missing replicate", r_idx, "for component", comp_cond))
        valid_replicate <- FALSE
        break
      }
      
      curr_comp_sample <- comp_samples[r_idx]
      expected_profile_cpm <- expected_profile_cpm + (data_cpm[, curr_comp_sample] * prop)
    }
    
    if (valid_replicate && !all(expected_profile_cpm == 0)) {
        lib_size <- sum(data_raw[, curr_mix_sample])
        raw_counts <- (expected_profile_cpm * lib_size) / 1e6
        synthetic_profiles[[as.character(curr_mix_sample)]] <- raw_counts
    }
  }
}

combined_profiles <- c(pure_profiles, synthetic_profiles)

if (length(combined_profiles) == 0) {
  stop("No valid expected profiles could be generated.")
}

# 2. Create combined matrix and ensure column names
synthetic_matrix <- as.matrix(do.call(cbind, combined_profiles))
if (ncol(synthetic_matrix) != length(combined_profiles)) {
    # This happens if there's only one sample and cbind makes it a vector
  synthetic_matrix <- matrix(combined_profiles[[1]], ncol=1)
  colnames(synthetic_matrix) <- names(combined_profiles)
}
# Double check names
if (is.null(colnames(synthetic_matrix))) {
  colnames(synthetic_matrix) <- names(combined_profiles)
}

factors_synthetic <- data.frame(
    sample = colnames(synthetic_matrix),
  condition = metadata[[cond_col]][match(colnames(synthetic_matrix), metadata[[sample_col]])],
    stringsAsFactors = FALSE
)
rownames(factors_synthetic) <- factors_synthetic$sample

# Extract feature data
mylength <- fData(raw_obj)$Length
names(mylength) <- rownames(fData(raw_obj))
mygc <- fData(raw_obj)$GC
names(mygc) <- rownames(fData(raw_obj))
mybiotype <- fData(raw_obj)$Biotype
names(mybiotype) <- rownames(fData(raw_obj))

# Biotype needs to be a data frame for readData
biotype_df <- data.frame(category = mybiotype, row.names = names(mybiotype))

new_noiseq <- readData(
    data = synthetic_matrix,
    factors = factors_synthetic,
    length = mylength,
    gc = mygc,
    biotype = biotype_df
)

saveRDS(new_noiseq, output_path)
