library(tidyverse)
xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")

# Arguments are:
# 1. Comma separated list of SQANTI global-count CSVs (one per tool)
# 2. Comma separated list of structural-category TSVs for the prefiltered transcriptomes (one per tool)
# 3. Output summary table path
args <- commandArgs(trailingOnly = TRUE)
sqanti_csvs <- args[1]
prefiltered_categories <- args[2]
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
  sq_data$Filtered <- TRUE

  l_res[[csv_file]] <- sq_data
}

categories_prefiltered <- unlist(strsplit(prefiltered_categories, ","))

# For each CSV file, get tool info from path and add it
for (delim_file in categories_prefiltered){
  # get tool from the path (e.g., results/bambu/NOIseq/...)
  meta <- parse_path_metadata(delim_file)
  tool <- meta$tool
  
  # Read structural categories data (isoform	structural_category)
  sq_data <- read.table(delim_file, header = TRUE, sep = "\t")
  sq_data$structural_category <- factor(sq_data$structural_category, levels = xaxislevelsF1, labels=xaxislabelsF1)

  sq_data_summary <- sq_data %>% 
    group_by(structural_category) %>% 
    summarise(n_trans = n(), .groups = "drop") %>%
    mutate(Tool = tool)

  # replace the structural_category coname by Biotype to match the other summary
  colnames(sq_data_summary)[colnames(sq_data_summary) == "structural_category"] <- "Biotype"
  sq_data_summary$Filtered <- FALSE
  l_res[[delim_file]] <- sq_data_summary
}

# Combine results
final_df <- do.call(rbind, l_res)

# Write full combined table
write.csv(final_df, output_summary_table, row.names = FALSE)
