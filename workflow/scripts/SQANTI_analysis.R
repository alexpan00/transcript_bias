library(NOISeq)
library(tidyverse)

xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")

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

source_script("sqanti_plots.R")


args <- commandArgs(trailingOnly = TRUE)
mydata <- readRDS(args[1])
output_prefix <- args[2]

# SQANTI3 categories plots
## Stacked frequency of SC per sample
res <- sqanti_cat_stacked(mydata)
output_cat_stacked_pct <- paste0(output_prefix, "_SQ_stacked_pct.png")
ggsave(filename = output_cat_stacked_pct, plot = res$plot_pct, height = 6, width = 8)

## Stacked total isoforms of SC per sample
output_cat_stacked_cnt <- paste0(output_prefix, "_SQ_stacked_cnt.png")
ggsave(filename = output_cat_stacked_cnt, plot = res$plot_cnt, height = 6, width = 8)

## Save data for combined plot
output_csv <- paste0(output_prefix, "_SQ_counts_per_sample.csv")
write.csv(res$data, file = output_csv, row.names = FALSE)

# Save total counts of each SQANTI3 category across all samples
global_summary <- fData(mydata) %>% 
    group_by(Biotype) %>% 
    summarise(n_trans = n())

write.csv(global_summary, file = paste0(output_prefix, "_SQ_counts_global.csv"), row.names = FALSE)