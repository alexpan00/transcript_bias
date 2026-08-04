library(NOISeq)

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

source_script("rep_plots.R")

args <- commandArgs(trailingOnly = TRUE)
long_obj <- args[1]
output_prefix <- args[2]


# Read long and short reads data
long_obj <- readRDS(long_obj)

#Select main condition from pData
cond <- if (ncol(pData(long_obj)) >= 2) colnames(pData(long_obj))[2] else colnames(pData(long_obj))[1]
for (plot_type in c("Expression", "Length", "GC")){
  out_plot <- paste0(output_prefix, "_replicability_", plot_type, ".png")
  p <- rep.plot(rep.dat(long_obj, factor= cond), plot_type = plot_type)
  ggsave(out_plot, p, height = 7, width = 5)
}