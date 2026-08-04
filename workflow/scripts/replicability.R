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