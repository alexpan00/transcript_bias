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

source_script("corplot.R")
source_script("lenplot.R")
source_script("pcaplot.R")
source_script("sqanti_plots.R")


args <- commandArgs(trailingOnly = TRUE)
mydata <- readRDS(args[1])
output_prefix <- args[2]
report_factors <- unlist(strsplit(args[4], ","))
analysis_type <- args[6]

# Get factors data frame
factors <- pData(mydata)

# Correlation plot
if (nrow(factors)> 1){
  p <- corplot(mydata)
} else {
  p <- ggplot2::ggplot()
}
output_corplot <- paste0(output_prefix, "_heatmap.png")
ggplot2::ggsave(filename=output_corplot, plot=p, height = 7, width = 8.5)

# create correlation matrix and save it as a table
if (inherits(mydata,"eSet") == FALSE)  stop("Error. You must give an eSet object\n")
if (!is.null(assayData(mydata)$exprs)){
  datos <- assayData(mydata)$exprs
} else {
  datos <- assayData(mydata)$counts
}
# Remove 0s
ceros = which(rowSums(datos) == 0)
hayceros = (length(ceros) > 0)
if (hayceros) {
  print(paste("Warning:", length(ceros), 
              "features with 0 counts in all samples are to be removed for this analysis."))
  datos0 = datos[-ceros,]
} else { datos0 = datos}
nsam <- NCOL(datos)
if (nsam == 1) {
  warning("You need at least 2 samples to compute correlation")
  cormat <- matrix(1, nrow = 1, ncol = 1, dimnames = list(colnames(datos), colnames(datos)))
} else {
  cormat <- cor(datos0, method = "spearman")
}
cor_table <- as.data.frame(as.table(cormat))
colnames(cor_table) <- c("Sample1", "Sample2", "Correlation")
output_cor_table <- paste0(output_prefix, "_correlation_table.tsv")
write.table(cor_table, file=output_cor_table, sep="\t", quote=FALSE, row.names=FALSE)
cat("Correlation table saved to:", output_cor_table, "\n")

for (factor in report_factors){
  # length plot 
  output_lenplot <- paste0(output_prefix, "_length_", factor, ".png")
  tryCatch({
    p <- mybias.plot(bias.dat(mydata, factor))
    ggsave(filename=output_lenplot, plot=p, height = 4, width = 6.5)
  }, error = function(e) {
    message("Warning: could not generate length bias plot for ", factor, ": ", e$message)
    ggsave(filename=output_lenplot, plot=ggplot2::ggplot() + ggplot2::ggtitle(paste("Length bias -", factor)), height = 4, width = 6.5)
  })
}

# PCA
if (nrow(pData(mydata))> 1){
  p <- pcaplot(mydata)
} else {
  p <- ggplot2::ggplot()
}
output_pcaplot <- paste0(output_prefix, "_pca.png")
ggsave(filename=output_pcaplot, plot=p, height = 3.5, width = 4.5)

for (factor in report_factors){
  #GC
  output_gcplot <- paste0(output_prefix, "_GC_", factor, ".png")
  tryCatch({
    p <- mybias.plot(bias.dat(mydata, factor, bias="GC"))
    ggsave(filename=output_gcplot, plot=p, height = 4, width = 6.5)
  }, error = function(e) {
    message("Warning: could not generate GC bias plot for ", factor, ": ", e$message)
    ggsave(filename=output_gcplot, plot=ggplot2::ggplot() + ggplot2::ggtitle(paste("GC bias -", factor)), height = 4, width = 6.5)
  })
}


# Expression vs bias wo bins
bias_v <- c("length", "GC")

for (factor in report_factors){
  for (bias in 1:2){
    for (hex in c(TRUE, FALSE)){
      scaling <- ceiling(nrow(unique(factors[factor]))/4)
      plot_type <- ifelse(hex, "hex", "point")
      bias_type <- bias_v[bias]
      output_len_exp <- paste0(output_prefix, "_", bias_type, "_", plot_type, "_", factor, ".png")
      p <- hex_bias(mydata, factor, bias, hex, bins = max(1, round(nrow(datos)/2000)))
      ggsave(filename=output_len_exp, plot=p, height = (4*scaling), width = (6.5+scaling))
    }
  }
}

if (analysis_type == "lr"){
    # ridgeline and boxplot
    for (factor in report_factors){
      for (bias_type in c("Length", "GC")){
        for (ridge in c(TRUE, FALSE)){
          scaling <- ceiling(nrow(unique(factors[factor]))/4)
          plot_type <- ifelse(ridge, "ridge", "boxplot")
          output_len_exp <- paste0(output_prefix, "_", bias_type, "_", plot_type, "_", factor, ".png")
          dat <- ridge.dat(mydata, factor = factor, bias = bias_type)
          p <- ridge.plot(dat, plot_type = plot_type)
          ggsave(filename=output_len_exp, plot=p, height = (6*scaling), width = (8.5+scaling))
        }
      }
    }

    # SQANTI3 categories plots
    ## Relative frequency of SC
    p <- sqanti_cat_relative(mydata)
    output_cat <- paste0(output_prefix, "_SQ_relative.png")
    ggsave(filename=output_cat, plot=p, height = 4, width = 6.5)

    ## Total number of SC
    p <- sqanti_cat_total(mydata)
    output_cat <- paste0(output_prefix, "_SQ_total.png")
    ggsave(filename=output_cat, plot=p, height = 4, width = 6.5)
}
