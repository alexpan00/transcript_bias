library(NOISeq)
library(ggplot2)
library(tidyr)
library(dplyr)


cpm <- function(df){
  df_cpm <- data.frame(t(10^6*t(df)/colSums(as.matrix(df))))
  return(df_cpm)
}

args <- commandArgs(trailingOnly = TRUE)
long_obj <- args[1]
class_file <- args[2]
output_prefix <- args[3]

# Read long expression and sqanti classification
long_obj <- readRDS(long_obj)
class_file <- read.delim(class_file)

# select relevant columns from the classification file to build transcript to
# gene relation
t2g <- class_file %>% 
  select(isoform, associated_gene) %>% 
  filter(!grepl("novel", .data$associated_gene))

# build the expression data
counts <- cpm(exprs(long_obj))
counts$isoform <- rownames(counts)

# add the gene
counts <- inner_join(counts, t2g, by = "isoform")

# remove trasncripts from samples in which they are not expressed
counts_long <- counts %>% 
  pivot_longer(!c(isoform, associated_gene), names_to = "sample", values_to = "Counts")
# add tissue/condition
cond <- if (ncol(pData(long_obj)) >= 2) colnames(pData(long_obj))[2] else colnames(pData(long_obj))[1]
pData_df <- pData(long_obj)
if (!"sample" %in% colnames(pData_df)) {
  pData_df$sample <- rownames(pData_df)
}
counts_long <- full_join(counts_long, pData_df[, unique(c("sample", cond)), drop = FALSE], by = "sample")

# mean by condition
mean_transcript_counts_long <- counts_long %>% 
  group_by(associated_gene, isoform, .data[[cond]]) %>% 
  summarise(mean_counts = mean(Counts)) 

# gene expression and number of isoforms per gene
gene_summary <- mean_transcript_counts_long %>% 
  filter(mean_counts > 0) %>% 
  group_by(.data[[cond]], associated_gene) %>% 
  summarise(mean_isoform_exp = mean(mean_counts),
            max_isoform_exp = max(mean_counts),
            total_gene_exp = sum(mean_counts),
            total_isoforms = n()) %>% 
  group_by(.data[[cond]], total_isoforms) %>% 
  summarise(median_mean_isoform_exp = median(mean_isoform_exp),
            sd_mean_isoform_exp = sd(mean_isoform_exp),
            median_total_gene_exp = median(total_gene_exp),
            sd_total_gene_exp = sd(total_gene_exp),
            median_max_isoform_exp = median(max_isoform_exp))

gene_summary[is.na(gene_summary)] <- 0

# Plot
p <- ggplot(gene_summary, aes(x= total_isoforms, y = log10(median_total_gene_exp))) +
  geom_line() +
  geom_point() +
  facet_wrap(~.data[[cond]], scales = "free_x") +
  xlab("Number of isoforms") +
  ylab("Median total gene expression") +
  theme_light()
output_box <- paste0(output_prefix, "_gene_expression_vs_n_isoforms.png")
ggsave(filename = output_box, plot = p, height = 4, width = 6.5)


p <- ggplot(gene_summary, aes(x = total_isoforms, y = log10(median_mean_isoform_exp))) +
  geom_line() +
  geom_point() +
  facet_wrap(~.data[[cond]], scales = "free_x") +
  xlab("Number of isoforms") +
  ylab("Median mean isoform expression") +
  theme_light()
output_box <- paste0(output_prefix, "_mean_isoform_expression_vs_n_isoforms.png")
ggsave(filename = output_box, plot = p, height = 4, width = 6.5)


p <- ggplot(gene_summary, aes(x = total_isoforms, y = log10(median_max_isoform_exp))) +
  geom_line() +
  geom_point() +
  facet_wrap(~.data[[cond]], scales = "free_x") +
  xlab("Number of isoforms") +
  ylab("Median max isoform expression") +
  theme_light()
output_box <- paste0(output_prefix, "_max_isoform_expression_vs_n_isoforms.png")
ggsave(filename = output_box, plot = p, height = 4, width = 6.5)