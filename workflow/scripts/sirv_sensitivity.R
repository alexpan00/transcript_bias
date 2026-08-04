library(NOISeq)
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)

args <- commandArgs(trailingOnly = TRUE)
long_obj <- args[1]
sirv_obj <- args[2]
output_csv <- args[3]

#show files paths
print("Input files:")
print(long_obj)
print(sirv_obj)


# Read long and short reads data
long_obj <- readRDS(long_obj)
sirv_obj <- readRDS(sirv_obj)

# create s dataframe to store sensitivity for each sample
sn_df <- data.frame(Sample = colnames(exprs(long_obj)), Detected_SIRVs = NA)

# create a selection of detected spike-ins in the long reads data, based on the SIRV object
sel <- rownames(fData(long_obj)) %in% rownames(fData(sirv_obj))

print("detection of SIRVs")
print(table(sel))

for (sample in colnames(exprs(long_obj))){
  # in the smaple column sbset the selection of sirvs and get the number of detected sirvs
  detected_sirvs <- sum(exprs(long_obj)[sel, sample] > 0)
  # create a new column in the pData with the number of detected sirvs
  sn_df[sn_df$Sample == sample, "Detected_SIRVs"] <- detected_sirvs
}

total_sirvs <- nrow(fData(sirv_obj))
sn_df <- sn_df %>% rowwise() %>% mutate(Sensitivity = if (total_sirvs > 0) Detected_SIRVs / total_sirvs else 0)
write.csv(sn_df, output_csv, row.names = FALSE)