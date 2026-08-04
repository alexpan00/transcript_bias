library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)
bam <- args[1]
gclen <- args[2]
out_name <- args[3]

# get the reference transcripts lengths
gclen <- read.table(gclen, header = F, sep = "\t", col.names = c("transcript_id", "length", "GC"))

# load the bam file
bam <- readGAlignments(bam)

# compute the median aligned length
median_len <- data.frame(bam) %>%
  group_by(seqnames) %>% 
  summarise(median_len = median(width))
median_len <- merge(median_len, gclen, by.x = "seqnames", by.y = "transcript_id")

p <- ggplot(median_len, aes(y = (median_len/length), x = log10(length))) +
  geom_hex() +
  scale_fill_viridis_c() + 
  theme_light()

ggsave(out_name, p, height = 3.5, width = 5)