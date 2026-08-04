library(NOISeq)

args <- commandArgs(trailingOnly = TRUE)
noiseq_object <- args[1]
output_file <- args[2]

mydata <- readRDS(noiseq_object)

if (!inherits(mydata, "eSet")) {
  stop("Input object must be an eSet/NOISeq object")
}

transcript_ids <- rownames(Biobase::fData(mydata))
if (is.null(transcript_ids) || length(transcript_ids) == 0) {
  transcript_ids <- rownames(Biobase::exprs(mydata))
}

transcript_ids <- transcript_ids[!is.na(transcript_ids) & transcript_ids != ""]

writeLines(transcript_ids, con = output_file)