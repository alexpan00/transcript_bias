library(ComplexUpset)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || length(args) > 3) {
  stop("Usage: plot_ujc_upset_from_manifest.R <file_list.txt> <output.png> [output.tsv]")
}

manifest_file <- args[1]
output_plot <- args[2]
output_table <- if (length(args) == 3) args[3] else NA_character_

if (!file.exists(manifest_file)) {
  stop(paste("Manifest file not found:", manifest_file))
}

input_files <- readLines(manifest_file, warn = FALSE)
input_files <- trimws(input_files)
input_files <- input_files[nzchar(input_files)]
input_files <- input_files[!startsWith(input_files, "#")]

if (length(input_files) == 0) {
  stop(paste("No input files found in manifest:", manifest_file))
}

missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(paste("The following input files do not exist:\n", paste(missing_files, collapse = "\n")))
}

read_ujc_summary <- function(path) {
  data <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  required_columns <- c("experiment", "tool", "UJC", "n_transcripts")
  missing_columns <- setdiff(required_columns, colnames(data))
  if (length(missing_columns) > 0) {
    stop(paste("Missing columns in", path, ":", paste(missing_columns, collapse = ", ")))
  }

  data <- data[!is.na(data$UJC) & data$UJC != "", , drop = FALSE]
  data$group <- paste(data$experiment, data$tool, sep = "_")
  data
}

all_data <- do.call(rbind, lapply(input_files, read_ujc_summary))

if (nrow(all_data) == 0) {
  stop("No UJC rows found across the provided files")
}

summary_table <- aggregate(
  n_transcripts ~ experiment + tool + group + UJC,
  data = all_data,
  FUN = sum
)

summary_table <- summary_table[order(summary_table$experiment, summary_table$tool, summary_table$UJC), ]

if (!is.na(output_table)) {
  write.table(
    summary_table,
    file = output_table,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

groups <- sort(unique(summary_table$group))
ujcs <- sort(unique(summary_table$UJC))

presence <- data.frame(UJC = ujcs, stringsAsFactors = FALSE, check.names = FALSE)

for (group_name in groups) {
  presence[[group_name]] <- presence$UJC %in% summary_table$UJC[summary_table$group == group_name]
}

if (length(groups) == 0) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No UJCs found") +
    theme_void()
} else {
  p <- ComplexUpset::upset(
    presence,
    intersect = groups,
    name = "UJC",
    width_ratio = 0.2,
    themes = list(
      `Intersection size` = list(
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
        )
      ),
      intersections_matrix = list(
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()
        )
      )
    )
  ) +
    labs(title = "UJC intersections by experiment and tool")
}

ggsave(
  output_plot,
  plot = p,
  width = max(8, 1.4 * length(groups)),
  height = 6.5,
  dpi = 300
)