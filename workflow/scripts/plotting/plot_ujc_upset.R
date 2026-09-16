library(ComplexUpset)
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: plot_ujc_upset.R <input1.tsv> <input2.tsv> ... <output.png> <output.tsv>")
}

input_files <- args[1:(length(args) - 2)]
output_plot <- args[length(args) - 1]
output_table <- args[length(args)]

tool_name_from_path <- function(path) {
  basename(dirname(dirname(path)))
}

experiment_name_from_path <- function(path) {
  sub("_UJC\\.tsv$", "", basename(path))
}

read_ujc_set <- function(path) {
  data <- read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (!"UJC" %in% colnames(data)) {
    stop(paste("Missing UJC column in", path))
  }
  unique(na.omit(data$UJC))
}

tool_names <- vapply(input_files, tool_name_from_path, character(1))
ujc_sets <- setNames(lapply(input_files, read_ujc_set), tool_names)

ujc_long <- do.call(
  rbind,
  lapply(seq_along(input_files), function(idx) {
    tool <- tool_names[[idx]]
    data <- read.delim(input_files[[idx]], sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    if (!"UJC" %in% colnames(data)) {
      stop(paste("Missing UJC column in", input_files[[idx]]))
    }

    # A tool whose transcripts are all mono-exonic yields a header-only file.
    # Recycling a scalar against a zero-length column would error, so build the
    # empty frame explicitly.
    if (nrow(data) == 0) {
      return(data.frame(experiment = character(0), tool = character(0),
                        transcript_id = character(0), UJC = character(0),
                        stringsAsFactors = FALSE))
    }

    data.frame(
      experiment = experiment_name_from_path(input_files[[idx]]),
      tool = tool,
      transcript_id = data$transcript_id,
      UJC = data$UJC,
      stringsAsFactors = FALSE
    )
  })
)

ujc_long <- ujc_long[!is.na(ujc_long$UJC) & ujc_long$UJC != "", ]
ujc_summary <- aggregate(
  transcript_id ~ experiment + tool + UJC,
  data = ujc_long,
  FUN = length
)
colnames(ujc_summary)[colnames(ujc_summary) == "transcript_id"] <- "n_transcripts"

write.table(
  ujc_summary,
  file = output_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

all_ujcs <- sort(unique(unlist(ujc_sets, use.names = FALSE)))

if (length(all_ujcs) == 0) {
  p <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No UJCs found across tools") +
    theme_void()
} else {
  presence <- data.frame(UJC = all_ujcs, stringsAsFactors = FALSE, check.names = FALSE)
  for (tool in names(ujc_sets)) {
    presence[[tool]] <- presence$UJC %in% ujc_sets[[tool]]
  }

  p <- ComplexUpset::upset(
    presence,
    intersect = names(ujc_sets),
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
    labs(title = "UJC intersections across tools")
}

ggsave(output_plot, plot = p, width = max(8, 1.4 * length(tool_names)), height = 6.5, dpi = 300)