cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")




sqanti_cat_relative <- function(input){
    p1 <- ggplot(data=fData(input), aes(x=Biotype)) +
    geom_bar(aes(y = (..count..)/sum(..count..)*100,fill=Biotype), color="black", size=0.3, width=0.7) +
    #geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
    scale_x_discrete(drop=FALSE) +
    xlab("") +
    ylab("Transcripts, %") +
    theme_light() +
    geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
    theme(axis.text.x = element_text(angle = 45)) +
    scale_fill_manual(values = cat.palette, guide='none') +
    theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
    scale_y_continuous(expand=expansion(mult = c(0,0.1)))

    return(p1)
}

sqanti_cat_total <- function(input){
  p1 <- ggplot(data=fData(input), aes(x=Biotype)) +
    geom_bar(aes(fill=Biotype), color="black", size=0.3, width=0.7) +
    #geom_text(aes(y = ((..count..)/sum(..count..)), label = scales::percent((..count..)/sum(..count..))), stat = "count", vjust = -0.25)  +
    scale_x_discrete(drop=FALSE) +
    xlab("") +
    ylab("Transcripts, #") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45)) +
    scale_fill_manual(values = cat.palette, guide='none') +
    theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
    scale_y_continuous(expand=expansion(mult = c(0,0.1)))
  p1
  return(p1)
}

sqanti_cat_stacked <- function(input){
  # Get expression and features
  if (!is.null(Biobase::assayData(input)$exprs)){
    counts <- Biobase::assayData(input)$exprs
  } else {
    counts <- Biobase::assayData(input)$counts
  }
  
  features <- Biobase::fData(input)
  
  # Melt for plotting, only keeping counts > 0
  plot_df <- data.frame()
  for (sample_name in colnames(counts)) {
    sample_counts <- counts[, sample_name]
    expressed_indices <- which(sample_counts > 0)
    
    if (length(expressed_indices) > 0) {
      temp_df <- data.frame(
        Sample = sample_name,
        Biotype = features$Biotype[expressed_indices]
      )
      plot_df <- rbind(plot_df, temp_df)
    }
  }
  
  # Calculate percentages per sample
  plot_summary <- plot_df %>%
    dplyr::group_by(Sample, Biotype) %>%
    dplyr::summarise(Count = dplyr::n(), .groups = 'drop') %>%
    dplyr::group_by(Sample) %>%
    dplyr::mutate(Percentage = Count / sum(Count) * 100)
    
  cond_col <- if (ncol(Biobase::pData(input)) >= 2) 2 else 1
  plot_summary$Condition <- Biobase::pData(input)[as.character(plot_summary$Sample), cond_col]
  
  p_pct <- ggplot(plot_summary, aes(x = Sample, y = Percentage, fill = Biotype)) +
    geom_bar(stat = "identity", color = "black", size = 0.3, width = 0.7) +
    xlab("") +
    ylab("Transcripts, %") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = cat.palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

  p_cnt <- ggplot(plot_summary, aes(x = Sample, y = Count, fill = Biotype)) +
    geom_bar(stat = "identity", color = "black", size = 0.3, width = 0.7) +
    xlab("") +
    ylab("Transcripts, #") +
    theme_light() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = cat.palette) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05)))
  
  return(list(plot_pct = p_pct, plot_cnt = p_cnt, data = plot_summary))
}