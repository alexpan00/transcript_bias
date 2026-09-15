library(pheatmap)
library(RColorBrewer)

# is_log_scale() lives in utils/scale_utils.R. corplot.R is always sourced by a
# script that has already defined source_script(); fall back to a direct path
# for the (unsupported) case of sourcing corplot.R on its own.
if (!exists("is_log_scale")) {
  if (exists("source_script")) {
    source_script("scale_utils.R")
  } else {
    source(file.path("scripts", "utils", "scale_utils.R"))
  }
}

miscolores <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")

cpm <- function(df){
  sums <- colSums(as.matrix(df))
  # Avoid division by zero if all counts are 0
  sums[sums == 0] <- 1
  df_cpm <- data.frame(t(10^6*t(df)/sums))
  return(df_cpm)
}

corplot <- function(mydata, method = "spearman"){
  if (inherits(mydata,"eSet") == FALSE)
    stop("Error. You must give an eSet object\n")
  
  
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
    stop("You need at least 2 samples to compute correlation")
    datos0 <- as.matrix(datos0)
  }
  
  cormat <- cor(datos0, method = method)
  cols_to_use <- intersect(c(2, 3), seq_len(ncol(Biobase::pData(mydata))))
  annot <- if (length(cols_to_use) > 0) Biobase::pData(mydata)[, cols_to_use, drop = FALSE] else Biobase::pData(mydata)[, 1, drop = FALSE]
  annot_colors <- generate_color_list(Biobase::pData(mydata))
  # This is just to avoid crashes with the testing sr-dataset
  if (all(cormat == 1)){
    cormat[1,1] <- 0.5
  }

  p <- pheatmap(cormat,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = FALSE,
           show_rownames = FALSE,
           annotation_col = annot,
           annotation_row = annot,
           annotation_colors = annot_colors,
           scale = "none",
           display_numbers = TRUE,
           fontsize = 8
  )
  return(p)
}


generate_color_list <- function(df) {
  palettes <- c("Set1", "Set2", "Set3", "Dark2", "Paired", "Pastel1", "Pastel2", "Accent")  # Different palettes
  color_list <- list()
  
  col_names <- colnames(df)[-1]  # Exclude the first column
  num_cols <- length(col_names)
  
  for (i in seq_along(col_names)) {
    col <- col_names[i]
    unique_levels <- unique(df[[col]])
    num_levels <- length(unique_levels)
    
    # Choose a palette, cycling through if there are more columns than palettes
    chosen_palette <- palettes[(i - 1) %% length(palettes) + 1]
    
    # Generate colors, falling back to rainbow if there aren't enough in the palette
    if (num_levels <= 8) {
      color_list[[col]] <- brewer.pal(max(num_levels, 3), chosen_palette)[1:num_levels]
    } else {
      color_list[[col]] <- rainbow(num_levels)
    }
    names(color_list[[col]]) <- unique_levels
  }
  
  return(color_list)
}

remove_zeros <- function(datos, condition = NULL){
  ceros = which(rowSums(datos) == 0)
  if (length(ceros) > 0 && length(ceros) != nrow(datos)) {
    if (!is.null(condition)){
      print(paste("Warning:", length(ceros), 
                  "features with 0 counts in", condition, "are to be removed for this analysis."))
    } else{
      print(paste("Warning:", length(ceros), 
                "features with 0 counts in all samples are to be removed for this analysis."))
    }
    datos = datos[-ceros,]
  }
  return(datos)
}

cor.dat <- function (input_long, input_short, factor = NULL, norm = FALSE, verbose = FALSE, full.join = FALSE)  {
  
  if (!inherits(input_long, "eSet") || !inherits(input_short, "eSet"))
    stop("Error. You must give an eSet object\n")
  
  if (!(is.null(assayData(input_long)$exprs) &is.null(assayData(input_short)$exprs))){
    datos_long <- assayData(input_long)$exprs
    datos_short <- assayData(input_short)$exprs
  } else{
    datos_long <- assayData(input_long)$counts
    datos_short <- assayData(input_short)$counts
  }
  # Remove 0s
  datos_long = remove_zeros(datos_long)
  datos_short = remove_zeros(datos_short)
  
  # Per condition
  if (is.null(factor)) {  # per sample  
    print("Correlation information is to be computed for:")
    print(colnames(datos_long))
    niveles = colnames(datos_long)
    
  } else {  # per condition
    mifactor = as.factor(pData(input_long)[,factor])
    niveles = levels(mifactor)
    print("Correlation information is to be computed for:")
    print(niveles)
    
    if (norm) {
      datos_long = sapply(niveles, 
                     function (k) {
                       rowMeans(as.matrix(datos_long[, mifactor == k]))
                     })
      datos_short = sapply(niveles, 
                     function (k) {
                       rowMeans(as.matrix(datos_short[, mifactor == k]))
                     }) 
    } else {
      datos_long = sapply(niveles, 
                     function (k) {
                       rowMeans(t(10^6*t(datos_long[, mifactor == k])/colSums(as.matrix(datos_long[, mifactor == k]))))                     
                     })
      datos_short = sapply(niveles, 
                     function (k) {
                       rowMeans(t(10^6*t(datos_short[, mifactor == k])/colSums(as.matrix(datos_short[, mifactor == k]))))                     
                     })    
    }    
    colnames(datos_long) = niveles
    colnames(datos_short) = niveles
  }
  
  ## Calculations for plot
  l_cor <- list()
  mismodelos <- list()
  for (nivel in niveles){
    tmp <- merge(remove_zeros(datos_long[, nivel, drop=F], nivel),
                            remove_zeros(datos_short[, nivel, drop = F], nivel),
                            by=0,
                            all=full.join)
    tmp[is.na(tmp)] <- 0 
    rownames(tmp) <- tmp[,1]
    tmp <- tmp[,-1]
    colnames(tmp) <- c("long", "short")
    # Put both columns on a comparable log-CPM scale. A column that is already
    # on a log scale (ratio_correction) is left untouched, but the other column
    # is still transformed -- otherwise the two axes of the regression below
    # would be on different scales and the reported R2 would be meaningless.
    long_is_log <- is_log_scale(tmp[, "long", drop = FALSE])
    short_is_log <- is_log_scale(tmp[, "short", drop = FALSE])
    if (!long_is_log && !short_is_log) {
      tmp <- remove_zeros(tmp, nivel)
      tmp <- log(cpm(tmp) + 1)
    } else {
      if (!long_is_log) {
        tmp[, "long"] <- log(cpm(tmp[, "long", drop = FALSE]) + 1)[, 1]
      }
      if (!short_is_log) {
        tmp[, "short"] <- log(cpm(tmp[, "short", drop = FALSE]) + 1)[, 1]
      }
    }
    model <- lm(short ~ long, tmp)
    l_cor[[nivel]] <- tmp
    mismodelos[[nivel]] <- model
  }
  if (verbose){
    for (nivel in niveles) {
      print(nivel)
      print(summary(mismodelos[[nivel]]))
    }
  }
  
  ## Results
  
  list("data2plot" = l_cor, "RegressionModels" = mismodelos, "Samples"= niveles)
}


mycor.plot <- function (dat, samples = NULL, toplot = "global")  {
  
  datos = dat[["data2plot"]]
  mismodelos = dat[["RegressionModels"]]
  
  if (is.null(samples)) samples <- dat[["Samples"]]
  if(length(samples) > 30) stop("Please select 30 samples or less to be plotted.")
  
  if (is.numeric(samples)) { samples = dat[["Samples"]][samples] }
  
  if (is.numeric(toplot)) {
    if (toplot == 1) { toplot = "global"} else { toplot = names(toplot)[toplot + 1] }
  }
  
  xlab <- "Long-reads CPM"
  ylab <- "Short-reads CPM"
  color_line <- "red"
  
  if ((toplot == "global") && (length(samples) <= 30)) {  ### DIAGNOSTIC PLOTS
    
    plots <- list()
    l_samples <- list()
    v_labels <- c()
    slopes <- list()
    for (i in 1:length(samples)) {
      sample_data <- datos[[i]]
      sample_data$sample <- samples[i]
      l_samples[[samples[i]]] <- sample_data
      r2 <- summary(mismodelos[[samples[i]]])$r.squared
      laF <- summary(mismodelos[[samples[i]]])$fstatistic
      p_value <- signif(pf(laF[1], df1 = laF[2], df2 = laF[3], lower.tail = FALSE), 2)
      v_labels[i] <- paste(paste("R² =", round(100 * r2, 2), "%"), paste("p-value:", p_value), sep = "\n")
      slopes[[samples[i]]] <- data.frame("intercept"=mismodelos[[samples[i]]]$coefficients[1], 
                                "slope"=mismodelos[[samples[i]]]$coefficients[2], 
                                "sample"=samples[i])
    }
    
    all_data <- do.call(rbind, l_samples)
    all_slopes <- do.call(rbind, slopes)
    p <- ggplot(all_data, aes(x = long, y=short)) +
      geom_point(size = 1) +
      geom_abline(data=all_slopes,
                  aes(intercept = intercept, slope = slope), 
                  color = color_line, 
                  linetype = "dashed", 
                  linewidth = 0.8) +
      geom_abline(aes(intercept=0, slope=1)) +
      labs(x = xlab, y = ylab) + facet_wrap(~sample, ncol=2) +
      theme_light()+ 
      theme(aspect.ratio = 1)
    p <- p + geom_text(data = data.frame(label=v_labels, sample= samples),
                       aes(x = -Inf, y = Inf, label=v_labels, group = sample),
                       inherit.aes = FALSE,
                       hjust = -1,
                       vjust = 1.2) +
      guides(color = guide_colourbar(barwidth=0.5))
    return(p)
    
  } else {  ### DESCRIPTIVE PLOTS
    
    long_data <- datos[[toplot]] %>%
      pivot_longer(cols = all_of(samples), names_to = "Sample", values_to = "Expression")
    
    ggplot(long_data, aes(x = .data[[1]], y = Expression, color = Sample)) +
      geom_line(size = 1) +
      scale_color_manual(values = miscolores[1:length(samples)]) +
      labs(title = toupper(toplot), x = "Length bins", y = "Mean expression") +
      theme_minimal() +
      theme(legend.position = "topright")
  }
  
}
