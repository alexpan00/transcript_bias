library(ggplot2)

miscolores <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")

cpm <- function(df){
  df_cpm <- data.frame(t(10^6*t(df)/colSums(as.matrix(df))))
  return(df_cpm)
}

font.size <- 2

rep.dat <- function (input, factor = NULL, norm = FALSE)  {
  
  # This plot shows the mean expression for each length bin, globally or for each biotype (if available).
  
  # datos: Count data matrix. Each column is a different biological sample.
  
  if (inherits(input,"eSet") == FALSE)    
    stop("Error. You must give an eSet object\n")
  
  if (!is.null(assayData(input)$exprs)){
    datos <- assayData(input)$exprs
  } else{
    datos <- assayData(input)$counts
  }
  ceros = which(rowSums(datos) == 0)
  
  if (length(ceros) > 0) {
    print(paste("Warning:", length(ceros), 
                "features with 0 counts in all samples are to be removed for this analysis."))
    datos = datos[-ceros,]
  }  
  
  nsam <- NCOL(datos)
  if (nsam == 1) datos <- as.matrix(datos)
  
  
  # Per condition
  if (is.null(factor)) {  # per sample  
    stop("Select a condtition to study the replicability between the samples of that condition")
    
    
  } else {  # per condition
    mifactor = as.factor(pData(input)[,factor])
    niveles = levels(mifactor)
    print("Length bias detection information is to be computed for:")
    print(niveles)
    l_df <- list()
    if (norm) {
      datos <- datos
    } else {
      datos <- cpm(datos)
    }    
    for (nivel in niveles){
      if (sum(mifactor == nivel) ==1){
        
        next
      } 
      df_level <- datos[, mifactor == nivel]
      ceros = which(rowSums(df_level) == 0)
      
      if (length(ceros) > 0) {
        print(paste("Warning:", length(ceros), 
                    "features with 0 counts in all samples are to be removed for", nivel))
        df_level = df_level[-ceros,]
      }  
      df_level <- data.frame(t(apply(df_level, 1, function(x){c("mean_exp"=sum(x)/sum(x>0), "n_samples"=sum(x>0))})))
      df_level <- merge(df_level, fData(input), by = 0)
      df_level$condition <- nivel
      df_level$n_samples <- factor(df_level$n_samples)
      l_df[[nivel]] <- df_level 
    }
  }
  final_df <- do.call(rbind, l_df)
  ## Results
  list("data2plot" = final_df, "Samples"= niveles)
}


rep.plot <- function(dat, samples = NULL, toplot = "global", toreport = FALSE, plot_type = c("Expression", "Length", "GC"))  {
  stat_box_data <- function(y) {
    return( 
      data.frame(
        y = 1.1*max(y), 
        label = paste('n =', length(y))
      )
    )
  }
  datos = dat[["data2plot"]]

  if (is.null(samples)) samples <- dat[["Samples"]]
  
  datos <- datos[datos$condition %in% samples,]
  if(length(samples) > 12) stop("Please select 12 samples or less to be plotted.")
  
  if (is.numeric(samples)) { samples = dat[["Samples"]][samples] }
  
  if (is.numeric(toplot)) {
    if (toplot == 1) { toplot = "global"} else { toplot = names(toplot)[toplot + 1] }
  }
  

  
  if ((toplot == "global") && (length(samples) <= 8)) {  ### DIAGNOSTIC PLOTS
    if (plot_type == "Expression"){
      p <- ggplot(datos, aes(n_samples, log10(mean_exp))) +
        geom_boxplot() + 
        geom_violin(alpha = 0.25) +
        stat_summary(
          fun.data = stat_box_data, 
          geom = "text", 
          hjust = 0.5,
          vjust = 0.9,
          size = font.size
        ) +
        theme_light() +
        facet_wrap(~ condition, ncol = 2)
    } else if (plot_type == "Length"){
      p <- ggplot(datos, aes(n_samples, log10(Length))) +
        geom_boxplot() + 
        geom_violin(alpha = 0.25) +
        stat_summary(
          fun.data = stat_box_data, 
          geom = "text", 
          hjust = 0.5,
          vjust = 0.9,
          size = font.size
        ) +
        theme_light() +
        facet_wrap(~ condition, ncol=2)
    } else if (plot_type == "GC"){
      p <- ggplot(datos, aes(n_samples, log10(GC))) +
        geom_boxplot() + 
        geom_violin(alpha = 0.25) +
        stat_summary(
          fun.data = stat_box_data, 
          geom = "text", 
          hjust = 0.5,
          vjust = 0.9,
          size = font.size
        ) +
        theme_light() +
        facet_wrap(~ condition, ncol = 2)
    }
    
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
  
  if((!toreport) && (length(samples) == 2)) layout(1)
  
  
}
