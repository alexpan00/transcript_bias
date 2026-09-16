library(ggridges)
library(dplyr)
library(tidyr)
library(ggplot2)
library(splines)

miscolores <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A")

## Data for gene length plot
remove_zeros <- function(datos, condition = NULL){
  ceros = which(rowSums(datos) == 0)
  if (length(ceros) > 0) {
    if (!is.null(condition)){
      print(paste("Warning:", length(ceros), 
                  "features with 0 counts in", condition, "are to be removed for this analysis."))
    } else{
      print(paste("Warning:", length(ceros), 
                  "features with 0 counts in all samples are to be removed for this analysis."))
    }
    datos = datos[-ceros,, drop=F]
  }
  return(datos)
}
trimmed_sd <- function(x, trim = 0.025) {
  n <- length(x)
  if (n == 0) return(NA)
  if (n == 1) return(0)
  if (n == 2) return(abs(diff(x)))
  
  x <- sort(x)
  ntrim <- floor(n * trim)
  x <- x[(ntrim + 1):(n - ntrim)]
  
  sd(x)
}

trimmed_cv <- function(x, trim = 0.025) {
  n <- length(x)
  if (n == 0) return(NA)
  if (n == 1) return(0)
  if (n == 2) return(abs(diff(x)))
  
  x <- sort(x)
  ntrim <- floor(n * trim)
  x <- x[(ntrim + 1):(n - ntrim)]
  
  sd(x) / abs(mean(x))
}

get_lims <- function(x){
  Q1 <- quantile(x, 0.25)
  Q3 <- quantile(x, 0.75)
  IQR <- Q3 - Q1

  lower_limit <- Q1 - 1.5 * IQR
  upper_limit <- Q3 + 1.5 * IQR

  return(c(lower_limit, upper_limit))
}

bias.dat <- function (input, factor = NULL, norm = FALSE, numXbin=200, verbose = FALSE, bias = "Length")  {
  
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
    print("Length bias detection information is to be computed for:")
    print(colnames(datos))
    mifactor = as.factor(colnames(datos))
    niveles = levels(mifactor)
    
  } else {  # per condition
    mifactor = as.factor(pData(input)[,factor])
    niveles = levels(mifactor)
    print("Length bias detection information is to be computed for:")
    print(niveles)
    
    if (norm) {
      datos = sapply(niveles, 
                     function (k) {
                       rowMeans(as.matrix(datos[, mifactor == k]))
                     })    
    } else {
      datos = sapply(niveles, 
                     function (k) {
                       rowMeans(t(10^6*t(datos[, mifactor == k])/colSums(as.matrix(datos[, mifactor == k]))))                     
                     })    
    }    
    colnames(datos) = niveles
  }
  
  
  
  # Length or GC content
  if (any(!is.na(fData(input)[,bias])) == FALSE)
    stop ("Feature length was not provided.\nPlease run addData() function to add 
          this information\n")
  
  long <- as.numeric(as.character(fData(input)[,bias]))
  names(long) <- rownames(fData(input))
  if (length(ceros) > 0) long = long[-ceros] 
  
  
  
  # Biotypes
  if (!is.null(featureData(input)$Biotype)) {  # read biotypes if they are provided
    infobio <- as.character(featureData(input)$Biotype)
    if (length(ceros) > 0) infobio = infobio[-ceros] 
    biotypes <- unique(infobio)
    names(biotypes) <- biotypes 
    # which genes belong to each biotype
    biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })
    names(biog) = biotypes
    bionum <- c(NROW(datos), sapply(biog, length))
    names(bionum) <- c("global", names(biotypes))   
    
  } else { infobio = NULL; biotypes = NULL; bionum = NULL }  
  
  
  
  ## Calculations for plot
  
  longexpr = vector("list", length = 1 + length(biotypes))
  names(longexpr) = c("global", names(biotypes))
  
  
  for (i in 1:length(longexpr))  {
    
    if (i == 1) {  # GLOBAL
      l_cond <- list()
      for (nivel in niveles){
        datos_cond <- datos[,nivel, drop = F]
        datos_cond <- remove_zeros(datos_cond, nivel)
        lengths_cond <- long[rownames(datos_cond)]
        numdatos = length(lengths_cond)      
        numbins = floor(numdatos / numXbin)
        misbins = quantile(lengths_cond, probs = seq(0,1,1/numbins), na.rm = TRUE)
        
        if (length(misbins) != length(unique(misbins))) {
          repes = names(table(misbins))[which(table(misbins) > 1)]
          for (rr in repes) {
            cuantos = length(which(misbins == rr))
            cuales = which(misbins == rr)
            sumo = (misbins[cuales[1]+cuantos] - misbins[cuales[1]])/cuantos
            for (j in cuales[-1]) misbins[j] = misbins[j-1] + sumo            
          }
        }
        misbins[1] <- misbins[1]-1
        miclasi = cut(lengths_cond, breaks = misbins, labels = FALSE)
        #misbins = sapply(1:numbins, function (i) mean(misbins[i:(i+1)]))
        l_class <- cbind(datos_cond,data.frame(lengths_cond, miclasi))
        l_cond[[nivel]] <- l_class %>% 
          group_by(miclasi) %>% 
          summarise(
            lengthbin = median(lengths_cond),
            across(all_of(nivel), list(mean = ~mean(.x, trim = 0.025), cv = ~trimmed_cv(.x, trim = 0.025)), .names = "{.col}_{.fn}")
          ) %>% 
          dplyr::select(-miclasi)
      }
      longexpr[[i]] <- l_cond
      
    } else {  # PER BIOTYPE
      
      datos2 = datos[biog[[i-1]],]
      long2 = long[biog[[i-1]]]
      
      if (bionum[i] >= numXbin*10) {  # more than numXbin*10 genes in the biotype
        
        numdatos = length(long2)      
        numbins = floor(numdatos / numXbin)
        misbins = quantile(long2, probs = seq(0,1,1/numbins), na.rm = TRUE)
        
        if (length(misbins) != length(unique(misbins))) {
          repes = names(table(misbins))[which(table(misbins) > 1)]
          for (rr in repes) {
            cuantos = length(which(misbins == rr))
            cuales = which(misbins == rr)
            sumo = (misbins[cuales[1]+cuantos] - misbins[cuales[1]])/cuantos
            for (j in cuales[-1]) misbins[j] = misbins[j-1] + sumo            
          }
        }
        
        miclasi = cut(long2, breaks = misbins, labels = FALSE)
        misbins = sapply(1:numbins, function (i) mean(misbins[i:(i+1)]))
        miclasi = misbins[miclasi]        
        longexpr[[i]] = aggregate(datos2, by = list("lengthbin" = miclasi), mean, trim = 0.025, na.rm = TRUE)        
        
      } else {   # less than numXbin*10 genes in the biotype
        
        longexpr[[i]] = cbind(long2, datos2)      
        
      }      
      
    }
  }
  
  
  
  ## SPLINES REGRESSION MODEL   
  library(splines)
  
  datos_list = longexpr[[1]]
  mismodelos = vector("list", length = (ncol(datos)))
  names(mismodelos) = niveles
  
  for (i in 1:length(datos_list)) {
    longi = if (!is.null(datos_list[[i]]$lengthbin)) datos_list[[i]]$lengthbin else datos_list[[i]][[1]]
    if (is.null(longi) || length(longi) < 2) next
    n_knots = max(1, round(length(longi)/10, 0))
    knots =  c(rep(longi[1],3), seq(longi[1], longi[length(longi)-1], length.out=n_knots), 
               rep(longi[length(longi)], 4))
    bx = splineDesign (knots, longi, outer.ok = TRUE)
    
    print(colnames(datos)[i])
    datos_vector <- pull(datos_list[[i]][2]) # expression of bins in a condition
    mismodelos[[i]] = lm(datos_vector ~ bx)
    
    if (verbose){
      print(summary(mismodelos[[i]]))
    }
  }
  
  
  
  
  ## Results
  
  list("data2plot" = longexpr, "RegressionModels" = mismodelos, "Samples"= niveles, "bias"= bias)
}


mybias.plot <- function (dat, samples = NULL, toplot = "global", toreport = FALSE, ylim= NULL,...)  {
  
  datos = dat[["data2plot"]]
  mismodelos = dat[["RegressionModels"]]
  
  if (is.null(samples)) samples <- dat[["Samples"]]
  if(length(samples) > 12) stop("Please select 12 samples or less to be plotted.")
  
  if (is.numeric(samples)) { samples = dat[["Samples"]][samples] }
  
  if (is.numeric(toplot)) {
    if (toplot == 1) { toplot = "global"} else { toplot = names(toplot)[toplot + 1] }
  }
  
  xlab <- if (dat[["bias"]] == "Length"){"Length bins"} else {"GC% bins"}
  color_line <- if (dat[["bias"]] == "Length"){"red"} else {"blue"}
  
  if ((toplot == "global") && (length(samples) <= 24)) {  ### DIAGNOSTIC PLOTS
    
    plots <- list()
    l_samples <- list()
    v_labels <- c()
    for (i in 1:length(samples)) {
      if (is.null(mismodelos[[samples[i]]])) {
        sample_data <- datos[[1]][[samples[i]]] %>%
          dplyr::select(Length = 1, Expression = 2, CV =3)%>%
          mutate(ModelFit = Expression, sample=samples[i])
        l_samples[[samples[i]]] <- sample_data
        v_labels[i] <- "N/A"
      } else {
        sample_data <- datos[[1]][[samples[i]]] %>%
          dplyr::select(Length = 1, Expression = 2, CV =3)%>%
          mutate(ModelFit = mismodelos[[samples[i]]]$fit, sample=samples[i])
        l_samples[[samples[i]]] <- sample_data
        r2 <- summary(mismodelos[[samples[i]]])$r.squared
        laF <- summary(mismodelos[[samples[i]]])$fstatistic
        p_value <- if (!is.null(laF)) signif(pf(laF[1], df1 = laF[2], df2 = laF[3], lower.tail = FALSE), 2) else "N/A"
        r2_val <- if (!is.null(r2)) round(100 * r2, 2) else 0
        v_labels[i] <- paste(paste("R² =", r2_val, "%"), paste("p-value:", p_value), sep = "\n")
      }
    }
    
    all_data <- do.call(rbind, l_samples)
    x_max <- get_lims(all_data$Length)[2]
    y_max <- get_lims(all_data$Expression)[2]
    p <- ggplot(all_data, aes(x = Length)) +
      geom_point(aes(y = Expression, color = CV), size = 1) +
      geom_line(aes(y = ModelFit), color = color_line, linetype = "dashed", size = 0.8) +
      labs(x = xlab, y = "Mean expression") + facet_wrap(~sample, ncol=2) +
      theme_light() + scale_color_viridis_c(option = "plasma") + 
      guides(color = guide_colourbar(barwidth=0.5)) + ylim(c(0,y_max))
    p <- p + geom_text(data = data.frame(label=v_labels, sample= samples),
                  aes(x = -Inf, y = Inf, label=v_labels, group = sample),
                  inherit.aes = FALSE,
                  hjust = -1,
                  vjust = 1.2,
                  size = 2)
    if (dat[["bias"]]  == "Length") {
      p <- p + xlim(c(0, x_max))
    }
    #results <- list(p = p, data = all_data)
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

hex_bias <- function(mydata, factor, bias = 1, hex = T, bins = 15, alpha = .25){
  # if hex is F plot all the points with alpha
  # bias 1 for length 2 for GC content
  # bins number of hexagons to divide the plot
  
  # selct the bias
  bias <- c("Length", "GC")[bias]
  
  # condtition and transcript info
  mifactor = as.factor(pData(mydata)[,factor])
  gclength <- fData(mydata)
  
  l_exp <- list()
  for (condition in unique(mifactor)) {
    # Normalize expresion(CPM) of the slected samples
    mean_expr <- rowMeans(rpkm(exprs(mydata)[, pData(mydata)[[factor]] == condition, drop = FALSE], long = 1000))
    
    # Combine into a data frame with length and GC
    df <- data.frame(
      mean_expr = mean_expr,
      bias = gclength[,bias],
      condition = condition
    )
    
    # Filter between 5th and 95th percentiles
    df_filtered <- df %>%
      filter(mean_expr != 0) %>% 
      filter(mean_expr > quantile(mean_expr, 0.05) & 
               mean_expr < quantile(mean_expr, 0.95)) %>% 
        filter(bias > quantile(bias, 0.05) &
               bias < quantile(bias, 0.95))
    if (nrow(df_filtered) == 0) {
      df_filtered <- df
    }
    
    
    # Save result to list
    l_exp[[condition]] <- df_filtered
  }
  df_exp <- bind_rows(l_exp)
  rm(l_exp)
  
  # wether to do point or hex plot
  if (hex == F){
    p <- ggplot(df_exp, aes(x = bias, y = mean_expr)) +
      geom_point(alpha = alpha) +
      facet_wrap(~condition, ncol = 2) +
      xlab(bias) +
      theme_light() 
  } else {
    p <- ggplot(df_exp, aes(x = bias, y = mean_expr)) +
      geom_hex(bins = bins) +
      xlab(bias) +
      scale_fill_viridis_c() + 
      facet_wrap(~condition, ncol = 2) +
      theme_light()
  }
  return(p)
}


ridge.dat <- function (input, factor = NULL, norm = FALSE, verbose = FALSE, bias = "Length", n_cuts = 15)  {
  
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
    print("Length bias detection information is to be computed for:")
    print(colnames(datos))
    
  } else {  # per condition
    mifactor = as.factor(pData(input)[,factor])
    niveles = levels(mifactor)
    print("Length bias detection information is to be computed for:")
    print(niveles)
    # sel <- sapply(niveles, 
    #               function(k){
    #                 apply(datos[, mifactor == k, drop =F], 1, function(x){all(x != 0)})
    #           })
    if (norm) {
      datos = sapply(niveles, 
                     function (k) {
                       rowMeans(as.matrix(datos[, mifactor == k]))
                     })    
    } else {
      datos = sapply(niveles, 
                     function (k) {
                       rowMeans(t(10^6*t(datos[, mifactor == k])/colSums(as.matrix(datos[, mifactor == k]))))                     
                     })    
    }    
    colnames(datos) = niveles
  }
  
  
  
  # Length or GC content
  if (any(!is.na(fData(input)[,bias])) == FALSE)
    stop ("Feature length was not provided.\nPlease run addData() function to add 
          this information\n")
  
  long <- as.numeric(as.character(fData(input)[,bias]))
  names(long) <- rownames(fData(input))
  if (length(ceros) > 0) long = long[-ceros] 
  
  
  
  # Biotypes
  if (!is.null(featureData(input)$Biotype)) {  # read biotypes if they are provided
    infobio <- as.character(featureData(input)$Biotype)
    if (length(ceros) > 0) infobio = infobio[-ceros] 
    biotypes <- unique(infobio)
    names(biotypes) <- biotypes 
    # which genes belong to each biotype
    biog <- lapply(biotypes, function(x) { which(is.element(infobio, x)) })
    names(biog) = biotypes
    bionum <- c(NROW(datos), sapply(biog, length))
    names(bionum) <- c("global", names(biotypes))   
    
  } else { infobio = NULL; biotypes = NULL; bionum = NULL }  
  
  
  
  ## Calculations for plot
  
  longexpr = vector("list", length = 1 + length(biotypes))
  names(longexpr) = c("global", names(biotypes))
  
  
  for (i in 1:length(longexpr))  {
    
    if (i == 1) {  # GLOBAL
      l_cond <- list()
      
      breaks <- quantile(long, probs = seq(0,1,1/n_cuts), na.rm = TRUE)
      breaks[1] <- breaks[1]-1
      for (nivel in niveles){
        datos_cond <- datos[,nivel, drop = F]
        datos_cond <- remove_zeros(datos_cond, nivel)
        #datos_cond <- datos_cond[sel[,nivel],,drop=F]
        lengths_cond <- long[rownames(datos_cond)]
        l_bins <- cut(lengths_cond, breaks = breaks)
        l_class <- cbind(datos_cond,data.frame(lengths_cond, l_bins))
        l_cond[[nivel]] <- l_class
      }
      longexpr[[i]] <- l_cond
      
    } else {  # PER BIOTYPE
      
      datos2 = datos[biog[[i-1]],]
      long2 = long[biog[[i-1]]]
      
      if (bionum[i] >= 200*10) {  # more than numXbin*10 genes in the biotype
        
        numdatos = length(long2)      
        numbins = floor(numdatos / 200)
        misbins = quantile(long2, probs = seq(0,1,1/numbins), na.rm = TRUE)
        
        if (length(misbins) != length(unique(misbins))) {
          repes = names(table(misbins))[which(table(misbins) > 1)]
          for (rr in repes) {
            cuantos = length(which(misbins == rr))
            cuales = which(misbins == rr)
            sumo = (misbins[cuales[1]+cuantos] - misbins[cuales[1]])/cuantos
            for (j in cuales[-1]) misbins[j] = misbins[j-1] + sumo            
          }
        }
        
        miclasi = cut(long2, breaks = misbins, labels = FALSE)
        misbins = sapply(1:numbins, function (i) mean(misbins[i:(i+1)]))
        miclasi = misbins[miclasi]        
        longexpr[[i]] = aggregate(datos2, by = list("lengthbin" = miclasi), mean, trim = 0.025, na.rm = TRUE)        
        
      } else {   # less than numXbin*10 genes in the biotype
        
        longexpr[[i]] = cbind(long2, datos2)      
        
      }      
      
    }
  }
  
  datos_list = longexpr[[1]]
  mismodelos = vector("list", length = (ncol(datos)))
  names(mismodelos) = niveles
  
  for (i in 1:length(datos_list)) {
    longi = datos_list[[i]]$l_bins
    datos_vector <- pull(datos_list[[i]][1]) # expression of bins in a condition
    l_models <- list()
    l_models[["kw"]] <- kruskal.test(datos_vector ~ longi)
    l_models[["fl"]] <- fligner.test(datos_vector ~ longi)
    mismodelos[[i]] <- l_models
    
    if (verbose){
      print(names(mismodelos)[i])
      print((mismodelos[[i]][[1]]))
      print((mismodelos[[i]][[2]]))
    }
  }
  list("data2plot" = longexpr, "Samples"= niveles, "bias"= bias, "RegressionModels" = mismodelos)
}


ridge.plot <- function (dat,plot_type="ridge", samples = NULL, toplot = "global", toreport = FALSE, ylim= NULL,...)  {
  stat_box_data <- function(y) {
    return( 
      data.frame(
        y = 1.1*max(y), 
        label = paste('n =', length(y))
      )
    )
  }
  
  datos = dat[["data2plot"]]
  mismodelos <- dat[["RegressionModels"]]
  if (is.null(samples)) samples <- dat[["Samples"]]
  if(length(samples) > 12) stop("Please select 12 samples or less to be plotted.")
  
  if (is.numeric(samples)) { samples = dat[["Samples"]][samples] }
  
  if (is.numeric(toplot)) {
    if (toplot == 1) { toplot = "global"} else { toplot = names(toplot)[toplot + 1] }
  }
  
  ylab <- if (dat[["bias"]] == "Length"){"Length bins"} else {"GC% bins"}
  color_line <- if (dat[["bias"]] == "Length"){"red"} else {"blue"}
  
  if ((toplot == "global") && (length(samples) <= 8)) {  ### DIAGNOSTIC PLOTS
    
    l_samples <- list()
    v_labels <- c()
    for (i in 1:length(samples)) {
      sample_data <- datos[[1]][[i]] %>%
        select(Expression = 1, Length = 3)%>%
        mutate(sample=samples[i])
      l_samples[[samples[i]]] <- sample_data
      KW_pval <- mismodelos[[samples[i]]][["kw"]]$p.value
      FK_pval <- mismodelos[[samples[i]]][["fl"]]$p.value
      v_labels[i] <- paste(paste("KW p-val:", round(KW_pval, 2)), paste("FK p-val:", round(FK_pval,2)), sep = "\n")
      l_samples[[samples[i]]] <- sample_data
    }
    
    all_data <- do.call(rbind, l_samples)
    intervals <- as.character(all_data$Length)
    lower_bounds <- as.numeric(sub("\\((.+),.*", "\\1", intervals))
    sorted_intervals <- unique(intervals[order(lower_bounds)])
    all_data$Length <- factor(all_data$Length, levels = sorted_intervals)
    all_data <- all_data %>% 
      group_by(sample, Length) %>% 
      mutate(n = n())
    medians <- all_data %>% 
      group_by(sample) %>% 
      summarise(log10_medians = log10(median(Expression)))
    if (plot_type == "boxplot"){
      p <- ggplot(all_data, aes(x= log10(Expression), y = Length)) +
        geom_boxplot(aes(color = n)) +
        geom_violin(alpha= 0.25)+
        stat_summary(
          fun.data = stat_box_data, 
          geom = "text", 
          hjust = 0.5,
          vjust = 0.9,
          size = 2
        ) +
        ylab(ylab)+
        scale_color_continuous(type = "viridis")+ 
        facet_wrap(~sample, ncol = 2)+ 
        geom_vline(data = medians,aes(xintercept=log10_medians)) +
        labs(color = "Transcripts per bin") +
        guides(color = guide_colourbar(barheight=0.5)) +
        theme_light()  +
        theme(legend.position = "bottom") + 
        theme(axis.text.x = element_text(size=8))
    } else {
      p <- ggplot(all_data, aes(x= log10(Expression), y = Length, fill = n)) +
        geom_density_ridges()+ 
        facet_wrap(~sample, ncol = 2) +
        geom_vline(data = medians,aes(xintercept=log10_medians)) +
        stat_summary(
          fun.data = stat_box_data, 
          geom = "text", 
          hjust = 1,
          vjust = -0.75,
          size = 2
        ) +
        ylab(ylab)+
        scale_fill_continuous(type = "viridis")+ 
        labs(fill = "Transcripts\nper bin") +
        guides(fill = guide_colourbar(barheight=0.5)) +
        theme_light()  +
        theme(legend.position = "bottom") + 
        theme(axis.text.x = element_text(size=8))
    }
    p <- p + geom_text(data = data.frame(label=v_labels, sample= samples),
                                aes(x = -Inf, y = Inf, label=v_labels, group = sample),
                                inherit.aes = FALSE,
                                hjust = -1,
                                vjust = 1.2,
                                size = 2)
    return(p)
  }
}
