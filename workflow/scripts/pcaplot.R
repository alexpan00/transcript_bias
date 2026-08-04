pcaplot <- function(mydata){
  counts <- exprs(mydata)
  pData_df <- pData(mydata)
  col1 <- if (ncol(pData_df) >= 2) colnames(pData_df)[2] else colnames(pData_df)[1]
  col2 <- if (ncol(pData_df) >= 3) colnames(pData_df)[3] else col1

  if (ncol(counts) < 2) {
    scores <- data.frame(PC1 = 0, PC2 = 0)
    scores[[col1]] <- pData_df[[col1]]
    scores[[col2]] <- pData_df[[col2]]
    if (col1 == col2) {
      p <- ggplot(scores, aes(PC1, PC2, color = .data[[col1]])) +
        geom_point() + 
        xlab("PC1: 0%") +
        ylab("PC2: 0%") +
        theme_light()
    } else {
      p <- ggplot(scores, aes(PC1, PC2, color = .data[[col1]], shape = .data[[col2]])) +
        geom_point() + 
        xlab("PC1: 0%") +
        ylab("PC2: 0%") +
        theme_light()
    }
    return(p)
  }

  sds <- apply(counts, 1, sd)
  sds[is.na(sds)] <- 0
  if (all(sds == 0)){
    counts <- matrix(rnorm(nrow(counts)*ncol(counts)),
                     nrow = nrow(counts),
                     ncol = ncol(counts))
    sds <- apply(counts, 1, sd)
    sds[is.na(sds)] <- 0
  }
  # remove rows with zero variance
  counts <- counts[sds != 0, , drop=FALSE]
  pca <- prcomp(t(counts), scale. = T, center = T)
  eigs <- pca$sdev^2
  pv <- round(eigs/sum(eigs), 4)*100
  scores <- as.data.frame(pca$x[,1:min(2, ncol(pca$x)), drop=FALSE])
  if (!"PC2" %in% colnames(scores)) scores$PC2 <- 0
  if (length(pv) < 2 || is.na(pv[2])) pv[2] <- 0
  
  scores[[col1]] <- pData_df[[col1]]
  scores[[col2]] <- pData_df[[col2]]
  
  if (col1 == col2) {
    p <- ggplot(scores, aes(PC1, PC2, color = .data[[col1]])) +
      geom_point() + 
      xlab(paste0("PC1: ",pv[1], "%")) +
      ylab(paste0("PC2: ",pv[2], "%")) +
      theme_light()
  } else {
    p <- ggplot(scores, aes(PC1, PC2, color = .data[[col1]], shape = .data[[col2]])) +
      geom_point() + 
      xlab(paste0("PC1: ",pv[1], "%")) +
      ylab(paste0("PC2: ",pv[2], "%")) +
      theme_light()
  }
  return(p)
}