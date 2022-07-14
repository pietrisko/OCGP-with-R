GPR_OCC_score <-  function (training,test){
  
  x <- as.matrix(training)
  y <- as.matrix(test)
  
  if (sum(grepl('^ENS',rownames(x)))==nrow(x)) {
    rownames(x) <- sapply(strsplit(rownames(x),'\\.'), function(x) x[1])
    rownames(y) <- sapply(strsplit(rownames(y),'\\.'), function(x) x[1])
  }
  
  cg <- intersect(rownames(x),rownames(y))
  
  if (length(cg)<1000) {
    stop('Common genes (rownames) are less than 1000',call. = T)
  }
  
  message(paste0('Number of common genes (rownames) between training and test sets: ',length(cg)))
  
  x <- x[cg,]
  y <- y[cg,]
  
  cor_xx <- cor(x,x,method = 'spearman',use = "complete.obs")
  cor_xy <- cor(x,y,method = 'spearman',use = "complete.obs")
  
  noise=0.01
  cor_xx=cor_xx+noise*diag(1,nrow(cor_xx),ncol(cor_xx))
  
  L = t(chol(cor_xx))
  alpha = solve(t(L),solve(L,matrix(1,nrow(cor_xx),1)))
  
  score = t(cor_xy)%*%alpha
  return(score)
}
