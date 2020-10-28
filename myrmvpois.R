myrmvpois<-function (n, mu, Sigma, ...) 
{
  Cor <- cov2cor(Sigma)
  
  for( i in 1:nrow(Cor)){
    for(j in 1:length(Cor[i,]))
      if (is.nan(Cor[i,j])){
        Cor[i,j]<-0
      }
  }
  # print(Cor)
  
  
  
  SDs <- sqrt(diag(Sigma))
  d <- length(SDs)
  if (length(mu) != length(SDs)) 
    stop("Sigma and mu/lambdas dimensions don't match")
  if (length(mu) == 1) 
    stop("Need more than 1 variable")
  normd <- myrmvnorm(n, rep(0, d), Cor)
  unif <- pnorm(normd)
  data <- t(matrix(qpois(t(unif), mu, ...), d, n))
  data <- myfixInf(data)
  return(data)
}
