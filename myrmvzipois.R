myrmvzipois<-function (n, mu, Sigma = diag(length(mu)), lambdas, ps, ...) 
{
  d <- ncol(Sigma)
  Cor <- cov2cor(Sigma)
  
  for( i in 1:nrow(Cor)){
    for(j in 1:length(Cor[i,]))
      if (is.nan(Cor[i,j])){
        Cor[i,j]<-0
      }
  }
  # print(Cor)
  
  
  
  
  
  SDs <- sqrt(diag(Sigma))
  if (missing(lambdas) || missing(ps)) {
    if (missing(mu)) 
      stop("Need to supply mu")
    if (length(mu) != length(SDs)) 
      stop("Sigma and mu dimensions don't match")
    lambdas <- unlist(lapply(1:length(SDs), function(i) .zipois_getLam(mu[i], 
                                                                       SDs[i])))
    ps <- unlist(lapply(1:length(SDs), function(i) .zipois_getP(mu[i], 
                                                                SDs[i])))
  }
  if (length(lambdas) != length(SDs)) 
    stop("Sigma and mu/lambdas dimensions don't match")
  if (length(lambdas) == 1) 
    stop("Need more than 1 variable")
  normd <- myrmvnorm(n, rep(0, d), Cor)
  unif <- pnorm(normd)
  data <- t(matrix(VGAM::qzipois(t(unif), lambdas, pstr0 = ps, 
                                 ...), d, n))
  data <- myfixInf(data)
  return(data)
}
