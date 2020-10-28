myfixInf<-function (data) 
{
  if (any(is.infinite(data))) {
    data <- apply(data, 2, function(x) {
      if (any(is.infinite(x))) {
        x[ind <- which(is.infinite(x))] <- NA
        x[ind] <- max(x, na.rm = TRUE) + 1
      }
      x
    })
  }
  data
}
