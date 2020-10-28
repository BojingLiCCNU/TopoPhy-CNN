myget_comm_params<-function (comm, mar = 2, distr, ...) 
{
  apply(comm, mar, function(x) {
    x <- as.numeric(x)
    ll <- c(list(), myfitdistr(x, distr, ...)$par)
    ll$mean <- mean(x)
    ll
  })
}