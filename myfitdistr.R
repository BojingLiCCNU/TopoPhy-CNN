myfitdistr<-function (x, densfun, start, control, ...) 
{
  if (class(x) != "numeric") 
    stop("Error: input must be numeric vector")
  Call <- match.call(expand.dots = TRUE)
  if (missing(start)) 
    start <- NULL
  if (missing(control)) 
    control <- list(fnscale = 1e+12, factr = 0.01, maxit = 10)
  dots <- names(list(...))
  dots <- dots[!is.element(dots, c("upper", "lower"))]
  if (missing(x) || length(x) == 0L || mode(x) != "numeric") 
    stop("'x' must be a non-empty numeric vector")
  if (any(!is.finite(x))) 
    stop("'x' contains missing or infinite values")
  if (missing(densfun) || !(is.function(densfun) || is.character(densfun))) 
    stop("'densfun' must be supplied as a function or name")
  n <- length(x)
  if (is.character(densfun)) {
    distname <- tolower(densfun)
    densfun <- switch(distname, zipois = VGAM::dzipois, 
                      zinegbin = VGAM::dzinegbin, negbin = dnbinom, pois = dpois, 
                      lognorm = "dlnorm", NULL)
    if (is.null(densfun)) 
      stop("unsupported distribution")
  }
  if (distname == "lognorm") {
    meanlog <- mean(log(x))
    sdlog <- sd(log(x))
    return(list(par = list(meanlog = meanlog, sdlog = sdlog)))
  }
  if (distname == "zipois") {
    if (!is.null(start)) 
      stop(gettextf("supplying pars for the %s distribution is not supported", 
                    "Poisson"), domain = NA)
    whichz <- which(x == 0)
    which1 <- which(x == 1)
    max <- abs(length(whichz) - length(which1))
    max <- max - max * 0.1
    zind <- na.omit(whichz[1:max])
    tempx <- x[-zind]
    pstr0 <- length(which(x == 0))/length(x)
    pstr0 <- abs(pstr0 - exp(-mean(tempx)))
    estimate <- mean(x)/(1 - pstr0)
    vars <- ((1 - pstr0) * (estimate^2 + estimate)) - ((1 - 
                                                          pstr0) * estimate)^2
    sds <- sqrt(vars)
    start <- c(lambda = estimate, pstr0 = pstr0)
    start <- start[!is.element(names(start), dots)]
    lower <- c(1e-04, 0)
    upper <- c(Inf, 0.99)
    loglikfn <- match.fun(logLikzip)
  }
  if (distname == "gamma" && is.null(start)) {
    if (any(x < 0)) 
      stop("gamma values must be >= 0")
    m <- mean(x)
    v <- var(x)
    start <- list(shape = m^2/v, rate = m/v)
    start <- start[!is.element(names(start), dots)]
    control <- c(control, list(parscale = c(1, start$rate)))
  }
  else if (distname == "zinegbin" && is.null(start)) {
    whichz <- which(x == 0)
    which1 <- which(x == 1)
    max <- abs(length(whichz) - length(which1))
    max <- max - max * 0.1
    zind <- na.omit(whichz[1:max])
    tempx <- x[-zind]
    pstr0 <- length(which(x == 0))/length(x)
    if (pstr0 != 0) 
      pstr0 <- abs(pstr0 - exp(-mean(tempx)))
    m <- mean(x)
    v <- var(x)
    if (pstr0 == 0) {
      if (v > n) 
        size <- m^2/(v - m)
      else size <- 100
      estimate <- m
    }
    else {
      v <- var(tempx)
      m <- mean(tempx)
      estimate <- m/(1 - pstr0)
      if (v < m) 
        size <- 100
      else size <- (m^2/(v - m)) * ((1 - pstr0) * estimate)
    }
    lower <- c(0.01, 0.01, 0)
    upper <- c(10000, 10000, 0.99)
    start <- c(size = size, munb = estimate, pstr0 = pstr0)
    start <- start[!is.element(names(start), dots)]
    loglikfn <- match.fun(logLikzinb)
  }
  else if (distname == "negbin") {
    m <- mean(x)
    v <- var(x)
    size <- if (v > m) 
      m^2/((v/2) - m)
    else 100
    start <- list(size = size, mu = m)
    start <- start[!is.element(names(start), dots)]
    lower <- c(0.01, 0.01)
    upper <- c(10000, 10000)
    loglikfn <- match.fun(logLiknb)
  }
  else if (distname == "pois") {
    m <- mean(x)
    return(list(par = list(lambda = m)))
  }
  start <- pmax(start, lower)
  start <- pmin(start, upper)
  names(upper) <- names(lower) <- names(start)
  res <- optim(start, loglikfn, x = x, method = "L-BFGS-B", 
               lower = lower, upper = upper, control = control, ...)
  return(res)
}