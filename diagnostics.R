m <- rbind(c(1, 1), c(2, 3))
layout(m)
par(mar = c(3, 3, 2, 1))

n <- length(theta[,1])

nwalkers <- 100
ncycles <- floor(n / nwalkers)
ndim <- NCOL(theta)

if (is.null(colnames(theta))) {
  names <- paste("parameter", 1:ndim)
} else {
  names <- colnames(theta)
}

if (nwalkers > 10) {
  walker.sample <- sample(1:nwalkers, 10)
} else {
  walker.sample <- 1:nwalkers
}

t <- 1:ncycles

 for (i in 1:ndim) {
  
  x <- theta[,i]
  dim(x) <- c(ncycles, nwalkers)
  x.mean <- rowMeans(x)
  
  plot(x[,1], bty = "n", main = names[i], type = "n",
       xlim = c(0, ncycles), ylim=range(x))

  k <- 0
  for (j in walker.sample) {
    k <- k + 1
    lines(t, x[,j], col=rainbow(10)[k], type = "s")
  }  
  lines(t, x.mean, lwd = 3)
  
  
  lag.max <- ncycles
  
  for (j in 1:nwalkers) {
    acf.i <- acf(x[,j], lag.max = lag.max, plot = FALSE)
    if (j == 1) acf.j <- rep(0, length(acf.i$lag))
    acf.j <- acf.j + acf.i$acf
  }
  acf.j <- acf.j / nwalkers
  lag.j <- acf.i$lag
  acf.mean <- acf(x.mean, lag.max = lag.max, plot = FALSE)
  
  plot(lag.j, acf.j, type = "l", bty = "n", ylim=c(-1,1))
  abline(h = 0, lty=2)
  lines(lag.j, acf.mean$acf, lwd=3)
  
  breaks = min(100, n/25)
  h <- hist(x, breaks = breaks, plot = FALSE)
  
  plot(h$mids, h$density, type = "s", bty = "n")
  
}
  