  n <- 10000
  v <- rchisq(n, df=4)
  w <- v + rnorm(n)
  x <- rchisq(n, df=2) + w
  y <- x + rnorm(n)
  z <- rnorm(n, sd=1)
  dat <- cbind(v,w,x,y,z)
  
  cont.pairs(dat, prob.levels=c(1,2), cex=0.5, npix=100, sigma=TRUE,
              dot.level=2, prob1d=c(1, 2), smooth1d=TRUE,
              cex.lab=1.3)
