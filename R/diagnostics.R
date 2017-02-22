# -----------------------------------------
# diagnostics.R
#
# functions to help diagnose the output of MCMC routines.
# -----------------------------------------

mcmc.diag.plot <- function(chain, max.chain = 20) {
  
  # mcmc.diag.plot - Generate output plots to help assess the convergence 
  #                  of MCMC output.
  # Inputs: 
  #   chain  - (list) containing output data from an MCMC, including
  #     theta  - (array) n * ndim array of posterior samples
  #             n samples of ndim vectors of parameters
  #     method - (string) name of MCMC method used
  #     nwalkers - number of walkers used (if method = 'gw.mcmc')
  #     nchains  - number of walkers used (if method = 'mh.mcmc')
  #   max.chains - maximum no. chains to overlay on trace plot
  #   
  # Value
  #   none
  #
  # Generate plots of chain traces (variable values vs. iteration)
  # auto-correlation plots and 1D densities (histograms).
  # History:
  #  13/07/16 - First working version
  #
  # Simon Vaughan, University of Leicester
  # Copyright (C) 2016 Simon Vaughan
  
  # check the input arguments
  if (missing(chain)) stop('Must specify chain input.')
  if (!"theta" %in% names(chain)) stop('** chain input list is missing theta.')
  
  m <- rbind(c(1, 1), c(2, 3))
  layout(m)
  par(mar = c(5, 5, 2, 1))
  
  # total length of chains
  n <- length(chain$theta[,1])
  
  # dimensions of density / no. variables
  ndim <- NCOL(chain$theta)
  
  nchains <- 1
  # no. walkers if using ensemble method
  if (chain$method == 'gw.mcmc') {
    nchains <- chain$nwalkers
  } 
  if (chain$method == 'mh.mcmc') {
    nchains <- chain$nchains
  } 
  
  # total no. cycles
  ncycles <- floor(n / nchains)
  
  # do the columns of theta come with names?
  if (is.null(colnames(chain$theta))) {
    names <- paste("parameter", 1:ndim)
  } else {
    names <- colnames(chain$theta)
  }
  
  # prepare an array for output    
  zeroACF <- array(NA, dim = ndim)
  
  # -----------------------------------------
  # trace plots of walkers
  
  if (nchains > max.chain) {
    walker.sample <- sample(1:nchains, max.chain)
  } else {
    walker.sample <- 1:nchains
  }
  
  t <- 1:ncycles
  
  for (i in 1:ndim) {
    
    # -----------------------------------------
    # Trace plot - for each walker plot theta[j] vs. iteration
    
    x <- chain$theta[,i]
    dim(x) <- c(ncycles, nchains)
    x.mean <- rowMeans(x)
    
    plot(x[,1], bty = "n", main = paste(names[i], "- trace"), type = "n",
         xlim = c(0, ncycles), ylim=range(x), 
         xlab = "iteration", ylab = names[i])
    
    n.col <- length(walker.sample)
    k <- 0
    for (j in walker.sample) {
      k <- k + 1
      lines(t, x[,j], col=rainbow(n.col)[k], type = "l")
    }  
    lines(t, x.mean, lwd = 3)
    
    # -----------------------------------------
    # Auto-correlation functions (ACFs) 
    # Compute the mean of the ACFs of each walkers, and
    # the ACF of the mean of all the walkers.
    
    lag.max <- ncycles
    
    for (j in 1:nchains) {
      acf.i <- acf(x[,j], lag.max = lag.max, plot = FALSE)
      if (j == 1) acf.j <- rep(0, length(acf.i$lag))
      acf.j <- acf.j + acf.i$acf
    }
    acf.j <- acf.j / nchains
    lag.j <- acf.i$lag
    acf.mean <- acf(x.mean, lag.max = lag.max, plot = FALSE)
    
    plot(lag.j, acf.j, type = "l", bty = "n", ylim=c(-1,1), main = "ACF", 
         xlab = "lag", ylab = "ACF")
    abline(h = 0, lty=2)
    lines(lag.j, acf.mean$acf, lwd=3)
    
    sign.change <- diff(sign(acf.mean$acf))
    if (!all(sign.change == 0)) zeroACF[i] <- min(which(sign.change != 0))
    
    # -----------------------------------------
    # 1D histogram
    
    if (exists('plot.dist')) {
      plot.dist(x, xlim = range(x), fill.col = "steelblue3", xlab = names[i],
                main = "density", breaks = 60)
      axis(1)
    } else {
      hist(x, breaks = 60, main = "density", border = NA, 
           col = "steelblue3", xlab = names[i])
    }
    
  }
}

# ------------------------------------------------
# for each parameter theta[1]...theta[M]
# calculate the R.hat statistic (Gelman & Rubin 1992)
# See also Gelman et al. (2004, sect 11.6)

Rhat <- function(theta) {
  
    L <- dim(theta)[1]
    sj <- apply(theta, 2, var)
    W <- mean(sj)
    mj <- apply(theta, 2, mean)
    B <- var(mj)*L
    v <- (L-1)/L*W + B/L
    R.hat <- sqrt(v/W)
    
  return(R.hat)
}

# ------------------------------------------------
# Perform checks for convergence of multiple
# Markov chains. Uses Gelman & Rubin's R.hat
# for each parameter, and also a visual check of
# the 80% regions for each parameter.

mcmc.conv <- function(chain) {

  # check the input arguments
  if (missing(chain)) stop('Must specify chain input.')
  if (!"theta" %in% names(chain)) stop('** chain input list is missing theta.')
  
  # total length of chains
  n <- length(chain$theta[,1])
  
  # dimensions of density / no. variables
  ndim <- NCOL(chain$theta)
  
  nchains <- 1
  # no. walkers if using ensemble method
  if (chain$method == 'gw.mcmc') {
    nchains <- chain$nwalkers
  } 
  if (chain$method == 'mh.mcmc') {
    nchains <- chain$nchains
  } 
  
  # total no. cycles
  ncycles <- floor(n / nchains)
  
  # do the columns of theta come with names?
  if (is.null(colnames(chain$theta))) {
    names <- paste("parameter", 1:ndim)
  } else {
    names <- colnames(chain$theta)
  }
  
  R.hat <- array(NA, dim = ndim)
  offset <- array(NA, dim = ndim)
  ci.lo <- array(NA, dim = c(ndim, nchains))
  ci.hi <- array(NA, dim = c(ndim, nchains))
  
  # loop over each parameter
  for (i in 1:ndim) {
  
    # Calculate R.hat (Gelman & Rubin 1992) for each
    # parameter as a check for convergence of chains
    theta <- chain$theta[, i]
    dim(theta) <- c(ncycles, nchains)
    R.hat[i] <- Rhat(theta)
  
  # calculate and plot the 80% intervals from each chain
    intv <- apply(theta, 2, quantile, prob = c(0.1, 0.9))
    ci.lo[i, ] <- intv[1, ]
    ci.hi[i, ] <- intv[2, ]
    offset[i] <- mean(theta)
    ci.lo[i, ] <- ci.lo[i, ] - offset[i] 
    ci.hi[i, ] <- ci.hi[i, ] - offset[i] 
    scale <- mean(ci.hi[i, ] - ci.lo[i, ]) / 2
    ci.lo[i, ] <- ci.lo[i, ] / scale
    ci.hi[i, ] <- ci.hi[i, ] / scale
  }

  layout(t(c(1, 2)), widths = c(0.3, 0.7))
  
  plot(R.hat, 1:ndim, xlim = c(0.7, 2), ylim = c(0.5, ndim+0.5), 
       bty = "n", pch = 16, yaxp = c(1, ndim, ndim-1),  
       xlab = "R.hat", ylab = "Parameter")
  
  # plot the R.hat results
  
  abline(v = 1.1, lty = 2)
  axis(3)
  
  # plot Credible Intervals (80%)
  
  plot(rep(0, ndim+1), 1:(ndim+1), ylim = c(0.5, ndim+0.5), xlim = c(-3, 3), 
       type = "n", bty = "n", ylab = "Parameter", 
       xlab = "80% region (scaled)", yaxp = c(1, ndim, ndim-1))
  axis(3)
  for (i in 1:ndim) {
    for (j in 1:nchains) {
      x <- i + (j-1) / nchains / 4
      segments(ci.lo[i, j], x, ci.hi[i, j], x, col=j)
    }
  }
}


