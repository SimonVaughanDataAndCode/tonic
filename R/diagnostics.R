# -----------------------------------------
# diagnostics.R
#
# functions to help diagnose the output of MCMC routines.
#
# History:
#  13/07/16 - First working version
#  25/02/17 - save/restore graphical parms, added documentation
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
# -----------------------------------------

#' Plots for assessing MCMC output.
#' 
#' \code{chain_diagnosis} Generate plots to help assess the convergence 
#'                  of MCMC output.
#'                  
#' Generate plots of chain traces (variable values vs. iteration), 
#' auto-correlation plots and 1D densities (histograms).
#'
#' @param chain (list) containing output data from an MCMC, including
#' @param max.chains maximum no. chains to overlay on trace plot
#'   
#' @return none
#' 
#' @section Further details:
#' The input \code{chain} should be a list such as produced by 
#'   \code{gw.mcmc} or \code{mh_sampler} that contains the following:
#' \describe{ 
#' \item{theta}{(array) n * ndim array of posterior samples
#'             n samples of ndim vectors of parameters}
#' \item{method}{(string) name of MCMC method used}
#' \item{nwalkers}{number of walkers used (if method = 'gw.mcmc')}
#' \item{nchains}{number of walkers used (if method = 'mh_sampler')}
#' }
#' 
#' @seealso \code{\link{gw_sampler}}, \code{\link{mh_sampler}}
#' 
#' @examples 
#' my_posterior <- function(theta) {
#'   cov <- matrix(c(1,0.98,0.8,0.98,1.0,0.97,0.8,0.97,2.0), nrow = 3)
#'   logP <- mvtnorm::dmvnorm(theta, mean = c(-1, 2, 0), sigma = cov, log = TRUE)
#'   return(logP)
#' }
#' chain <- gw.mcmc(my_posterior, theta.0 = c(0,0,0), nsteps = 10e4, burn.in = 1e4) 
#' chain_diagnosis(chain)
#'
#' @export
chain_diagnosis <- function(chain, max.chain = 20) {
  
  # check the input arguments
  if (missing(chain)) stop('Must specify chain input.')
  if (!"theta" %in% names(chain)) stop('** chain input list is missing theta.')

  # store the current graphical parameters  
  retire <- par(no.readonly = TRUE)
  
  # now change the graphical parameters to make room for 3 panels
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
  if (chain$method == 'mh_sampler') {
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
  
  # restore graphical parameters
  par(retire)
}

# ------------------------------------------------

#' compute the Gelman & Rubin R.hat statistic.
#' 
#' \code{Rhat} computes the Gelman & Rubin R.hat statistic.
#' Given an array of \code{theta} values, for some parameter, produced by
#' several chains, compute the R.hat statistic as a check for convergence. The
#' \code{R.hat} statistic (Gelman & Rubin 1992) should be close to zero if the
#' chains are converging. See also Gelman et al. (2004, sect 11.6)
#' 
#' @param theta (array) of samples of one parameters. \code{ncycles} rows by
#' \code{nchains} columns. We compare the within-chain and between-chain
#' variances.
#' 
#' @return
#'  The R.hat statistic (scalar).
#'  
#' @seealso \code{\link{chain_convergence}}
#'
#' @export
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

#' Perform checks for convergence of multiple Markov chains. 
#' 
#' \code{chain_convergence} checks convergence of multiple Markov chains.
#' 
#' Uses Gelman & Rubin's R.hat for each parameter, and also a visual check of 
#' the 80\% regions for each parameter.
#' 
#' @param chain (array) of MCMC samples, \code{N} rows (sampled) by
#' \code{M} columns (variables). 
#' 
#' @return
#'  The values of R.hat for each of the \code{M} variables.
#'  
#' @seealso \code{\link{Rhat}}
#'
#' @export
chain_convergence <- function(chain) {

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
  if (chain$method == 'mh_sampler') {
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
  return(R.hat)
}


