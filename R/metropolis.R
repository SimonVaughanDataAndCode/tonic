# ------------------------------------------------
# metropolis.R
# R functions for using the random walk Metropolis-Hastings
# algorithm for MCMC
#
# History:
#  14/07/16 - v0.1 - First working version
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
# 
# See e.g. Roberts (2015) http://arxiv.org/abs/1504.01896 
# ------------------------------------------------


#' Random Walk Metropolis-Hastings sampler using multiple chains.
#' 
#' \code{mh_sampler} returns posterior samples from a RW-MH sampler.
#' 
#' A simple, pure R function to carry generate MCMC output using the random walk
#' Metropolis-Hastings algorithm. Uses a Normal or t distribution for the
#' proposal.
#' 
#' @param nchains (integer) number of distinct 'chains' to run.
#' @param adapt (logical) if \code{TRUE} then adjust the covariance matrix using the 
#'              samples from the burn in period.
#' @param cov (array/matrix) covariance matrix for proposal distribution.
#' @param proposal (string) Use "normal" proposal, or else "t" distribution.
#' @param merge.chains (logical) combine output from all chains into one
#' @inheritParams gw_sampler
#' 
#' @return
#'  A list with components
#'   \item{theta}{(array) \code{nsteps} samples from M-dimensional posterior [nsteps rows, M columns]}
#'   \item{func}{(string) name of posterior function sampled}
#'   \item{lpost}{(vector) nsteps values of the LogPosterior density at each sample position}
#'   \item{method}{(string) sampling method used ("mh_sampler")}
#'   \item{nchains}{(integer) number of chains used}
#'   \item{accept.rate}{(float) the fraction of proposals accepted.}
#' If \code{merge.chains = FALSE} then \code{theta} will be a 3D array with
#' dimensions \code{nchains * (nsteps/nchains) * M}.
#' 
#' @section Notes:
#' Generate samples from some 'target' density function specified by the 
#' \code{posterior} function. Uses the random walk Metropolis-Hastings
#' algorithm. The 'proposal' distribution is either a multivariate normal
#' (default) or a Student's t distribution (with \code{df=3}). The latter has
#' slightly fatter tails, so might be preferred for non-normal target densities.
#' 
#' This function generates a total of \code{nsteps} samples after a 'burn in'
#' period of \code{burn.in} steps. Allowing for a burn in period helps remove
#' memory of the starting position (important if this is not well chosen). We
#' run \code{nchains} in parallel, from slightly different starting positions.
#' Each chain only needs to run for \code{nsteps / nchain} iterations (after the
#' burn in is complete). Having multiple chains helps assess convergence: are
#' they 'well mixed' meanng that samples from each walker appear to be drawn
#' from the same distribution? We can examine the moments and histograms of the 
#' chains to check this.
#' 
#' The proposal density (\code{normal} or \code{t}) requires a covariance
#' matrix. This is specified on input by the \code{cov} argument. If not
#' specified this will default to an identity matrix. Chosing a sensible
#' covariance matrix for the proposal is part of the 'art' of MCMC. If the
#' proposal is badly chosen - and nothing like the target density - then
#' convergence (in distribution) will be very slow.
#' 
#' MCMC practitioners often aim for an acceptance rate of 0.1-0.5. If the 
#' acceptance rate is very low it means there will be many repeats of the same 
#' values before new ones are accepted, and so the random walk will explore the 
#' target density slowly. This is usually because the proposal density is too 
#' large compared to the target density - try reducing the scale of the 
#' covariance matrix. If the acceptance rate is very high (approaching 100\%) 
#' then probably the proposal distribution is too small, each new proposed 
#' position is very close to the current position (and so of similar density and
#' highly likely to be accepted). Again, this means the random walk will move
#' only slowly through the target density. In this case, making the covariance
#' larger usually helps.
#' 
#' The target density function: The target density should is specified by the 
#' \code{posterior} function. In fact, this function should compute the log
#' density given parameters theta and any other arguments, i.e. log(p) = 
#' \code{posterior(theta, ...)} where the vector of parameters \code{theta} is
#' the first argument of the posterior function. Where the density is zero or
#' not defined, e.g. because prior = 0 for certain values of the parameters, it 
#' should return \code{-Inf}. Otherwise, the output of \code{posterior(theta,
#' ...)} should be a real, scalar value.
#'  
#' @seealso \code{\link{chain_convergence}}, \code{\link{gw_sampler}},
#' \code{\link{chain_diagnosis}}, \code{\link{contour_matrix}}
#' 
#' @examples
#' my_posterior <- function(theta) {
#'   cov <- matrix(c(1,0.98,0.8,0.98,1.0,0.97,0.8,0.97,2.0), nrow = 3)
#'   logP <- mvtnorm::dmvnorm(theta, mean = c(-1, 2, 0), sigma = cov, log = TRUE)
#'   return(logP)
#' }
#' chain <- mh_sampler(my_posterior, theta.0 = c(0,0,0), nsteps = 10e4, burn.in = 1e4) 
#'
#' @export
mh_sampler <- function(posterior, 
                   theta.0, 
                   nsteps = 1E4,
                   nchains = 5,
                   burn.in = 2000,
                   update = 5,
                   chatter = 0,
                   cov = NULL,
                   thin = NULL, 
                   merge.chains = TRUE,
                   adapt = FALSE,
                   proposal = "normal", ...) {

   
  # check the input arguments
  if (missing(theta.0)) stop('Must specify theta.0 start position.')
  if (missing(posterior)) stop('Must specify name of posterior function')
  if (!exists('posterior')) {
    stop('The specified log density function does not exist.')
  }

  # dimensions of the PDF
  M <- length(theta.0)
  
  # default covariance matrix if none supplied
  cov <- diag(M)

  # ensure the number of steps, walkers, etc. are integers  
  nsteps <- as.integer(nsteps)
  nchains <- as.integer(nchains)
  burn.in <- as.integer(burn.in)
  
  # prepare a working array
  p.prop <- array(NA, dim = nchains)
  
  # number of interations needed
  nrows.keep <- ceiling(nsteps / nchains)
  if (nrows.keep < 10) stop('Make nsteps larger')
  nrows.burnin <- ceiling(burn.in / nchains)
  ncycles <- nrows.keep + nrows.burnin
  
  # initialise the array for output. There are two additional 'columns': The M+1
  # column stores the accept/reject flag (0=rejected, 1=accepted proposal at 
  # each update). This is useful for tracking the acceptance rate. The M+2 
  # column stores the log posterior (PDF) values at the current position of the 
  # chain. This saves recomputing the posterior density at the current position 
  # when evaluating the accept probability.
  theta <- array(NA, dim=c(ncycles, nchains, M+2))

  # starting locations of chains
  for (j in 1:nchains) {
     z <- mvtnorm::rmvnorm(1, mean = theta.0, sigma = cov/4)
     theta[1, j, 1:M] <- z  # randomized start position
     theta[1, j, M+1] <- 1  # accept
     theta[1, j, M+2] <- posterior(z, ...)
  }

  # check that the density is positive at the starting positions
  start.densities <- array(NA, dim = M)
  for (j in 1:nchains) {
    start.densities[j] <- posterior(theta[1, j, 1:M], ...)
  }
  if (!all(is.finite(start.densities))) {
    stop('** Non-finite values of log posterior at initial values.')
  }
  
  # loop over the chain
  
  i.count <- 1
  for (i in 2:ncycles) {
    
    # carry forward the previous theta value
    # will over-write if update occurs
    theta[i, , ] <- theta[i-1, , ]
    theta[i, , M+1] <- 0
    
    # draw a value from the proposal distribution
    # either Normal or Student's t distribution (df=3)
    if (proposal == "normal") {
      z <- mvtnorm::rmvnorm(nchains, sigma = cov)
    } else {
      z <- mvtnorm::rmvt(nchains, sigma = cov/3, df = 3)
    }
    # add these random vectors to the previous position
    # to produce a new 'proposal' position
    theta.prop <- theta[i, 1:nchains, 1:M] + z
    
    # now for each chain compute the log posterior density
    # at the proposed position.
    for (j in 1:nchains) {
      p.prop[j] <- posterior(theta.prop[j, ], ...)
    }
    
    # compute ratio of posteriors at new and old locations
    #   r = p(theta.prop) / p(theta[t-1])
    # in terms of the log posterior function 
    # this is log[p(theta.old)] - log[p(theta.new)]
    # (log.r is a vector with nchains elements.)
    p.pres <- theta[i, , M+2] 
    log.r <- p.prop - p.pres
    
    # decide whether or not to update theta.
    # theta is updated with probability min(r,1)
    # otherwise it is left as before.
    r <- exp( pmin(log.r, 0) )
    u <- runif(nchains)
    mask <- (r >= u)
    theta[i, mask, 1:M] <- theta.prop[mask, 1:M]
    theta[i, mask, M+1] <- 1
    theta[i, mask, M+2] <- p.prop[mask]

    # progress report to user if requested
    if (chatter > 0) {
      if (i %% update == 0) {
        accept.rate <- mean( theta[i.count:i, , M+1], na.rm=TRUE ) 
        if (i <= nrows.burnin) {
          cat("\r-- Burn in cycle", i, "of", nrows.burnin)
        } else {
          cat("\r-- Cycle", i-nrows.burnin, "of", ncycles-nrows.burnin)
        }
        cat(". Acceptance rate:", signif(accept.rate*100, 2), "%")
      }
      if (i == ncycles) cat('', fill=TRUE)
      if (i >= nrows.burnin & i.count == 1) {
        cat(' - Finished burn-in', fill=TRUE)
        i.count <- nrows.burnin + 1
      }
    }
    
    # if adapt == TRUE then adjust the covariance matrix
    # using the samples from the burn in period.
    if (i == nrows.burnin) {
      if (adapt == TRUE) {
        theta.burn <- matrix(theta[1:nrows.burnin, , 1:M], 
                             nrow = nchains * nrows.burnin, 
                             ncol = M, byrow = FALSE)
        cov <- cov(theta.burn)
      }
    }
    
  } # end of loop over i = 2, ncycles
  
  # strip off the burn-in period and keep only nsteps 
  nrows <- nrows.keep
  theta <- theta[(1:nrows) + nrows.burnin, , ]

  # thin the output by keeping only every few rows
  if (!is.null(thin)) {
    nrow.keep <- floor(nrows / thin)
    mask <- (1:nrow.keep) * thin
    theta <- theta[mask, , ]
  }
  
  # Strip off the acceptance and log(posterior) columns 
  accept <- theta[, , M+1]
  lpost <- as.vector(theta[, , M+2])
  theta <- theta[, , 1:M]
  
  # check the acceptance rate (column M+1).
  # Also strip off the log(posterior) values which are no longer needed.
  accept.rate <- mean(accept, na.rm = TRUE)
  if (chatter > 0) {
    cat('\n-- Final acceptance rate: ', accept.rate, fill = TRUE)
    if (accept.rate < 0.05) {
      cat('-- Low acceptance rate. Consider the following suggestions:', 
          fill = TRUE)
      cat('-- 1. Adjust the start position: theta.0,', fill = TRUE)
      cat('-- 2. Decrease the covariance for proposal distribution.', fill = TRUE
      )
    }
  }
  
  # reshape the array from [nrows, M, nwalkers] to [nrows*nwalkers, M]
  # so each column is one variable, each row is one sample from the M
  # M-dimensional distribution.
  if (merge.chains == TRUE) {
    nrows <- dim(theta)[1]
    theta <- matrix(theta, nrow = nchains * nrows, ncol = M, byrow = FALSE)
  }
  
  # return the final array
  return(list(theta = theta,
              func = deparse(substitute(posterior)),
              lpost = lpost,
              method = "mh_sampler",
              nchains = nchains,
              accept.rate = accept.rate))
}


