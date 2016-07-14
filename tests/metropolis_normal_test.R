# ------------------------------------------------

# a 3D un-normalised distribution comprising 3 separate Normal densities
# separated in all 3 dimensions

my.pdf <- function(theta) {
  n <- 3
  cov <- diag(n) * 0.1
  P1 <- mvtnorm::dmvnorm(theta, mean=c(0,1,2), sigma=cov)
  cov <- diag(n) * 0.03
  P2 <- mvtnorm::dmvnorm(theta, mean=c(-1,0,1), sigma=cov)
  cov <- diag(n) * 0.03
  P3 <- mvtnorm::dmvnorm(theta, mean=c(1,-1,-1), sigma=cov)
  logP <- log(0.9*P1 + 0.0*P2 + 0.9*P3)
  return(logP)
}

# ------------------------------------------------
# test MCMC function mh.mcmc

chain <- mh.mcmc(my.pdf, theta.0=c(0,0,0), nsteps=2e4, cov = 15*diag(3),
                  chatter=1, burn.in = 1e4, adapt = TRUE)

source("diagnostics.R")
mcmc.conv(chain)

# examine diagnostic plots
source('~/R/tonic/plot_contour.R')
mcmc.diag.plot(chain)

# make matrix of 1D and 2D density plots
cont.pairs(chain$theta, prob.levels=c(0.683, 0.9, 0.99), smooth1d=TRUE)
