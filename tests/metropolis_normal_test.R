# ------------------------------------------------


# a 3D un-normalised distribution comprising 3 separate Normal densities
# separated in all 3 dimensions


my.pdf <- function(theta) {
  n <- 3
  cov <- diag(n) * 0.03
  P1 <- mvtnorm::dmvnorm(theta, mean=c(0,0,0), sigma=cov)
  cov <- diag(n) * 0.03
  P2 <- mvtnorm::dmvnorm(theta, mean=c(1,1,1), sigma=cov)
  cov <- diag(n) * 0.03
  P3 <- mvtnorm::dmvnorm(theta, mean=c(1,-1,-1), sigma=cov)
  logP <- log(0.5*P1 + 0.5*P2 + 0.5*P3)
  return(logP)
}

# ------------------------------------------------
# test MCMC function mh.mcmc

chain <- mh.mcmc(my.pdf, theta.0=c(0,0,0), nsteps=5e4, cov = diag(3),
                  chatter=1, burn.in = 1e4)

# examine diagnostic plots
mcmc.diag.plot(chain)

# make matrix of 1D and 2D density plots
#source('~/R/tonic/plot_contour.R')
cont.pairs(chain$theta, prob.levels=c(0.683, 0.9, 0.99), smooth1d=TRUE)
