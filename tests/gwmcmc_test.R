# ------------------------------------------------


# a 3D un-normalised distribution comprising 3 separate Normal densities
# separated in all 3 dimensions

my.pdf <- function(theta) {
  cov <- matrix(c(1,0.98,0.8,
                  0.98,1.0,0.97,
                  0.8,0.97,2.0), nrow = 3)
  logP <- mvtnorm::dmvnorm(theta, mean=c(-1,2,0), 
                         sigma = cov, log = TRUE)
  return(logP)
}

# ------------------------------------------------
# load the MCMC and plotting functions
source("gwmcmc.R")
source("plot_contour.R")
source("diagnostics.R")

# ------------------------------------------------
# test

chain <- gw.mcmc(my.pdf, theta.0=c(0,0,0), nsteps=10e4, burn.in = 1e4)

mcmc.diag.plot(chain)

cont.pairs(chain$theta, prob.levels=c(0.683, 0.9, 0.99), smooth1d=TRUE)
