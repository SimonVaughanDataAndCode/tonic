# ------------------------------------------------


# a 3D un-normalised distribution comprising 3 separate Normal densities
# separated in all 3 dimensions

my.pdf <- function(theta) {
  n <- 3
  cov <- diag(n) * 0.03
  P1 <- mvtnorm::dmvnorm(theta, mean=c(0,0,0), sigma=cov)
  cov <- diag(n) * 0.03
  P2 <- mvtnorm::dmvnorm(theta, mean=c(2,2,2), sigma=cov)
  cov <- diag(n) * 0.03
  P3 <- mvtnorm::dmvnorm(theta, mean=c(3,-2,3), sigma=cov)
  logP <- log(0.5*P1 + 0.5*P2 + 0.5*P3)
  return(logP)
}

# ------------------------------------------------
# test

chain <- gw.mcmc(my.pdf, theta.0=c(0,0,0), nsteps=10e4, cov.init = diag(3)*2,
                  chatter=1, burn.in = 1e4, walk.rate=5) #, merge.walkers=FALSE)

mcmc.diag.plot(chain)
theta <- chain$theta

x <- c(0,0,0)
print(mean(sqrt((theta[,1]-x[1])^2+(theta[,2]-x[2])^2+(theta[,3]-x[3])^2) < 2))
x <- c(2,2,2)
print(mean(sqrt((theta[,1]-x[1])^2+(theta[,2]-x[2])^2+(theta[,3]-x[3])^2) < 2))
x <- c(3,-2,3)
print(mean(sqrt((theta[,1]-x[1])^2+(theta[,2]-x[2])^2+(theta[,3]-x[3])^2) < 2))

#source('~/R/tonic/plot_contour.R')
cont.pairs(theta, prob.levels=c(0.683, 0.9, 0.99), smooth1d=TRUE)
