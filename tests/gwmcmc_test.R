# ------------------------------------------------


# a 3D un-normalised distribution comprising 3 separate Normal densities
# separated in all 3 dimensions

my.pdf <- function(theta) {
  n <- 3
  cov <- diag(n) * 0.05
  P1 <- mvtnorm::dmvnorm(theta, mean=c(0,0,0), sigma=cov)
  cov <- diag(n) * 0.05
  P2 <- mvtnorm::dmvnorm(theta, mean=c(2,2,2), sigma=cov)
  cov <- diag(n) * 0.05
  P3 <- mvtnorm::dmvnorm(theta, mean=c(3,-2,3), sigma=cov)
  logP <- log(0.5*P1 + 0.5*P2 + 0.5*P3)
  return(logP)
}

# ------------------------------------------------
# test

theta <- gw.mcmc(my.pdf, theta.0=c(0,0,0), nsteps=4e5, cov.init = diag(3)*2,
                  chatter=1, walk.rate=4)

#plot(result[,1], -result[,2]-result[,1])
#plot(theta[,1], theta[,2])

x <- c(0,0,0)
print(mean(sqrt((theta[,1]-x[1])^2+(theta[,2]-x[2])^2+(theta[,3]-x[3])^2) < 2))
x <- c(2,2,2)
print(mean(sqrt((theta[,1]-x[1])^2+(theta[,2]-x[2])^2+(theta[,3]-x[3])^2) < 2))
x <- c(3,-2,3)
print(mean(sqrt((theta[,1]-x[1])^2+(theta[,2]-x[2])^2+(theta[,3]-x[3])^2) < 2))

source('~/R/PlotContourPairs/plot_contour.R')
cont.pairs(theta, prob.levels=c(0.683, 0.9, 0.99), smooth1d=TRUE)
