# ------------------------------------------------


# a 3D un-normalised distribution comprising 3 separate Normal densities
# separated in all 3 dimensions

my.pdf <- function(theta) {
  n <- 2
  xlim <- c(0,1)
  ylim <- c(0,1)
  P1 <- (theta[1] <= xlim[2]) & (theta[1] >= xlim[1]) &
    (theta[2] >= ylim[1]) & (theta[2] <= ylim[2])
  xlim <- c(0,1) + 2
  ylim <- c(0,1) + 2
  P2 <- (theta[1] <= xlim[2]) & (theta[1] >= xlim[1]) &
    (theta[2] >= ylim[1]) & (theta[2] <= ylim[2])
  logP <- log(0.2*P1 + 0.8*P2)
  return(logP)
}

# ------------------------------------------------
# test

theta <- gw.mcmc(my.pdf, theta.0=c(0.5,0.5), nsteps=5e4,
                  chatter=1, walk.rate=1, atune=2) #, merge.walkers=FALSE)

plot(theta[,1], theta[,2], xlim = c(0,3), ylim = c(0,3))
