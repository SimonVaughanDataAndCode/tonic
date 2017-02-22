# ------------------------------------------------
# a mixture of two square (uniform) distribution with different masses.
# This is a tough test: although the distribution have very simple shapes,
# they are separated by a 'valley' of zero probability, that is hard for a
# chain to cross.

my.pdf <- function(theta) {
  n <- 2
  xlim <- c(0,1)
  ylim <- c(0,1)
  P1 <- (theta[1] <= xlim[2]) & (theta[1] >= xlim[1]) &
    (theta[2] >= ylim[1]) & (theta[2] <= ylim[2])
  xlim <- c(0,1) + 1.5
  ylim <- c(0,1) + 1.5
  P2 <- (theta[1] <= xlim[2]) & (theta[1] >= xlim[1]) &
    (theta[2] >= ylim[1]) & (theta[2] <= ylim[2])
  logP <- log(0.2*P1 + 0.8*P2)
  return(logP)
}

# ------------------------------------------------
# test
# In this case we need nsteps >~ 2e5 or the second mass of probability is
# not properly explored.

chain <- gw.mcmc(my.pdf, theta.0=c(0.5,0.5), nsteps=2e5,
                  chatter=1, walk.rate=5, atune=2) #, merge.walkers=FALSE)

plot(chain$theta[,1], chain$theta[,2], xlim = c(0,3), ylim = c(0,3))
hist(chain$theta[,1])