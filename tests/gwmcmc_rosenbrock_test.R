# ------------------------------------------------
# the 2D Rosenbrock density

rosenbrock <- function(theta) {
  logP <- -5 * (theta[2] - theta[1]^2)^2 - 
           (1 / 20) * (1 - theta[1])^2
  return(logP)
}


# ------------------------------------------------
# plot contours of density

  x <- seq(-7, 7, by = 0.02)
  y <- seq(-5, 50, by = 0.1)
  z <- array(NA, dim = c(length(x), length(y)))
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      z[i, j] <- rosenbrock(c(x[i], y[j]))
    }
  }
  
  plot(0, 0, xlim = range(-5, 6), ylim = range(-2, 35),
       main = "", type = "n", bty = "n")
  contour(x, y, exp(z), levels = c(0.1,0.3,0.8), add = TRUE,
          drawlabels = FALSE)
  
# ------------------------------------------------
# test MCMC
  
# generate a chain
  chain <- gw.mcmc(rosenbrock, theta.0 = c(0,0), 
                   nsteps = 1e5, burn.in = 1e5,
                   chatter = 1, walk.rate = 5,
                   thin = 1)
  
# use diagnostic plots
  mcmc.diag.plot(chain)

# plot the points
  layout(1)
  plot(chain$theta[,1], chain$theta[,2], bty = "n", 
       pch = 1, cex = 0.5)

  cont.pairs(chain$theta)
  
# ------------------------------------------------
# use the plotting function from plot_contour.R

# plot density contours of chain output
  plot.cont(chain$theta[,1], chain$theta[,2], npix = 100, smooth2d = FALSE,
            xlim = range(-5, 6), ylim = range(-2, 35),
            prob.levels = c(0.9, 0.99))
  axis(1)
  axis(2)

# overlay true density contours
  contour(x, y, exp(z), levels = c(0.1, 0.3, 0.8), 
          col = "red", add = TRUE,
          drawlabels = FALSE, lwd = 1)
  
# ------------------------------------------------
  cat('-- Effective sample size (from integrated autocorrelation).', fill = TRUE)
  print( coda::effectiveSize(chain$theta) )
  
  mcmc.diag.plot(chain)
  