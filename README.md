# tonic
MCMC with R

Tonic is a collection of pure R tools for generating and manipulating MCMC output. 

The current version includes:

 gw.mcmc    - sample a target density using the ensemble sampler of Goodman & Weare (2010)

 cont.pairs - plot an NxN matrix summarising pairwise relationships between variables

## Installation

Tonic is not (yet) and R package. To set up the R functions source the .R files

```R
source("gwmcmc.R")
source("plot_contour.R")
```

## Getting started

Given an R function that takes a vector of parameters as its first argument, and returns a (scalar) log density (up to some normalisation constant), we can generate a sample using e.g.

```R
theta <- gw.mcmc(my.pdf, theta.0 = c(0,0,0), nsteps = 1e5)
```

and plot the result using

```R
  cont.pairs(theta)
```

Inside the .R scripts is more detailed information on the input parameters and output formats.

## Example output

Use the .R scripts in the tests directory to check these work. Below is an example using

```R
  cont.pairs(theta, prob.levels = c(1, 2), cex = 0.5, npix = 100, sigma = TRUE,
              smooth1d = TRUE, cex.lab = 1.3)
```

![example](figures/ContPairs_test.png)
