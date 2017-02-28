# tonic

Tonic is a collection of pure R tools for generating and manipulating MCMC output. 

The current version includes the top-level functions:

 * gw.mcmc    - sample a target density using the ensemble sampler of Goodman & Weare (2010)
 * mh.mcmc    - sample a target density using the random walk M-H method 
 * contour_matrix - plot an NxN matrix summarising pairwise relationships between variables
 * mcmc.diag.plots - make diagnostic plots of chain outputs

## Installation

Tonic is an R package, but is still in development. To set up from GitHub first install (if you haven't already) Hadley Wickham's devtools package.
```
   install.packages("devtools")
```
Now you can install tonic straight from GitHub:
```
   devtools::install_github("svdataman/tonic")
```
It will also install the mvtnorm package which it depends on. Now you're good to go.

## Getting started

Let's define an R function that takes a vector of parameters as its first argument, and returns a (scalar) log density (up to some normalisation constant)

```R
my_posterior <- function(theta) {
  cov <- matrix(c(1,0.98,0.8,0.98,1.0,0.97,0.8,0.97,2.0), nrow = 3)
  logP <- mvtnorm::dmvnorm(theta, mean = c(-1, 2, 0), sigma = cov, log = TRUE)
  return(logP)
}
```
Now we can generate a sample using e.g.

```R
   chain <- tonic::gw.mcmc(my_posterior, theta.0 = c(0,0,0), nsteps = 1e4)
```

and plot the result using

```R
  tonic::contour_matrix(chain$theta)
```

For more help on each comment, try the in-built help, e.g.
```R
   ?gw.mcmc
   ?contour_matrix
```

Below is an example for a four-variable problem.

![example](figures/ContPairs_test.png)

## To do

* complete documentation of functions
* annealing? (with evidence calculation)
* evidence from nested sampling?
* periodic saving of output
* configure input/output so that can pick-up chain(s) where left off.
* add Laplace evidence function, MAP finder
 
## Referencing tonic

If you find tonic useful in your work, please cite the following paper for
which tonic was developed:

[S. Vaughan et al., 2016, MNRAS, v461, pp3145-3152](http://adsabs.harvard.edu/abs/2016MNRAS.461.3145V)
