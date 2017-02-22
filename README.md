# tonic

Tonic is a collection of pure R tools for generating and manipulating MCMC output. 

The current version includes the top-level functions:

 gw.mcmc    - sample a target density using the ensemble sampler of Goodman & Weare (2010)
 
 mh.mcmc    - sample a target density using the random walk M-H method 

 cont.pairs - plot an NxN matrix summarising pairwise relationships between variables
 
 mcmc.diag.plots - make diagnostic plots of chain outputs

## Installation

Tonic is an R package, but is still in development. To set up from GitHub first install (if you haven't already) Hadley Wickham's devtools package.
```
   install.packages("devtools")
```
Now you can install tonic straight from GitHub:
```
   devtools::install_github("svdataman/tonic")
```
It will also install the mvtnorm package which it depends on. Now, load into your R session with
```
   require(tonic)
```
and you're good to go.

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

## To do

* annealing? (with evidence calculation)
* evidence from nested sampling?
* periodic saving of output
* configure input/output so that can pick-up chain(s) where left off.
* add Laplace evidence function, MAP finder
 
