# --------------------------------------------------
# plot_contour.R
# Collection of R functions for plotting an NxN matrix of
# contour or scatter plots. A little like the pairs()
# function in R for scatter plots, this plots contours of 
# the 2D density (in lower right panels) and 1D density
# (histogram or smoothed density estimate) on diagonals
# --------------------------------------------------

# --------------------------------------------------
#
# History:
#  04/04/16 - First working version
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
#
# --------------------------------------------------

#' Interpolate an even grid of z(x,y) onto arbitrary (x,y) points
#' 
#' \code{interp_image} 
#' 
#' From an image \code{z(x,y)} at evenly spaced grid points \code{x = i*dx} 
#' \code{y = j*dy}, return the values of \code{z} at arbirary points
#' \code{(x, y)} between grid points. The values are calculated using bilinear
#' interpolation. I.e. for an output point \code{(x, y)} the corners of the
#' surrounding grid square are identified, and the result is a weighted
#' combination of these four \code{z} values: \code{z(x1, y1)}, \code{z(x2, y1)},
#' \code{z(x1, y2)}, \code{z(x2, y2)} where the weights depend on the distances
#' from the \code{x} and \code{y} boundaries of the grid square.
#' 
#' @param img (list) containing three components: \code{x, y, z}
#' @param x   (vector) \code{x} positions at which to compute \code{z}
#' @param y   (vector) \code{y} positions at which to compute \code{z}
#'
#' @return
#'  The output is a vector of \code{z} values at the position of each
#'  input \code{(x, y)} coordinate.
#'  
#' @section Notes:
#' The input \code{img} should be a list with three components:
#'  \code{img$z} is a vector of \code{nx} evenly spaced \code{x} positions,
#'  \code{img$y} is a vector of \code{ny} evenly spaced \code{y} position,
#'  \code{img$z} is an \code{nx*ny} array of \code{z} values at the 
#'  \code{(x,y)} coordinates.
#'  
#' @seealso \code{\link{chain_convergence}}, \code{\link{contour_matrix}}
#'
#' @export
interp_image <- function(image, x, y) {
  
  # check the inputs
  
  if (missing(image)) stop('Missing IMAGE input.')
  if (missing(x)) stop('Missing X input vector.')
  if (missing(y)) stop('Missing Y input vector.')
  
  if (!exists("z", where=image)) stop('IMAGE input does not contain z array.')
  if (!exists("x", where=image)) stop('IMAGE input does not contain x vector.')
  if (!exists("y", where=image)) stop('IMAGE input does not contain y vector.')
  
  if (length(x) != length(y)) stop('X and Y lengths differ.')
  
  # parameters of the image
  
  dx <- diff(image$x[1:2])
  dy <- diff(image$y[1:2])
  nx <- NROW(image$z)
  ny <- NCOL(image$z)
  x0 <- image$x[1]
  y0 <- image$y[1]
  
  # find pixels nearest to the points (x,y).
  # x is between image$x[pix.x1) and image$x[pix.x2]
  # y is between image$y[pix.y1) and image$y[piy.x2]
  
  pix.x1 <- floor( (x - x0) / dx ) + 1
  pix.y1 <- floor( (y - y0) / dy ) + 1
  pix.x2 <- floor( (x - x0) / dx ) + 2
  pix.y2 <- floor( (y - y0) / dy ) + 2
  
  pix.x2 <- pmin( pix.x2, nx )
  pix.y2 <- pmin( pix.y2, ny )
  
  # weighting factors for the weighted sum
  
  w.x <- (x-image$x[pix.x1])/dx 
  w.y <- (y-image$y[pix.y1])/dy
  
  # value of the image at each of the four corners
  
  z11 <- image$z[pix.x1 + nx*(pix.y1-1)]
  z12 <- image$z[pix.x1 + nx*(pix.y2-1)]
  z21 <- image$z[pix.x2 + nx*(pix.y1-1)]
  z22 <- image$z[pix.x2 + nx*(pix.y2-1)]
  
  # the weighted sum of the four corner values
  
  zout <- (1-w.x)*(1-w.y)*z11 + 
              (1-w.x)*(w.y)*z12 + 
              (w.x)*(1-w.y)*z21 + 
              (w.x)*(w.y)*z22
  
  # return the vector of interpolated values
  
  return(zout)
}

# --------------------------------------------------
# hist2d
#
# History:
#  04/04/16 - First working version
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
#
# --------------------------------------------------

#' Bin a set of (x,y) points into a 2D histogram
#' 
#' \code{hist2d} produce a 2D histogram from \code{x, y} points.
#' 
#' From an set of \code{(x, y)} positions, compute a 2D histogram.
#' Count the number of input points falling within each pixel on
#' and evenly spaced grid points of \code{pix.x = i*dx}, \code{pix.y = j*dy}. 
#' return this values as \code{z}. 
#' 
#' @param npix number of bins on each axis (pixels = \code{npix^2})
#' @param lims Four-element vector for \code{x, y} limits of the histogram
#' @param prob (logical) make the area under the histogram equal to 1?
#' @inheritParams contour_matrix
#'
#' @return
#'  List containing components \code{x, y, z}.
#'  
#' @section Notes:
#' The output \code{img} is a list with three components:
#'  \code{z} is a vector of \code{nx} evenly spaced \code{x} positions,
#'  \code{y} is a vector of \code{ny} evenly spaced \code{y} position,
#'  \code{z} is an \code{nx*ny} array of \code{z} values at the 
#'  \code{(x, y)} coordinates.
#'  
#' @seealso \code{\link{chain_convergence}}, \code{\link{contour_matrix}}
#' 
#' @examples 
#' n <- 1000
#' x <- rnorm(n)
#' y <- x + rnorm(n)
#' img <- hist2d(x, y, npix=11)
#' image(img$x, img$y, log10(img$z))
#' points(x,y)
#'
#' @export
hist2d <- function(x, y, 
                   npix = 25, 
                   lims = c(range(x), range(y)),
                   prob = FALSE) {
  
  if (missing(x)) stop('Missing X input vector.')
  if (missing(y)) stop('Missing Y input vector.')
  if (length(x) != length(y)) stop('X and Y have different lengths.')
  
  # determine x, y limits of the 2D histogram
  
  xlim <- lims[1:2]
  ylim <- lims[3:4]
  
  # determine resolution (pixel size) of the 2D histogram
  
  dx <- (xlim[2] - xlim[1]) / (npix - 1)
  dy <- (ylim[2] - ylim[1]) / (npix - 1)
  
  # define the pixel locations in term of x and y
  
  x.bin <- seq(xlim[1]-dx/2, xlim[2]+dx/2, by=dx)
  y.bin <- seq(ylim[1]-dy/2, ylim[2]+dy/2, by=dy)
  
  # for each point (x,y) find which pixel it falls into
  
  x.pix <- floor( (x-x.bin[1]) / dx ) + 1
  y.pix <- floor( (y-y.bin[1]) / dy ) + 1

  # compute contingency table, showing counts in pixel x.bin, y.bin
  # Note: we need the explicit factor conversion here, otherwise
  # rows/columns with zero counts will be dropped from the table.
  
  freq2D <- table( factor(x.pix, levels=1:npix), factor(y.pix, levels=1:npix) )
    
  # if prob is TRUE then normalise so that thr area under the
  # histogram is 1.
  
  norm <- sum(freq2D)*dx*dy
  if (prob == TRUE)  freq2D <- freq2D / norm

  # freq2D values are for the middle of the bins. Shift the x.bin
  # y.bin positions to the bin centre.
  
  x.bin <- x.bin + dx/2
  y.bin <- y.bin + dy/2

  # output the pixel x and y ordinates, and the 2D histogram (z)

  img <- list(x=x.bin[1:npix], y=y.bin[1:npix], z=freq2D)
  
  return(img)
  
}

#
# History:
#  04/04/16 - First working version
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
#
# ---------------------------------------

#' Compute contours encircling specified masses.
#' 
#' \code{get_levels2d} returns density levels enclosing specified masses. 
#' 
#' Function to compute the density levels at which to draw contour lines in
#' order to enclose a fraction \code{prob} of the mass.
#' 
#' @param prob (vector) probability masses at which to define levels.
#' @inheritParams interp_image
#'
#' @return
#'  Vector of the density levels (one for each element of \code{prob}). 
#'  If \code{prob = 0.95} then a contour a level contains 95\% of the mass.
#'  
#' @section Notes:
#' The input \code{img} should be a list with three components:
#'  \code{img$z} is a vector of \code{nx} evenly spaced \code{x} positions,
#'  \code{img$y} is a vector of \code{ny} evenly spaced \code{y} position,
#'  \code{img$z} is an \code{nx*ny} array of \code{z} values at the 
#'  \code{(x, y)} coordinates.
#'  
#' @seealso \code{\link{chain_convergence}}, \code{\link{contour_matrix}}
#' 
#' @examples
#' n <- 1000
#' x <- rnorm(n)
#' y <- x + rnorm(n)
#' img <- hist2d(x, y, npix=11)
#' get_levels2d(img, prob = c(0.5, 0.9))
#'
#' @export
get_levels2d <- function(img, prob=0.95) {
  
  # check the inputs

  if (missing(img)) stop('Missing data IMG.')

  if (!exists("z", where=img)) stop('IMG input does not contain z array.')
  if (!exists("x", where=img)) stop('IMG input does not contain x vector.')
  if (!exists("y", where=img)) stop('IMG input does not contain y vector.')
  
  # probability should be between 0 and 1  
  
  if (max(prob) > 1) stop('PROB > 1 error.')
  if (min(prob) < 0) stop('PROB < 0 error.')
  
  dx <- diff(img$x[1:2])
  dy <- diff(img$y[1:2])
  norm <- sum(img$z) 
  sz <- sort(img$z)
  c1 <- cumsum(sz) * dx * dy
  clevels <- approx(c1, sz, xout = 1 - prob)$y
  
  return( clevels )
}

# ------------------------------------------------
#
# History:
#  04/04/16 - First working version
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
#
# -------------------------------------------------

#' Find lower, upper intervals enclosing certain masses of a histogram.
#' 
#' \code{get_levels1d} returns values bounding specified masses.
#' 
#' Function to compute the lower and upper limits of intevals containing a
#' fraction (\code{prob}) of the mass of the points \code{x}. If \code{x} are
#' input data points, we compute the lower and upper quantiles. For for
#' \code{prob = 0.95} we find the 0.025 and 0.975 quantiles, between which are 
#' 95\% of the data points.
#' 
#' @param x (vector) x data values.
#' @inheritParams get_levels2d
#'
#' @return
#'  Vector of the x values (two for each element of \code{prob}). 
#'  If \code{prob = 0.95} then the lower, upper values enclose 95\% of the mass
#'  of x.
#'  
#' @seealso \code{\link{chain_convergence}}, \code{\link{contour_matrix}}
#' 
#' @examples
#' n <- 1000
#' x <- rnorm(n)
#' get_levels1d(x, prob = c(0.5, 0.9))
#'
#' @export
get_levels1d <- function(x, prob = 0.95) {

  # check the inputs
  if (missing(x)) stop('Missing X input vector.')
  if (length(x) < 2) stop('Not enough data in X')
  if (max(prob) > 1) stop('PROB > 1 error.')
  if (min(prob) < 0) stop('PROB < 0 error.')
  
  prob.levels <- c( (1-prob)/2, 1-(1-prob)/2)
  prob.levels <- sort(prob.levels)
  quants <- quantile(x, probs=prob.levels)
  return( quants )
}

# --------------------------------------------------
#
# History:
#  04/04/16 - First working version
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
# --------------------------------------------------

#' Produce a 2D density/histogram plot given data values \code{x, y}.
#' 
#' \code{plot_density_contours} generates a 2D density/histogram plot.
#' 
#' Given vectors \code{x, y} specifying the coordinates of data points,
#' \code{plot_density_contours} generates a 2D density plot. What is plotted are contours
#' enclosing specified masses of the distribution, and points outside of the
#' outermost contour are marked individually. The contours are based on
#' either a smooth 2D kernel estimate or an a (rough) 2D histogram.
#' 
#' @param x,y (vectors) x and y data values (equal length vectors).
#' @param xlim,ylim (vectors) The (x, y) limits of the axes (x1, x2), (y1, y2).
#' @param plot.dots (logical) Plot points outside of outer contour?
#' @param plot.contours (logical) Plot the contours?
#' @param xlab,ylab (strings) Labels for x, y axes.
#' @param pch (integer) Plotting 'character', i.e. symbol to use. Same as 
#'         'base' graphics.
#' @inheritParams contour_matrix
#'
#' @return
#'  None.
#'  
#' @seealso \code{\link{chain_convergence}}, \code{\link{contour_matrix}}
#' 
#' @examples
#' n <- 1000
#' x <- rnorm(n)
#' y <- x + rnorm(n)
#' img <- hist2d(x, y, npix = 40)
#' plot_density_contours(x,y, xlim = c(-4,4), ylim = c(-4, 4))
#' axis(1); axis(2)
#'
#' @export
plot_density_contours <- function(x, y, 
                      xlim, ylim, 
                      npix = 25, 
                      prob.levels = c(0.683, 0.90),
                      plot.dots = TRUE,
                      plot.contours = TRUE,
                      plot.image = FALSE,
                      dot.level = NULL,
                      jittr = FALSE,
                      col = "black",
                      pch = 1,
                      cex = 1.0, 
                      smooth2d = TRUE, 
                      xlab = "",
                      ylab = "",
                      ...) {
  
  # set up the plot window
  
  plot(1, 1, type="n", main="", axes=FALSE, xlim=xlim, ylim=ylim,
       xlab = xlab, ylab = ylab, ...)
  
  # compute the 2D smoothed density, and density levels at which 
  # to draw the contours. Make sure the (X,Y) range of the smoothed
  # "image" is large enough to include all data.
  
  lims <- c(min(xlim[1], min(x)), max(xlim[2], max(x)),
            min(ylim[1], min(y)), max(ylim[2], max(y)))
  
  if (smooth2d == TRUE) {
    z <- MASS::kde2d(x,y, n=npix, lims=lims)
  } else {
    z <- hist2d(x, y, npix=npix, lims, prob=TRUE)
  }
  
  clevels <- get_levels2d(z, prob=prob.levels)
  
  # add contours to the plot
  if (plot.image == TRUE) image(z$x, z$y, sqrt(z$z), add=TRUE, col=terrain.colors(20))
  if (plot.contours == TRUE) contour(z, levels=clevels, drawlabels=FALSE, add=TRUE) 
  
  # find the density at the position of each points
  
  if (plot.dots == TRUE) {
    density.xy <- interp_image(z, x, y)
    
    # define the contour outside of which to draw individual points
    
    if (is.null(dot.level)) { dot.level <- length(clevels) }
    dot.level <- min(dot.level, length(clevels))
    
    # find which points lie outside the dot.level contour
    
    if (dot.level > 0) mask <- (density.xy < clevels[dot.level])
    if (dot.level == 0) mask <- 1:length(x)
    
    # mask = where (x,y) is outside of the dot.level contour 
    # add these as points to the plot
    
    x.p <- x[mask]
    y.p <- y[mask]
    if (jittr == TRUE) {
      x.p <- jitter(x.p)
      y.p <- jitter(y.p)
    }
    points(x.p, y.p, pch=pch, cex=cex, col=col)
  }
  
}

# --------------------------------------------------
# NAME: 
#     plot_density
#
# PURPOSE:
#     Produce a single histogram or density plot
#
# AUTHOR:
#     Simon Vaughan
#
# CALLING SEQUENCE:
#     plot_density(x)
#
# INPUTS:
#   x         - vector of values
#
# OPTIONAL INPUTS:
#   breaks    - specify the break positions for a histogram
#   xlim      - 2 element vector of limits for x axis
#   prob      - mark probability intervals (NULL or vector)
#   smooth    - FALSE = plot histogram, TRUE = plot smooth density
#   fill.col  - colour for filling histogram/density plot
#   xlab, xlab, main - plot titles, same as plot()
#
# OUTPUT:
#   clevels   - the quantiles for the input probability levels
#
# History:
#  04/04/16 - v0.1 - First working version
#  28/04/16 - v0.2 - Fixed bug in quantile estimation
#  16/06/16 - v0.3 - Return NULL when prob not given at input
#                     added inputs for xlab, ylab, main and '...'.
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
# --------------------------------------------------

#' Produce a 1D density/histogram plot given data values \code{x}.
#' 
#' \code{plot_density} generates a 1D density/histogram plot.
#' 
#' Produce a single histogram or density plot, for use with \code{contour_matrix}.
#' 
#' @param x (vector) x data values.
#' @param xlim (vector) The limits of the x axis (x1, x2).
#' @param xlab,ylab (strings) Labels for x, y axes.
#' @inheritParams contour_matrix
#'
#' @return
#'  None.
#'  
#' @seealso \code{\link{chain_convergence}}, \code{\link{contour_matrix}}
#' 
#' @examples
#' n <- 1000
#' x <- rnorm(n)
#' y <- x + rnorm(n)
#' img <- hist2d(x, y, npix = 40)
#' junk <- plot_density(x, xlim = c(-4,4))
#' axis(1)
#'
#' @export
plot_density <- function(x, 
                      breaks = 30, 
                      xlim, 
                      prob = NULL, 
                      smooth = FALSE, 
                      fill.col = NULL,
                      plot.lines = TRUE,
                      xlab = "",
                      ylab = "",
                      main = "",
                      ...) {

  # compute a histogram and a smooth density curve
  
  h <- hist(x, breaks=breaks, plot=FALSE)
  s <- density(x, n=2^14)
  
  # prepare the plot window
  
  ylim <- c(0, max(h$density))
  if (smooth == TRUE) ylim[2] <- max(s$y)
  plot(0, 0, type = "n", xlab = xlab, ylab = ylab, main = main, xaxt = "n", 
       yaxt = "n", bty = "n", xlim = xlim, ylim = ylim, ...)
  
  # plot histogram or density curve

  nbin <- length(h$mids)
  if (smooth == FALSE) {
    x.plot <- c(xlim[1], rep(h$breaks, each=2), xlim[2])
    y.plot <- c(0, 0, rep(h$density, each=2), 0, 0)
    polygon(x.plot, y.plot, col = fill.col)
    lines(c(h$breaks[1], h$breaks[1:nbin], h$breaks[nbin+1]),
          c(0, h$density, 0), type="s", lwd=1)
  } else {
    if (!is.null(fill.col))
        polygon(c(xlim[1], s$x, xlim[2]), c(0, s$y, 0), col=fill.col)
    lines(s$x, s$y, lwd=1)
  }

  # mark mean and 1D intervals on histogram/density
  
  if (!is.null(prob)) {
    nlevels <- length(prob)
    clevels <- get_levels1d(x, prob=prob)
    if (plot.lines == TRUE) {
      abline(v=mean(x), lty=1, lwd=2, col="grey80")
      abline(v = clevels, lty=2, lwd=2, col="grey80")
    }
    clevels <- c( clevels[1:nlevels], mean(x), 
                  clevels[(nlevels+1):(2*nlevels)] )
  } else {
    clevels <- NULL
  }
  
  return(clevels)
}

# --------------------------------------------------
#
# History:
#  04/04/16 - First working version
#
# Simon Vaughan, University of Leicester
# Copyright (C) 2016 Simon Vaughan
# --------------------------------------------------

#' Produce a matrix of contours plots from a data array
#' 
#' \code{contour_matrix} generates an \code{M*M} matrix of scatter/density plots.
#'                  
#' Given an array of data with \code{M} columns (variables) and \code{N} rows 
#' (observations), produce a matrix of plots showing contours for each pair of 
#' parameters. The \code{pairs()} function shows a matrix of scatter plots for 
#' each pair of variables, here the scatter plots are replaced by contour plots.
#' The contours are produced by using \code{MASS:kde2d} to produce a 2D Gaussian
#' kernal density estimate, and finding the contours that encolse a fraction
#' (prob) of the total density.
#'
#' @param theta        (array) \code{N} by \code{M} array of parameter values
#' @param ranges      (array) \code{M} by \code{2} array of ranges for the plots
#' @param cex         (float) character expansion factor.
#' @param prob.levels (float) probability levels, e.g. to plot contours enclosing
#'                       90\% and 95\% of the mass on 2D distributions
#' @param labels      (array of strings) names of the variables
#' @param upper       (string) plot also in upper-right triangle. See notes.
#' @param dot.level   (integer) draw dots outside of which contour? (1,2,...)
#' @param breaks      (integer) how many bins for 1D histograms
#' @param smooth1d    (logical) plot smoothed distributions, rather than 
#'                     histograms on the leading diagonal?
#' @param smooth2d    (logical) use smooth, kernel density estimates, rather 
#'                      than 2D histograms?                       
#' @param npix        (integer) use \code{npix*npix} grid of pixels for 
#'                      computing 2D smoothed distributions (resolution of 
#'                      contour maps)
#' @param cex.lab     (float) expansion factor for axis labels.
#' @param cex.axis    (float) expansion factor for axis names.
#' @param prob1d      (float array) probability levels at which to mark 
#'                      intervals on 1D distributions, e.g. c(0.683, 0.90).
#' @param sigma       (logical) are probability levels \code{prob1d} and 
#'                      \code{prob.levels} given in terms of sigmas?
#'                      (If \code{FALSE} - the default - then these are 
#'                      probabilities.)
#' @param jittr       (logical) add 'jitter' to points to reduce overlap.
#' @param thin        (integer) a factor by which to 'thin out' data input data
#'                      i.e. plot only a fraction \code{1/thin} of the data 
#'                      points.
#' @param col         (string) colour of the data points.
#' @param plot.image  (logical) Plot a colour image underneath the contours?
#' @param plot.1dlines (logical) Mark intervals on 1D plots (histogram/density)?
#' @param fill.col    (string) Colour to use under the histogram.
#' @param ...         (anything) any other graphical keywords to be passed to 
#'                        \code{plot(...)}
#'   
#' @return Array of the 1D intervals for each parameter, or NULL.
#' 
#' @section Notes:
#' The input \code{chain} should be a list such as produced by 
#'   \code{gw_sampler} or \code{mh_sampler} that contains the following:
#' \describe{ 
#' \item{theta}{(array) n * ndim array of posterior samples
#'             n samples of ndim vectors of parameters}
#' \item{method}{(string) name of MCMC method used}
#' \item{nwalkers}{number of walkers used (if \code{method = 'gw.mcmc'})}
#' \item{nchains}{number of walkers used (if \code{method = 'mh.mcmc'})}
#' }
#' 
#' The \code{upper} parameter determines what is to be shown in the upper-right
#' corner of the plot. The options are: \itemize{
#'                      \item{\code{"NA"} - no plot} 
#'                      \item{\code{"points"} - scatter plot of points}
#'                      \item{\code{"image"} - intensity map} 
#'                      \item{\code{"contour"} - contour map}
#'                      }
#'                      
#' You may notice that the fraction of points plotted outside the 
#' \code{dot.levels} contour line is often smaller than
#' \code{prob.levels[dot.levels]}. This is usually true for finite \code{N} and
#' is a natural consequence of plotting contours that enclose a fraction of the
#' smoothed density rather than a fraction of the points.
#' 
#' @seealso \code{\link{gw_sampler}}, \code{\link{mh_sampler}}, 
#   \code{\link{get_levels2d}}, \code{\link{get_levels1d}}, 
#   \code{\link{interp_image}}, \code{\link{hist2d}}, 
#   \code{\link{plot_density}}, \code{\link{plot_density_contours}}
#' 
#' @examples 
#' n <- 10000
#' x <- rchisq(n, df=2)
#' y <- x + rnorm(n)
#' z <- rnorm(n, sd=1)
#' dat <- cbind(x,y,z)
#' contour_matrix(dat, smooth1d=TRUE)
#' contour_matrix(dat, prob.levels=c(1,2), cex=0.5, npix=100, sigma=TRUE,
#'           dot.level=2, prob1d=c(1, 2), smooth1d=TRUE, cex.lab=1.3,
#'           plot.1dlines = TRUE)
#'
#' @export
contour_matrix <- function(theta, 
                       ranges = NULL, 
                       cex = 1.0, 
                       prob.levels = c(0.9, 0.95),
                       labels = colnames(theta), 
                       upper = NA, 
                       dot.level = NULL,
                       breaks = 30, 
                       smooth1d = FALSE, 
                       smooth2d = TRUE,
                       npix = 100, 
                       cex.lab = 1, 
                       cex.axis = 1, 
                       prob1d = NULL, 
                       sigma = FALSE, 
                       jittr=FALSE, 
                       thin = NULL, 
                       pch = 1,
                       col = "black", 
                       plot.image = FALSE, 
                       plot.1dlines = FALSE,
                       fill.col = "steelblue3", ...) {
  
  # check the inputs

  if (missing(theta)) stop('Missing theta input data array.')

  # get the number of columns (M) and rows (n)
  
  M <- NCOL(theta)
  n <- NROW(theta)
  
  if (M < 2) stop('** Only one column in your data array.')
  if (n < 10) stop('** Too few points in the input data array.')

  # probability should be between 0 and 1  
  
  if (sigma == FALSE) {
    if (max(prob.levels) > 1) stop('** PROB.LEVELS > 1 error.')
    if (min(prob.levels) < 0) stop('** PROB.LEVELS < 0 error.')
  }
  
  # store the current graphical parameters  
  retire <- par(no.readonly = TRUE)
  
  # if requested, 'thin' the input array
  if (!is.null(thin)) {
    thin <- max(thin, 1)
    n.samp <- floor(n / thin)
    keep <- sample(1:n, n.samp)
    theta <- theta[keep,]
    n <- NROW(theta)
  }
  
  # define 1D probability levels
  if (is.null(prob1d)) {
    prob1d <- prob.levels
  }
  
  # get the ranges for each variable

  if (is.null(ranges)) {
    ranges <- array(NA, dim=c(M,2))
    for (i in 1:M) {
      ranges[i,] <- range(theta[,i])    
    }
  }
  
  if (NCOL(ranges) != 2 | NROW(ranges) != M) {
    cat('** RANGES not properly defined in contour_matrix', fill=TRUE)
  }

  # if sigma is specified on input, use sigma to (re)define the prob.levels
  # contour levels.
    
  if (sigma == TRUE) {
    prob.levels = 1 - 2 * (1 - pnorm(prob.levels))
    if (!is.null(prob1d)) {
      prob1d <- 1 - 2 * (1 - pnorm(prob1d))
    }
  }
  
  # levels1d is used to store the intervals for the 1D distribution, if used
  
  levels1d <- NULL
   
  # define a plotting array comprising M*M regions
  
  par(mfcol = c(M, M))
  
  # define a plotting array comprising M*M regions, leave some room around 
  # the edges of the plot
  
  par(mar = c(0,0,0,0), 
      oma = c(6,6,6,6), 
      mgp = c(3,1,0),
      mfcol = c(M, M))
  
  # loop over an M*M square of parameter pairs theta_i, theta_j
  
  for (i in 1:M) {      # i is the rows from top to bottom
    for (j in 1:M) {    # j is the columns from left to right
      
      x <- theta[,i]
      y <- theta[,j]
 
     # define ranges for the plot
      
      xlim <- ranges[i,]
      ylim <- ranges[j,]
      
      # if in the bottom-left triangle
 
      if (i < j) {

        plot_density_contours(x, y, xlim, ylim, npix=25, prob.levels, 
                  plot.image=plot.image,
                  dot.level=dot.level, jittr=jittr, 
                  col=col, cex=cex, pch=pch, ...)

        # add axis labels to X axis (if on bottom row, j=M) or
        # on left column (i=1).
        
        if (j == M) { mtext( labels[i], 1, line=3, cex=cex.lab ) }
        if (i == 1) { mtext( labels[j], 2, line=3, cex=cex.lab ) }
        
        # draw a box around the plot region
        
        box(which="plot")
        
        # add axes labels and tick marks only to outer
        # plots (as in the PAIRS function)
        
        if (i == 1) {
          axis(2, label=TRUE, cex.axis=cex.axis)
        } 
        if (j == M) {
          axis(1, label=TRUE, cex.axis=cex.axis)
        } 
      } 
      
      # if we're on the diagonal do not plot any data, just 
      # give the name of the ith parameter
      
      if (j == i) {
         levels1d.i <- plot_density(x, breaks, xlim, 
                                prob=prob1d, smooth = smooth1d,
                                fill.col = fill.col, plot.lines = plot.1dlines)
        
        if (i == 1) {
          levels1d <- levels1d.i
        } else {
          levels1d <- rbind(levels1d, levels1d.i)
        }
      
      # add a title/axis to the X axis if at the bottom-right panel

        if (j == M) { mtext( labels[i], 1, line=3, cex=cex.lab ) }
        if (i == M) { axis(1, labels=TRUE, cex.axis=cex.axis) }
      }
      
      # if we're in the upper right half of the plot array add
      # the points to make a scatter plot.
      
      if (i > j) {
        if (is.na(upper)) { 
          plot(1, 1, type="n", xaxt="n", yaxt="n", bty="n")
          next 
        }
 
        # set up the plot window
        
        if (upper == "points") {
          plot_density_contours(x, y, xlim, ylim, 
                    plot.dots=TRUE, dot.level=0, plot.contours=FALSE, plot.image=FALSE,
                    col=col, cex=cex, ...)
        }
        if (upper == "image") {
          plot_density_contours(x, y, xlim, ylim, npix=npix,
                    plot.dots=FALSE, dot.level=0, plot.contours=FALSE, plot.image=TRUE,
                    col=col, cex=cex, ...)
        }
        if (upper == "contour") {
          plot_density_contours(x, y, xlim, ylim, prob.levels, npix=npix,
                    plot.dots=FALSE, dot.level=0, plot.contours=TRUE, plot.image=FALSE,
                    col=col, cex=cex, ...)
        }
        
       # draw a box around the plot region
        
        box(which="plot")
      }            # end if (i > j)

    }              # end loop over j
  }                # end loop over i
  
  # output the 1D intervals
  
  theta.intervals <- NULL
  if (!is.null(levels1d)) {
    nlevels <- length(prob1d)
    rownames(levels1d) <- paste(labels)
    colnames(levels1d)[nlevels+1] <- "mean"
    theta.intervals <- signif(levels1d,4)
  }
  
  # restore graphical parameters
  par(retire)
  
  # return the intervals
  return(theta.intervals)

}
