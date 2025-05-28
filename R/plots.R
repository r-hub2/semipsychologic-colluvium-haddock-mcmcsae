#' Trace, density and autocorrelation plots for (parameters of a) draws
#' component (dc) object
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=50)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' plot(sim$u)
#' }
#'
#' @export
## @method plot dc
#' @param x a draws component object.
#' @param nrows number of rows in plot layout.
#' @param ncols number of columns in plot layout.
#' @param ask ask before plotting the next page; default is \code{FALSE}.
#' @param ... arguments passed to \code{\link[stats]{density}}.
plot.dc <- function(x, nrows, ncols, ask=FALSE, ...) {
  if (missing(nrows)) nrows <- min(5L, n_vars(x))
  if (missing(ncols)) ncols <- 3L * min(2L, ceiling(n_vars(x) / nrows))
  oldpar <- par(no.readonly=TRUE)
  on.exit({
    safe_params <- intersect(names(oldpar), names(par()))
    safe_params <- setdiff(safe_params, c("pin", "plt", "usr", "fig"))
    par(oldpar[safe_params])
  })
  par(mar=c(2.4,2.8,0.9,0.8))
  par(mfrow=c(nrows, ncols))
  par(ask=ask)
  labs <- labels(x)
  for (i in seq_len(n_vars(x))) {
    drawsi <- get_from(x, i)
    drawsi <- do.call("cbind", drawsi)  # niter x nchain matrix
    lab <- if (is.null(labs)) "" else labs[i]
    trace_plot(drawsi, ylab=lab)
    drawsi <- as.vector(drawsi)  # concatenate the chains
    plot(density(drawsi, ...), main=lab)
    acf(drawsi, main="")  # note that the concatenation of chains disturbs the acf a little, especially if the chains are short
  }
}

#' Trace, density and autocorrelation plots
#' 
#' Trace, density and autocorrelation plots for selected components of an \code{mcdraws} object.
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=50)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' plot(sim, c("beta", "u", "u_sigma", "v_sigma"), ask=TRUE)
#' }
#'
#' @export
## @method plot mcdraws
#' @param x an object of class \code{mcdraws}.
#' @param vnames optional character vector to select a subset of parameters.
#' @param nrows number of rows in plot layout.
#' @param ncols number of columns in plot layout.
#' @param ask ask before plotting the next page; default is \code{FALSE}.
#' @param ... arguments passed to \code{\link[stats]{density}}.
plot.mcdraws <- function(x, vnames, nrows, ncols, ask=FALSE, ...) {
  obj <- NULL
  labs <- NULL
  for (v in vnames) {
    vars <- get_name_range(v, FALSE)
    if (all(vars$name != par_names(x))) {
      warn("parameter '", vars$name, "' not in MCMC output")
      next
    }
    if (is.null(vars$range)) vars$range <- seq_len(n_vars(x[[vars$name]]))
    if (is.null(obj)) {
      obj <- get_from(x[[vars$name]], vars$range)
    } else {
      for (i in seq_along(obj))
        obj[[i]] <- cbind(obj[[i]], as.matrix.dc(get_from(x[[vars$name]], vars$range, i), colnames=FALSE))
    }
    labs <- c(labs, paste0(vars$name, labels(x[[vars$name]])[vars$range]))
  }
  if (!is.null(obj)) {
    class(obj) <- "dc"
    attr(obj, "labels") <- labs
    plot(obj, nrows, ncols, ask, ...)  # plot dc object
  }
}

#' Plot trace of MCMC output for a single (univariate) parameter
#'
#' @noRd
#' @param dc1 a matrix with MCMC output where different chains are stored in different columns.
#' @param xlab x-axis label.
#' @param ylab y-axis label.
trace_plot <- function(dc1, xlab="iterations", ylab="") {
  if (!is.matrix(dc1)) dc1 <- do.call("cbind", dc1)  # niter x nchain matrix
  matplot(seq_len(nrow(dc1)), dc1, xlab=xlab, ylab=ylab, type="l", col=seq_len(ncol(dc1)))
}

#' Plot a set of model coefficients or predictions with uncertainty intervals
#' based on summaries of simulation results or other objects.
#'
#' This function plots estimates with error bars. Multiple sets of
#' estimates can be compared. The error bars can either be based on
#' standard errors or on explicitly specified lower and upper bounds.
#' The function is adapted from function \code{plot.sae} in package
#' \pkg{hbsae}, which in turn was adapted from function
#' \code{coefplot.default} from package \pkg{arm}.
#'
#' @examples
#' \donttest{
#' # create artificial data
#' set.seed(21)
#' n <- 100
#' dat <- data.frame(
#'   x=runif(n),
#'   f=factor(sample(1:20, n, replace=TRUE))
#' )
#' model <- ~ reg(~ x, prior=pr_normal(precision=1), name="beta") + gen(factor=~f, name="v")
#' gd <- generate_data(model, data=dat)
#' dat$y <- gd$y
#' # fit a base model
#' model0 <- y ~ reg(~ 1, name="beta") + gen(factor=~f, name="v")
#' sampler <- create_sampler(model0, data=dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' (summ0 <- summary(sim))
#' # fit 'true' model
#' model <- y ~ reg(~ x, name="beta") + gen(factor=~f, name="v")
#' sampler <- create_sampler(model, data=dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' (summ <- summary(sim))
#' # compare random effect estimates against true parameter values
#' plot_coef(summ0$v, summ$v, list(est=gd$pars$v), n.se=2, offset=0.2,
#'   maxrows=10, est.names=c("base model", "true model", "true"))
#' }
#'
#' @export
#' @param ... \code{dc_summary} objects (output by the \code{summary} method for
#'  simulation objects of class \code{dc}), \code{sae} objects (output by the
#'  functions of package \pkg{hbsae}), or lists. In case of a list the components
#'  used are those with name \code{est} for point estimates, \code{se}
#'  for standard error based intervals or \code{lower} and \code{upper} for
#'  custom intervals. Instead of \code{dc_summary} objects matrix objects are
#'  also supported as long as they contain columns named "Mean" and "SD" as do
#'  \code{dc_summary} objects. Named parameters of other types that do not match any
#'  other argument names are passed to lower-level plot functions.
#' @param n.se number of standard errors below and above the point estimates
#'  to use for error bars. By default equal to 1. This only refers to the
#'  objects of class \code{dc_summary} and \code{sae}.
#' @param est.names labels to use in the legend for the components of the \code{...} argument
#' @param sort.by vector by which to sort the coefficients, referring to the first object passed.
#' @param decreasing if \code{TRUE}, sort in decreasing order (default).
#' @param index vector of names or indices of the selected areas to be plotted.
#' @param maxrows maximum number of rows in a column.
#' @param maxcols maximum number of columns of estimates on a page.
#' @param offset space used between plots of multiple estimates for the same area.
#' @param cex.var the font size for the variable names, default=0.8.
#' @param mar a numeric vector of the form \code{c(bottom, left, top, right)},
#'  specifying the number of lines of margin on each of the four sides of the plot.
plot_coef <- function(..., n.se=1, est.names, sort.by=NULL, decreasing=FALSE,
                      index=NULL, maxrows=50L, maxcols=6L,
                      offset=0.1, cex.var=0.8, mar=c(0.1,2.1,5.1,0.1)) {

  dotargs <- list(...)
  toplot <- b_apply(dotargs, \(obj) any(class(obj)[1L] == c("dc_summary", "matrix", "sae", "list")))
  grpar <- dotargs[!toplot]  # other graphical parameters
  if (!length(grpar)) grpar <- NULL

  x <- dotargs[toplot]
  if (!length(x)) stop("nothing to plot")

  # extract point estimates to plot
  ests <- lapply(x, \(obj) {
    switch(class(obj)[1L],
      dc_summary=, matrix = obj[, "Mean"],
      sae=, list = obj$est,
      stop("unsupported class '", class(obj)[1L], "'")
    )
  })

  ints <- lapply(x, \(obj) {
    # for dc_summary and sae objects by default use mean +/- se intervals
    switch(class(obj)[1L],
      dc_summary=, matrix = cbind(obj[, "Mean"] - n.se * obj[, "SD"], obj[, "Mean"] + n.se * obj[, "SD"]),
      sae = cbind(obj$est - n.se * sqrt(obj$mse), obj$est + n.se * sqrt(obj$mse)),
      list = {
        if (!is.null(obj$se)) {
          if (is.null(obj$est)) stop("list with 'se' but no 'est' component")
          cbind(obj$est - n.se * obj$se, obj$est + n.se * obj$se)
        } else if (is.null(obj$lower)) {
          if (!is.null(obj$upper)) stop("'upper' specified without 'lower'")
          if (is.null(obj$est)) stop("object without 'est' and 'lower'")
          NULL
        } else if (is.null(obj$upper)) {
          stop("'lower' specified without 'upper'")
        } else {
          cbind(obj$lower, obj$upper)
        }
      }
    )
  })

  M <- length(ests[[1L]])
  lab <- names(ests[[1L]])
  if (is.null(lab)) lab <- as.character(seq_len(M))

  # make all objects compatible with the first, i.e. of the same length, if possible
  if (length(ests) > 1L) {
    for (i in 2:length(ests)) {
      obj.est <- ests[[i]]
      obj.int <- ints[[i]]
      if (is.null(names(obj.est))) {
        if (length(obj.est) != M) stop("unable to match components of input objects")
        names(ests[[i]]) <- lab
        if (!is.null(ints[[i]])) {
          if (nrow(obj.int) != M) stop("unable to match components of input objects")
          rownames(ints[[i]]) <- lab
        }
      } else {
        ind <- fmatch(names(obj.est), lab)
        ests[[i]] <- NA_real_ * ests[[1L]]
        ests[[i]][ind[!is.na(ind)]] <- obj.est[!is.na(ind)]
        if (!is.null(ints[[i]])) {
          ints[[i]] <- NA_real_ * ints[[1L]]
          ints[[i]][ind[!is.na(ind)], ] <- obj.int[!is.na(ind), ]
        }
      }
    }
  }

  # allow plotting for a subset only
  if (is.null(index)) {
    o <- lab
  } else {
    if (is.character(index)) {
      if (anyNA(fmatch(index, lab))) stop("some elements of index cannot be matched")
      o <- index
    } else {
      o <- lab[index]
    }
    M <- length(o)
  }

  if (is.null(sort.by)) {
    sort.by <- seq_along(ests[[1L]])
  } else {
    if (length(sort.by) != M) stop("'sort.by' must have same length as first plot object")
  }
  names(sort.by) <- lab

  o <- o[order(sort.by[o], decreasing=decreasing)]

  maxrows <- min(maxrows, M)
  cols <- ceiling(M/maxrows)
  pages <- ceiling(cols/maxcols)

  compute.xlim <- all("xlim" != names(grpar))

  oldpar <- par(no.readonly=TRUE)
  on.exit({
    safe_params <- intersect(names(oldpar), names(par()))
    safe_params <- setdiff(safe_params, c("pin", "plt", "usr", "fig"))
    par(oldpar[safe_params])
  })
  for (page in seq_len(pages)) {
    colrange <- ((page - 1L) * maxcols + 1L):min((page * maxcols), cols)
    plot.new()
    par(mfrow=c(1L, min(maxcols, cols)))
    for (co in colrange) {
      par(mfg=c(1L, co - (colrange[1L] - 1L)))  # force plot in the intended column
      rowrange <- ((co - 1L) * maxrows + 1L):min((co * maxrows), M)
      offsets <- offset * (seq_along(x) - 1L)  # - ((length(x) + 1) %/% 2))
      # first determine xlim, if not set manually
      if (compute.xlim) {
        xlim <- range(ests[[1L]][o[rowrange]][1L], na.rm=TRUE)
        for (i in seq_along(x)) {
          if (is.null(ints[[i]]))
            xlim <- range(xlim, ests[[i]][o[rowrange]], na.rm=TRUE)
          else
            xlim <- range(xlim, ints[[i]][o[rowrange], ], na.rm=TRUE)
        }
        grpar$xlim <- xlim
      }
      for (i in seq_along(x)) {
        if (is.null(ints[[i]]))
          interv <- NULL
        else
          interv <- ints[[i]][o[rowrange], ]
        if (i == 1L) {
          do.call(cplot, c(list(ests[[i]][o[rowrange]], interv, varnames=o[rowrange],
            cex.var=cex.var, offset=offsets[i], mar=mar), grpar))
        } else {
          cplot(ests[[i]][o[rowrange]], interv, col.pts=i, add=TRUE, offset=offsets[i], mar=mar)
        }
      }
    }  # END for co
    if (!missing(est.names) || length(ests) > 1L) {  # add a legend
      par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
      plot(0, 0, type="n", bty="n", xaxt="n", yaxt="n")
      if (missing(est.names))
        est.names <- paste0("est", seq_along(ests))
      dots.only <- b_apply(x, \(obj) is.list(obj) && is.null(obj[["se"]]) && is.null(obj[["lower"]]))
      legend("topright", est.names, xpd=TRUE, horiz=TRUE, inset=c(0, 0),
        bty="n", pch=rep.int(20L, length(ests)), lty=ifelse(dots.only, 0, 1),
        col=seq_along(ests)
      )
    }
  }  # END for page
}

# Plot confidence or credible intervals for vectors of coefficients
# 
# This function is adapted from coefplot.default of package \pkg{arm}.
cplot <- function (coefs, intervals=NULL,
                   varnames=NULL, v.axis=TRUE, h.axis=TRUE,
                   cex.var=0.8, cex.pts=0.9, col.pts=1, pch.pts=20, 
                   var.las=2, xlab="", ylab="", main="",
                   mar=c(0.1,2.1,5.1,0.1), plot=TRUE,
                   add=FALSE, offset=0.1, ...) {
  
  if (is.list(coefs)) coefs <- unlst(coefs)
  m <- length(coefs)
  id <- seq_len(m)
  if (is.null(intervals)) intervals <- cbind(coefs, coefs)
  if (is.null(varnames))
    maxchar <- 0
  else
    maxchar <- max(i_apply(varnames, nchar))
  k <- 1/m
  oldmar <- par("mar")
  on.exit(par(mar=oldmar))
  mar[2L] <- max(oldmar[2L], trunc(mar[2L] + maxchar/10)) + 0.1
  par(mar=mar)
  if (plot) {
    if (add) {
      id <- id + offset
      points(coefs, id, pch=pch.pts, cex=cex.pts, col=col.pts)
      segments(intervals[, 1L], id, intervals[, 2L], id, lwd=1, col=col.pts)
    } else {
      plot(intervals, c(id + k, id - k), type="n",
        axes=FALSE, main=main, xlab=xlab, ylab=ylab, ...)
      if (h.axis) axis(3L)
      if (v.axis) {
        axis(2L, rev(id), varnames[rev(id)], las=var.las,
          tck=FALSE, lty=0, cex.axis=cex.var)
      }
      abline(v=0, lty=2L)
      points(coefs, id, pch=pch.pts, cex=cex.pts, col=col.pts)
      segments(intervals[, 1L], id, intervals[, 2L], id, lwd=1, col=col.pts)
    }
  } else {  # do not plot
    plot(intervals, c(id + k, id - k), type="n",
      axes=FALSE, main="", xlab=xlab, ylab=ylab, ...)
  }
}
