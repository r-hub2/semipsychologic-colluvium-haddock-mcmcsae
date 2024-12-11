#' Run a Markov Chain Monte Carlo simulation
#'
#' Given a sampler object this function runs a MCMC simulation and stores the
#' posterior draws. A sampler object for a wide class of multilevel models
#' can be created using \code{\link{create_sampler}}, but users can also define
#' their own sampler functions, see below.
#' \code{MCMCsim} allows to choose the parameters for which simulation results
#' must be stored. It is possible to define derived quantities that will also
#' be stored. To save memory, it is also possible to only store Monte Carlo
#' means/standard errors for some large vector parameters, say. Another
#' way to use less memory is to save the simulation results of large vector
#' parameters to file.
#' For parameters specified in \code{plot.trace} trace plots or pair plots of
#' multiple parameters are displayed during the simulation.
#'
#' A sampler object is an environment containing data and functions to use
#' for sampling. The following elements of the sampler object are used by
#' \code{MCMCsim}:
#' \describe{
#'   \item{start}{function to generate starting values.}
#'   \item{draw}{function to draw samples, typically from a full conditional
#'     posterior distribution.}
#'   \item{rprior}{function to draw from a prior distribution.}
#'   \item{coef.names}{list of vectors of parameter coefficient names, for
#'     vector parameters.}
#   next two elements undocumented; used by create_sampler, but awkward
#   \item{store_default}{function that returns a character vector of parameter
#     names that should be stored by default.}
#   \item{store_mean_default}{function that returns a character vector of
#     parameter names whose means should be stored by default.}
#'   \item{MHpars}{vector of names of parameters that are sampled using a
#'     Metropolis-Hastings (MH) sampler; acceptance rates are kept for these
#'     parameters.}
#'   \item{adapt}{function of acceptance rates of \code{MHpars} to adapt
#'     MH-kernel, called every 100 iterations during the burn-in period.}
#' }
#'
#' @examples
#' # 1. create a sampler function
#' sampler <- new.env()
#' sampler$draw <- function(p) list(x=rnorm(1L), y=runif(1L))
#' # 2. do the simulation
#' sim <- MCMCsim(sampler, store=c("x", "y"))
#' str(sim)
#' summary(sim)
#'
#' # example that requires start values or a start function
#' sampler$draw <- function(p) list(x=rnorm(1L), y=p$x * runif(1L))
#' sampler$start <- function(p) list(x=rnorm(1L), y=runif(1L))
#' sim <- MCMCsim(sampler, store=c("x", "y"))
#' summary(sim)
#' plot(sim, c("x", "y"))
#'
#' # example using create_sampler; first generate some data
#' n <- 100
#' dat <- data.frame(x=runif(n), f=as.factor(sample(1:4, n, replace=TRUE)))
#' gd <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1), name="beta"), data=dat)
#' dat$y <- gd$y
#' sampler <- create_sampler(y ~ x + f, data=dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400, n.chain=2)
#' (summary(sim))
#' gd$pars
#'
#' @export
#' @param sampler sampler object created by \code{\link{create_sampler}}.
#' @param from.prior whether to sample from the prior. By default \code{from.prior=FALSE}
#'  and samples are taken from the posterior.
#' @param n.iter number of draws after burnin.
#' @param n.chain number of independent chains.
#' @param thin only every \code{thin}'th draw is kept.
#' @param burnin number of draws to discard at the beginning of each chain.
#' @param start an optional function to generate starting values or a list containing for each chain
#'  a named list of starting values. It may be used to provide starting values for some or all parameters.
#'  The sampler object's own start function, if it exists, is called to generate any starting values not
#'  provided by the user.
#' @param store vector of names of parameters to store MCMC draws for. By default, simulations are
#'  stored for all parameters returned by \code{sampler$store_default}.
#' @param store.all if \code{TRUE} simulation vectors of all parameters returned by the sampling
#'  function of \code{sampler} will be stored. The default is \code{FALSE}, and in that case
#'  only simulations for the parameters named in \code{store} are stored.
#' @param pred list of character strings defining derived quantities to be computed (and stored) for each draw.
#' @param store.mean vector of names of parameters for which only the mean (per chain) is to be stored.
#'  This may be useful for large vector parameters (e.g. regression residuals) for which storing complete
#'  MCMC output would use too much memory. The function \code{sampler$store_mean_default}
#'  exists it provides the default.
# By default this is used only for residuals. If the data-level precision is modeled, then
# the means of the square roots of the precisions are stored as well. These are required for the
# computation of DIC. If \code{compute.weights=TRUE} means of the weights are stored.
#' @param store.sds if \code{TRUE} store for all parameters in \code{store.mean}, besides the mean, also
#'  the standard deviation. Default is \code{FALSE}.
#' @param to.file vector of names of parameters to write to file.
#' @param filename name of file to write parameter draws to.
#'  Each named parameter is written to a separate file, named \code{filename_parametername}.
#' @param write.single.prec Whether to write to file in single precision. Default is \code{FALSE}.
#' @param verbose if \code{FALSE} no output is sent to the screen during the simulation. \code{TRUE} by default.
#' @param n.progress update diagnostics and plots after so many iterations.
#' @param trace.convergence vector of names of parameters for which Gelman-Rubin R-hat diagnostics are printed to the screen every \code{n.progress} iterations.
#' @param stop.on.convergence if \code{TRUE} stop the simulation if the R-hat diagnostics for all parameters in \code{trace.convergence} are less than \code{convergence.bound}.
#' @param convergence.bound threshold used with \code{stop.on.convergence}.
#' @param plot.trace character vector of parameter names for which to plot draws
#'  during the simulation. For one or two parameters trace plots will be shown,
#'  and if more parameters are specified the results will be displayed in a pairs
#'  plot. For vector parameters a specific component can be selected using brackets,
#'  e.g. \code{"beta[2]"}.
#' @param add.to.plot if \code{TRUE} the plot is updated every \code{n.progress} iterations,
#'  otherwise a new plot (with new scales) is created after every \code{n.progress} iterations.
#' @param plot.type default is "l" (lines).
#' @param n.cores the number of cpu cores to use. Default is 1, i.e. no parallel computation.
#'  If an existing cluster \code{cl} is provided, \code{n.cores} will be set to the number
#'  of workers in that cluster.
#' @param cl an existing cluster can be passed for parallel computation. If \code{NULL} and
#'  \code{n.cores > 1}, a new cluster is created.
#' @param seed a random seed (integer). For parallel computation it is used to independently
#'  seed RNG streams for all workers.
#' @param export a character vector with names of objects to export to the workers. This may
#'  be needed for parallel execution if expressions in \code{pred} depend on global variables.
#' @returns An object of class \code{mcdraws} containing posterior draws as well as some meta information.
MCMCsim <- function(sampler, from.prior=FALSE, n.iter=1000L, n.chain=3L, thin=1L,
                    burnin=if (from.prior) 0L else 250L, start=NULL,
                    store, store.all=FALSE, pred=NULL,
                    store.mean, store.sds=FALSE,
                    to.file=NULL, filename="MCdraws_", write.single.prec=FALSE,
                    verbose=TRUE, n.progress=n.iter %/% 10L,
                    trace.convergence=NULL, stop.on.convergence=FALSE, convergence.bound=1.05,
                    plot.trace=NULL, add.to.plot=TRUE, plot.type="l",
                    n.cores=1L, cl=NULL, seed=NULL, export=NULL) {

  started <- proc.time()

  # checks
  n.iter <- as.integer(round(n.iter))
  if (n.iter < 0L) stop("'n.iter' must be nonnegative")
  burnin <- as.integer(round(burnin))
  if (burnin < 0L) stop("'burnin' must be nonnegative")
  n.chain <- as.integer(round(n.chain))
  if (n.chain <= 0L) stop("'n.chain' must be positive")
  thin <- as.integer(round(thin))
  if (thin < 1L) stop("'thin' must be at least 1 (= no thinning)")
  n.progress <- max(as.integer(round(n.progress)), thin)
  if (missing(sampler)) stop("argument 'sampler' cannot be empty")
  if (n.chain < 1L) stop("'n.chain' must be at least 1")
  chains <- seq_len(n.chain)

  if (is.null(cl))
    n.cores <- min(as.integer(n.cores)[1L], n.chain)
  else
    n.cores <- length(cl)
  if (n.cores > 1L) {
    if (n.cores > parallel::detectCores()) stop("too many cores")
    if (is.null(cl)) {
      cl <- setup_cluster(n.cores, seed)
    } else {
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    }
    if (!is.null(export)) parallel::clusterExport(cl, export, parent.frame())
    chains.per.worker <- rep(n.chain %/% n.cores, n.cores) + rep(1:0, c(n.chain %% n.cores, n.cores - n.chain %% n.cores))
    distr.list <- vector("list", n.cores)
    for (i in seq_len(n.cores)) distr.list[[i]]$n.chain <- chains.per.worker[i]
    if (is.null(start) || is.function(start)) {
      start.list <- rep(list(start), n.cores)
    } else {
      if (!is.list(start) || (length(start) != n.chain) || !all(b_apply(start, is.list)))
        stop("start has wrong format; see the help for possible ways to specify initial values")
      start.list <- split(start, rep(seq_len(n.cores), chains.per.worker))
    }
    for (i in seq_len(n.cores)) distr.list[[i]]$start <- start.list[[i]]
    message("distribution of chains over workers: ", paste0(chains.per.worker, collapse=", "))
    parallel::clusterExport(cl, "sampler", envir=environment())
    mc <- match.call()
    mc$sampler <- quote(sampler)
    mc$n.chain <- quote(worker$n.chain)
    mc$start <- quote(worker$start)
    mc$verbose <- FALSE
    mc$n.cores <- 1L
    mc$cl <- mc$seed <- mc$export <- NULL
    run_chains <- function(worker, ...) eval(mc)
    message("running the simulation")
    # TODO allow to write to file using separate files per worker, to combine afterwards
    out <- do.call(combine_chains, parallel::parLapply(cl, distr.list, run_chains))
    out[["_cluster"]] <- cl  # store the cluster for parallel post-processing
    if (verbose) {
      cat("\n")
      print(proc.time() - started)
    }
    return(out)
  } else {
    if (!is.null(seed)) set.seed(seed)
  }

  # derived quantities
  if (is.null(pred)) {
    do.pred <- FALSE
  } else {
    if (!is.list(pred) || is.null(names(pred))) stop("'pred' must be a named list")
    predict <- function(p) {}
    for (j in seq_along(pred))
      predict <- add(predict, substitute(p[[comp]] <- evalq(expr, env=p), list(expr=str2lang(pred[[j]]), comp=names(pred)[j])))
    predict <- add(predict, quote(p))
    do.pred <- TRUE
  }

  # set up starting values
  p <- list()
  if (!is.null(start) && !is.function(start))
    if (!is.list(start) || (length(start) != n.chain) || !all(b_apply(start, is.list)))
      stop("start has wrong format; see the help for possible ways to specify initial values")
  for (ch in chains) {
    # user-defined function or values can be provided for some or all of the parameters
    if (is.null(start)) {
      p[[ch]] <- list()
    } else {
      if (is.function(start)) {
        p[[ch]] <- start()
      } else {  # start is a list
        p[[ch]] <- start[[ch]]
      }
    }
    # if sampler provides a start function, use it to generate start values
    #   for parameters for which no values have been provided by the user
    if (!is.null(sampler$start)) {
      # NB currently we rely on sampler$start to not overwrite existing (user-provided) start values
      #    might also skip assignment of user-provided parameters here (TODO)
      p[[ch]] <- sampler$start(p[[ch]])
    }
  }

  if (from.prior)
    draw <- sampler$rprior  # function that returns a draw from the prior
  else {
    draw <- sampler$draw  # function that returns a draw from the posterior
    if (is.null(draw)) {
      if (sampler$prior.only)
        stop("sampler does not have a function 'draw'; please use 'from.prior=TRUE' if you want to sample from priors")
      else
        stop("sampler does not have a function 'draw'")
    }
  }

  # do a test draw both for timing and storage information
  timing <- system.time({
    test_draw <- draw(p[[1L]])
  })
  if (do.pred) {
    if (any(names(pred) %in% names(test_draw)))
      stop("names in 'pred' clashing with sampler's parameter names: ",
           paste0(intersect(names(pred), names(test_draw)), collapse=", "))
    timing <- timing + system.time(test_draw <- predict(test_draw))
  }
  if (verbose && (burnin + n.iter) * n.chain * timing["elapsed"] > 1)
    cat("\nEstimated computation time:", format((burnin + n.iter) * n.chain * timing["elapsed"]), "s\n")

  # derive names and ranges of variables to trace
  trace.convergence <- lapply(trace.convergence, get_name_range)
  plot.trace <- lapply(plot.trace, get_name_range)
  if (any(b_apply(plot.trace, function(x) length(x$range) > 1L))) stop("for plots only single elements of vector parameters may be specified")
  
  # remove components to trace or plot that are not in the sampler's output
  trace.convergence <- Filter(function(x) any(x$name == names(test_draw)), trace.convergence)
  plot.trace <- Filter(function(x) any(x$name == names(test_draw)), plot.trace)
  if (length(plot.trace)) {
    oldpar <- par(no.readonly=TRUE)
    on.exit(par(oldpar))
    par(mar=c(4, 3.8, 0.9, 0.8), cex.axis=0.8, cex.lab=0.8)
  }
  # check that ranges, if specified, are allowed
  for (v in trace.convergence)
    if (!all(v$range %in% seq_along(test_draw[[v$name]]))) stop("illegal range in trace.convergence for parameter ", v$name)
  for (v in plot.trace)
    if (!all(v$range %in% seq_along(test_draw[[v$name]]))) stop("illegal range in plot.trace for parameter ", v$name)
  trace.convergence.names <- unlst(lapply(trace.convergence, function(v) if (length(test_draw[[v$name]]) == 1L) v$name else paste0(v$name, "[", v$range, "]")))
  plot.trace.names <- unlst(lapply(plot.trace, function(v) if (length(test_draw[[v$name]]) == 1L) v$name else paste0(v$name, "[", v$range, "]")))
  if (missing(store.mean)) {  # set parameters whose mean is stored by default
    if (is.null(sampler$store_mean_default))
      store.mean <- NULL
    else
      store.mean <- sampler$store_mean_default(from.prior)
  }

  if (missing(store)) {
    if (is.null(sampler$store_default))
      store <- NULL
    else
      store <- sampler$store_default(from.prior)
    if (store.all)
      store <- union(store, Filter(function(x) !endsWith(x, "_") && all(x != store.mean), names(test_draw)))
  }
  # always store derived quantities and 'diagnostic' parameters
  if (do.pred) store <- union(store, names(pred))
  if (length(trace.convergence)) store <- union(store, s_apply(trace.convergence, `[[`, "name"))
  if (length(plot.trace)) store <- union(store, s_apply(plot.trace, `[[`, "name"))

  out <- list()

  out[["_means"]] <- list()

  for (v in store.mean) {
    if (is.null(test_draw[[v]])) {
      warn("parameter '", v, "' not in sampler output")
      store.mean <- setdiff(store.mean, v)
      next
    }
    out[["_means"]][[v]] <- list()
    if (store.sds) out[["_sds"]][[v]] <- list()
    for (ch in chains) {
      out[["_means"]][[v]][[ch]] <- 0 * test_draw[[v]]
      if (store.sds) {
        out[["_sds"]][[v]][[ch]] <- 0 * test_draw[[v]]
      }
    }
  }

  # set up vectors and matrices to store MCMC draws
  n.draw <- n.iter %/% thin  # number of draws retained
  list.store <- NULL  # only used in case non-vector objects should be stored
  for (v in store) {
    if (is.null(test_draw[[v]])) {
      warn("parameter '", v, "' not in sampler output")
      store <- setdiff(store, v)
      next
    }
    if (!is.numeric(test_draw[[v]])) {
      if (is.list(test_draw[[v]]) && substring(v, nchar(v) - 6L) == "_trees_") {
        # allow storing non-vector objects in a list, e.g. BART trees
        out[[v]] <- list()
        for (ch in chains)
          out[[v]][[ch]] <- vector(mode="list", length=n.draw)
        list.store <- c(list.store, v)
      } else {
        warn("non-numeric parameter '", v, "'")
      }
      store <- setdiff(store, v)
      next
    }
    if (!length(test_draw[[v]])) {
      warn("parameter '", v, "' has length 0")
      store <- setdiff(store, v)
      next
    }
    out[[v]] <- list()
    for (ch in chains)
      out[[v]][[ch]] <- matrix(NA_real_, nrow=n.draw, ncol=length(test_draw[[v]]))
  }

  out[["_info"]] <- list(n.iter=n.iter, n.chain=n.chain, thin=thin, burnin=burnin, n.draw=n.draw,
                         from.prior=from.prior, parnames=store, list.pars=list.store,
                         call=match.call())
  out[["_state"]] <- p  # current state of the system
  out[["_model"]] <- sampler  # pointer to the sampler's scope

  if (verbose) cat("\nOutput size:", format(c(object.size(out))/1e6, digits=3L), "MB\n\n")
  if (is.null(to.file)) {
    write.to.file <- FALSE
  } else {
    outfile <- list()
    for (v in to.file) {  # one output file per parameter (multiple chains in the same file)
      if (is.null(test_draw[[v]])) {
        warn("parameter '", v, "' not in sampler output")
        to.file <- setdiff(to.file, v)
      } else if (!is.numeric(test_draw[[v]])) {
        warn("non-numeric parameter '", v, "'")
        to.file <- setdiff(to.file, v)
      } else if (!length(test_draw[[v]])) {
        warn("parameter '", v, "' has length 0")
        to.file <- setdiff(to.file, v)
      } else {
        outfile[[v]] <- file(paste0(filename, v, ".dat"), "wb")
        write_header(outfile[[v]], n.iter, n.chain, length(test_draw[[v]]),
          if (is.null(sampler$coef.names)) NULL else sampler$coef.names[[v]], write.single.prec)
      }
    }
    write.size <- if (write.single.prec) 4L else NA_integer_
    write.to.file <- length(to.file) >= 1L
  }

  rm(test_draw, timing)

  # accept/reject vectors for MH parameters
  # MH parameters are assumed to be scalar parameters
  MH <- length(sampler$MHpars) >= 1L
  if (MH) {
    MHpars <- sampler$MHpars
    acc.rates <- rep.int(list(0), length(MHpars))
    names(acc.rates) <- MHpars
    out[["_accept"]] <- rep.int(list(rep.int(list(0L), n.chain)), length(MHpars))
    names(out[["_accept"]]) <- MHpars
    previous.values <- rep.int(list(rep.int(list(0), n.chain)), length(MHpars))
    adapt <- is.function(sampler$adapt)
  } else {
    adapt <- FALSE
  }

  # change to points if there is only one point to plot in each update
  if (n.progress %/% thin == 1L) plot.type <- "p"

  # do the simulation
  # part 1: burnin, possibly adaptation for MH
  for (i in seq_len(burnin)) {
    for (ch in chains) p[[ch]] <- draw(p[[ch]])

    if (adapt) {  # adaptive MH
      for (a in seq_along(MHpars))
        for (ch in chains) {
          if (any(previous.values[[a]][[ch]] != p[[ch]][[MHpars[a]]]))
            out[["_accept"]][[a]][[ch]] <- out[["_accept"]][[a]][[ch]] + 1L
          previous.values[[a]][[ch]] <- p[[ch]][[MHpars[a]]]
        }

      # adaptation of MH proposals
      if (i %% 50L == 0L) {  # only adapt during burnin
        # adapt for each chain in the same way, depending on average acceptance rates
        for (a in seq_along(MHpars)) {
          acc.rates[[a]] <- Reduce("+", out[["_accept"]][[a]]) / (n.chain * 50L)
          out[["_accept"]][[a]] <- rep.int(list(0L), n.chain)  # reset
        }
        sampler$adapt(acc.rates)  # adapt MH proposals
      }
    }

    if (i %% 10L == 0L && verbose) cat("\rburnin iteration", i)
  }  # END for (i in seq_len(burnin))
  if (verbose) cat("\r                         ")  # blank the output

  if (adapt) {  # reset, only keep acceptance rates after burnin
    for (a in seq_along(MHpars))
      out[["_accept"]][[a]] <- rep.int(list(0L), n.chain)
  }

  # part 2: sampling after burnin, store parameters, compute predicted quantities
  for (i in seq_len(n.iter)) {
    for (ch in chains) p[[ch]] <- draw(p[[ch]])
    if (i %% thin == 0L) {  # store
      index <- i %/% thin  # storage index
      if (do.pred) for (ch in chains) p[[ch]] <- predict(p[[ch]])
      if (!is.null(list.store)) for (ch in chains) for (v in list.store) out[[v]][[ch]][[index]] <- p[[ch]][[v]]
      for (ch in chains) {
        for (v in store) out[[v]][[ch]][index, ] <- p[[ch]][[v]]
        for (v in store.mean) {
          v_update(out[["_means"]][[v]][[ch]], plus=TRUE, p[[ch]][[v]])
          if (store.sds)
            v_update(out[["_sds"]][[v]][[ch]], plus=TRUE, p[[ch]][[v]]^2)
        }
      }  # END for (ch in chains)
      if (write.to.file)
        for (v in to.file)
          writeBin(unlst(lapply(p, `[[`, v)), con=outfile[[v]], size=write.size)
    }  # END if (i %% thin == 0L)

    if (MH)
      for (a in seq_along(MHpars))
        for (ch in chains) {
          if (any(previous.values[[a]][[ch]] != p[[ch]][[MHpars[a]]]))
            out[["_accept"]][[a]][[ch]] <- out[["_accept"]][[a]][[ch]] + 1L
          previous.values[[a]][[ch]] <- p[[ch]][[MHpars[a]]]
        }

    if (i %% 10L == 0L && verbose) cat("\riteration", i)
    if (i %% n.progress == 0L) {  # show progress
      if (length(trace.convergence) && (n.chain > 1L) && (index > 1L)) {
        diagnostic <- unlst(lapply(trace.convergence, function(v) R_hat(get_from(out[[v$name]], v$range, draws=seq_len(index)))))
        names(diagnostic) <- trace.convergence.names
        if (verbose) {
          cat("\n")
          print(diagnostic, digits=3L)
        }
        if (stop.on.convergence && max(diagnostic) < convergence.bound) {
          # TODO remove unused part of data objects
          break
        }
      }
      if (verbose && length(plot.trace)) {
        i0 <- max(index - (n.progress %/% thin) + 1L, 1L)
        if (length(plot.trace) == 1L) {  # 1d traceplot
          yvalues <- get_from(out[[plot.trace[[1L]]$name]], vars=plot.trace[[1L]]$range, draws=i0:index)
          xvalues <- thin * (i0:index)
          for (ch in chains) {
            if (add.to.plot) {
              if (i0 == 1L && ch == 1L)
                plot(xvalues, yvalues[[ch]], type=plot.type, xlim=thin*c(1L, n.draw), ylim=range(unlst(yvalues)), xlab="iteration", ylab=plot.trace.names)
              else
                lines(xvalues, yvalues[[ch]], type=plot.type, col=ch)
            } else {
              if (ch == 1L)
                plot(xvalues, yvalues[[ch]], type=plot.type, ylim=range(unlst(yvalues)), xlab="iteration", ylab=plot.trace.names)
              else
                lines(xvalues, yvalues[[ch]], type=plot.type, col=ch)
            }
          }  # END for (ch in chains)
        } else if (length(plot.trace) == 2L) {  # 2d traceplot
          xvalues <- get_from(out[[plot.trace[[1L]]$name]], vars=plot.trace[[1L]]$range, draws=i0:index)
          yvalues <- get_from(out[[plot.trace[[2L]]$name]], vars=plot.trace[[2L]]$range, draws=i0:index)
          for (ch in chains) {
            if (ch == 1L && (!add.to.plot || (add.to.plot && i0 == 1L)))
              plot(xvalues[[ch]], yvalues[[ch]], type=plot.type, xlim=range(unlst(xvalues)), ylim=range(unlst(yvalues)), xlab=plot.trace.names[1L], ylab=plot.trace.names[2L])
            else
              lines(xvalues[[ch]], yvalues[[ch]], type=plot.type, col=ch)
          }
        } else {  # pairs plot for >=3 variables
          plot.matrix <- matrix(NA_real_, n.chain * (index - i0 + 1L), length(plot.trace))
          for (v in seq_along(plot.trace))
            plot.matrix[, v] <- unlst(get_from(out[[plot.trace[[v]]$name]], vars=plot.trace[[v]]$range, draws=i0:index))
          pairs(plot.matrix, labels=plot.trace.names, cex.labels=1, gap=0, pch=20L, col=rep_each(chains, index - i0 + 1L))
        }
        Sys.sleep(0)  # forces the plot to update
      }  # END if (!is.null(plot.trace))
    }  # END if (i %% n.progress == 0L)
  }  # END for (i in seq_len(n.iter))

  # finalize output
  if (MH) {
    for (a in seq_along(MHpars))
      out[["_accept"]][[a]] <- lapply(out[["_accept"]][[a]], function(x) x / n.iter)
  }
  if (write.to.file) for (v in to.file) close(outfile[[v]])
  tryCatch({
    out[["_state"]] <- p  # store the final state(s) p
    if (!is.null(store.mean)) {
      out[["_means"]] <- lapply(out[["_means"]], function(x) lapply(x, function(y) y / n.draw))
      if (store.sds) {
        out[["_sds"]] <- lapply(out[["_sds"]], function(x) lapply(x, function(y) y / n.draw))
        for (v in store.mean) {
          for (ch in seq_along(out[["_sds"]][[v]]))
            out[["_sds"]][[v]][[ch]] <- sqrt(out[["_sds"]][[v]][[ch]] - out[["_means"]][[v]][[ch]]^2)
        }
      }
    }
    # class attributes set at the end of the simulation for faster assignment
    class(out) <- "mcdraws"
    for (v in store) {
      if (!is.null(sampler$coef.names))
        attr(out[[v]], "labels") <- sampler$coef.names[[v]]
      if (is.null(attr(out[[v]], "labels")))
        if (n_vars(out[[v]]) == 1L)
          attr(out[[v]], "labels") <- v
        else
          attr(out[[v]], "labels") <- as.character(seq_len(n_vars(out[[v]])))
      class(out[[v]]) <- "dc"  # draws component class
    }
    sim.time <- proc.time() - started
    out[["_info"]]$time <- as.numeric(sim.time[1L] + if (!is.na(sim.time[4L])) sim.time[4L] else 0)
    if (verbose) {
      cat("\n\n")
      print(sim.time)
    }
  }, error=function(e) {
    print(e)
    cat("Error during finalization of output. Raw results returned.\n")
  })
  out
}

#' Convert a draws component object to another format
#'
#' Use \code{to_mcmc} to convert a draws component to class \code{\link[coda]{mcmc.list}},
#' allowing one to use MCMC diagnostic functions provided by package \pkg{coda}.
#' Use \code{as.array} to convert to an array of dimension \code{(draws, chains, parameters)}.
#' The array format is supported by some packages for analysis or visualisation of MCMC
#' simulation results, e.g. \pkg{bayesplot}.
#' Use \code{as.matrix} to convert to a matrix, concatenating the chains.
#' Finally, use \code{to_draws_array} to convert either a draws component or
#' (a subset of components of) an mcdraws object to a \code{draws_array} object
#' as defined in package \pkg{posterior}.
#'
#' @examples
#' \donttest{
#' data(iris)
#' sampler <- create_sampler(Sepal.Length ~ reg(~ Petal.Length + Species, name="beta"), data=iris)
#' sim <- MCMCsim(sampler, burnin=100, n.chain=2, n.iter=400)
#' summary(sim)
#' if (require("coda", quietly=TRUE)) {
#'   mcbeta <- to_mcmc(sim$beta)
#'   geweke.diag(mcbeta)
#' }
#' if (require("posterior", quietly=TRUE)) {
#'   mcbeta <- to_draws_array(sim$beta)
#'   mcbeta
#'   draws <- to_draws_array(sim)
#'   str(draws)
#' }
#' str(as.array(sim$beta))
#' str(as.matrix(sim$beta))
#'
#' # generate some example data
#' n <- 250
#' dat <- data.frame(x=runif(n), f=as.factor(sample(1:5, n, replace=TRUE)))
#' gd <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1), name="beta"), data=dat)
#' dat$y <- gd$y
#' sampler <- create_sampler(y ~ reg(~ x + f, name="beta"), data=dat)
#' sim <- MCMCsim(sampler, n.chain=2, n.iter=400)
#' str(sim$beta)
#' str(as.array(sim$beta))
#' bayesplot::mcmc_hist(as.array(sim$beta))
#' bayesplot::mcmc_dens_overlay(as.array(sim$beta))
#' # fake data simulation check:
#' bayesplot::mcmc_recover_intervals(as.array(sim$beta), gd$pars$beta)
#' bayesplot::mcmc_recover_hist(as.array(sim$beta), gd$pars$beta)
#'
#' ex <- mcmcsae_example()
#' plot(ex$dat$fT, ex$dat$y)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, n.chain=2, n.iter=400, store.all=TRUE)
#' str(sim$beta)
#' str(as.matrix(sim$beta))
#' # fake data simulation check:
#' bayesplot::mcmc_recover_intervals(as.matrix(sim$beta), ex$pars$beta)
#' bayesplot::mcmc_recover_intervals(as.matrix(sim$u), ex$pars$u)
#' }
#'
#' @param x a component of an mcdraws object corresponding to a scalar or vector model parameter.
#' @param components optional character vector of names of draws components in an mcdraws object.
#'  This can be used to select a subset of components to convert to
#'   \code{\link[posterior]{draws_array}} format.
#' @param colnames whether column names should be set.
#' @param ... currently ignored.
#' @returns The draws component(s) coerced to an \code{\link[coda]{mcmc.list}} object,
#'  a \code{\link[posterior]{draws_array}} object, an array, or a matrix.
#' @name MCMC-object-conversion
NULL

#' @export
#' @rdname MCMC-object-conversion
to_mcmc <- function(x) {
  if (!inherits(x, "dc")) stop("'x' should be an object of class 'dc'")
  n.draw <- n_draws(x)
  out <- lapply(x, `attr<-`, which="mcpar", value=c(1L, n.draw, 1L))
  out <- lapply(out, `dimnames<-`, value=list(NULL, labels(x)))
  # turn each chain into mcmc object
  out <- lapply(out, `attr<-`, which="class", value="mcmc")
  class(out) <- "mcmc.list"
  out
}

#' @export
#' @rdname MCMC-object-conversion
to_draws_array <- function(x, components=NULL) {
  if (inherits(x, "mcdraws")) {
    f <- function(x, name) {
      out <- x[[name]]
      if (n_vars(out) > 1L) labels(out) <- paste(name, labels(out), sep="_")
      to_draws_array(out)
    }
    if (is.null(components)) components <- par_names(x)
    out <- f(x, components[1L])
    for (i in seq_along(components[-1L]))
      out <- posterior::bind_draws(out, f(x, components[i + 1L]), along="variable")
    out
  } else if (inherits(x, "dc")) {
    posterior::as_draws_array(as.array.dc(x))
  } else {
    stop("not an object of class 'mcdraws' or 'dc'")
  }
}

#' @method as.array dc
#' @export
#' @rdname MCMC-object-conversion
as.array.dc <- function(x, ...) {
  chains <- seq_len(n_chains(x))
  out <- array(NA_real_, dim=c(n_draws(x), n_chains(x), n_vars(x)), dimnames=list(NULL, paste("chain", chains, sep=":"), labels(x)))
  for (ch in chains) out[, ch, ] <- x[[ch]]
  out
}

#' @method as.matrix dc
#' @export
#' @rdname MCMC-object-conversion
as.matrix.dc <- function(x, colnames=TRUE, ...) {
  dimx <- dim(x[[1L]])
  nv <- dimx[2L]
  if (nv < 100L) {
    out <- do.call("rbind", x)
  } else {
    nch <- length(x)
    nit <- dimx[1L]
    out <- matrix(NA_real_, nch*nit, nv)
    for (ch in seq_len(nch))
      out[((ch - 1L)*nit + 1L):(ch*nit), ] <- x[[ch]]
  }
  if (colnames) dimnames(out) <- list(NULL, attr(x, "labels"))
  out
}

#' Get the parameter names from an mcdraws object
#' 
#' @examples
#' data(iris)
#' sampler <- create_sampler(Sepal.Length ~ 
#'     reg(~ Petal.Length + Species, name="beta"), data=iris)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=400)
#' (summary(sim))
#' par_names(sim)
#'
#' @export
#' @param obj an mcdraws object.
#' @returns The names of the parameters whose MCMC simulations are stored in \code{obj}.
par_names <- function(obj) obj[["_info"]][["parnames"]]

#' Read MCMC draws from a file
#'
#' Read draws written to file by \code{\link{MCMCsim}} used with argument \code{to.file}.
#'
#' @examples
#' \dontrun{
#' # NB this example creates a file "MCdraws_e_.dat" in the working directory
#' n <- 100
#' dat <- data.frame(x=runif(n), f=as.factor(sample(1:5, n, replace=TRUE)))
#' gd <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1), name="beta"), data=dat)
#' dat$y <- gd$y
#' sampler <- create_sampler(y ~ reg(~ x + f, name="beta"), data=dat)
#' # run the MCMC simulation and write draws of residuals to file:
#' sim <- MCMCsim(sampler, n.iter=500, to.file="e_")
#' summary(sim)
#' mcres <- read_draws("e_")
#' summary(mcres)
#' }
#'
#' @export
#' @param name name of the parameter to load the corresponding file with posterior draws for.
#' @param filename name of the file in which the draws are stored.
#' @returns An object of class \code{dc} containing MCMC draws for a (vector) parameter.
# TODO read subsets
read_draws <- function(name, filename=paste0("MCdraws_", name, ".dat")) {
  con <- file(filename, "rb")
  info <- read_header(con)
  n.iter <- info[["n.iter"]]
  n.chain <- info[["n.chain"]]
  n.par <- info[["n.par"]]
  read.size <- if (info[["single"]]) 4L else NA_integer_
  res <- array(readBin(con, what="numeric", n=n.iter*n.chain*n.par, size=read.size), dim=c(n.par, n.chain, n.iter))
  close(con)
  out <- list()
  for (i in seq_len(n.chain))
    out[[i]] <- matrix(res[, i, ], ncol=n.par, byrow=TRUE)
  if (!is.null(info$parnames)) attr(out, "labels") <- info$parnames
  class(out) <- "dc"
  out
}

#' Select chains, samples and parameters from a draws component (dc) object
#' 
#' @noRd
#' @param dc an object of class \code{dc}.
#' @param vars an integer vector indicating which parameters to select.
#' @param chains an integer vector indicating which chains to select.
#' @param draws an integer vector indicating which samples to select.
#' @returns The selected part of the draws component. The output object's class is \code{dc}.
# for internal use; see below subset.dc for external use
get_from <- function(dc, vars=NULL, chains=NULL, draws=NULL) {
  if (is.null(draws) && is.null(vars)) {
    if (is.null(chains))
      dc
    else
      dc[chains]
  } else if (is.null(vars)) {
    if (is.null(chains))
      lapply(dc, function(x) x[draws, , drop=FALSE])
    else
      lapply(dc[chains], function(x) x[draws, , drop=FALSE])
  } else if (is.null(draws)) {
    if (is.null(chains))
      lapply(dc, function(x) x[, vars, drop=FALSE])
    else
      lapply(dc[chains], function(x) x[, vars, drop=FALSE])
  } else {
    if (is.null(chains))
      lapply(dc, function(x) x[draws, vars, drop=FALSE])
    else
      lapply(dc[chains], function(x) x[draws, vars, drop=FALSE])
  }
}

#' Select a subset of chains, samples and parameters from a draws component (dc) object
#'
#' @examples
#' n <- 300
#' dat <- data.frame(x=runif(n), f=as.factor(sample(1:7, n, replace=TRUE)))
#' gd <- generate_data(~ reg(~ x + f, prior=pr_normal(precision=1), name="beta"), data=dat)
#' dat$y <- gd$y
#' sampler <- create_sampler(y ~ reg(~ x + f, name="beta"), data=dat)
#' sim <- MCMCsim(sampler)
#' (summary(sim$beta))
#' (summary(subset(sim$beta, chains=1)))
#' (summary(subset(sim$beta, chains=1, draws=sample(1:n_draws(sim), 100))))
#' (summary(subset(sim$beta, vars=1:2)))
#'
#' @export
#' @method subset dc
#' @param x a draws component (dc) object.
#' @param chains an integer vector indicating which chains to select.
#' @param draws an integer vector indicating which samples to select.
#' @param vars an integer vector indicating which parameters to select.
#' @param ... not used.
#' @returns The selected part of the draws component as an object of class \code{dc}.
subset.dc <- function(x, chains=NULL, draws=NULL, vars=NULL, ...) {
  labs <- attr(x, "labels")
  out <- get_from(x, vars=vars, chains=chains, draws=draws)
  if (!is.null(labs)) attr(out, "labels") <- if (is.null(vars)) labs else labs[vars]
  class(out) <- "dc"
  out
}

#' Return variable name and range of a composite name
#'
#' @noRd
#' @param composite.name a character string such as \code{beta[2:5]}
#' @param default.first if \code{TRUE} the range will default to 1 and otherwise to \code{NULL}
#' @returns A list with name and range as separate components.
get_name_range <- function(composite.name, default.first=TRUE) {
  spl <- strsplit(composite.name, "\\[[[:space:]]*")[[1L]]
  if (length(spl) > 1L)
    range <- as.integer(eval(str2lang(strsplit(spl[2L], "\\]")[[1L]])))
  else
    range <- if (default.first) 1L else NULL
  list(name=spl[1L], range=range)
}

#' Summarize a draws component (dc) object
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' summary(sim$u)
#' }
#'
#' @export
#' @method summary dc
#' @param object an object of class \code{dc}.
#' @param probs vector of probabilities at which to evaluate quantiles.
#' @param na.rm whether to remove NA/NaN draws in computing the summaries.
# TODO allow removing NA/NaNs also from n_eff and R-hat computation
#' @param time MCMC computation time; if specified the effective sample size per unit of time
#'  is returned in an extra column labeled 'efficiency'.
#' @param abbr if \code{TRUE} abbreviate the labels in the output.
#' @param batch.size number of parameter columns to process simultaneously. A larger batch size may speed things up a little,
#'  but if an out of memory error occurs it may be a good idea to use a smaller number and try again. The default is 100.
#' @param ... arguments passed to \code{\link{n_eff}}.
#' @returns A matrix with summaries of class \code{dc_summary}.
summary.dc <- function(object, probs=c(0.05, 0.5, 0.95), na.rm=FALSE, time=NULL, abbr=FALSE, batch.size=100L, ...) {
  col_names <- c("Mean", "SD", "t-value", "MCSE", paste0("q", probs), "n_eff")
  if (!is.null(time)) col_names <- c(col_names, "efficiency")
  if (n_chains(object) >= 2L) col_names <- c(col_names, "R_hat")
  nv <- n_vars(object)
  labs <- labels(object)
  if (is.null(labs) && nv > 1L) labs <- seq_len(nv)
  if (length(labs) != nv) labs <- NULL  # be error tolerant
  if (abbr && !is.null(labs) && (max(nchar(labs)) > 10L)) labs <- abbreviate(labs, minlength=8L)
  out <- matrix(NA_real_, nrow=nv, ncol=length(col_names), dimnames=list(labs, col_names))
  # for some summary statistics we use as.matrix, which can easily use too much memory for large vector parameters
  # n_eff may also use too much memory for large vector parameters
  #   --> split into chunks of 100 columns
  batch <- seq_len(batch.size)
  n.batch <- nv %/% batch.size + (nv %% batch.size > 0L)
  for (i in seq_len(n.batch)) {
    if ((i == n.batch) && (nv %% batch.size))
      ind <- (nv - (nv %% batch.size) + 1L):nv
    else
      ind <- batch + (i - 1L) * batch.size
    sim <- get_from(object, vars=ind)
    out[ind, "n_eff"] <- n_eff(sim, ...)
    if (any("R_hat" == col_names)) out[ind, "R_hat"] <- R_hat(sim)
    sim <- as.matrix.dc(sim, colnames=FALSE)
    out[ind, "Mean"] <- .colMeans(sim, nrow(sim), ncol(sim), na.rm=na.rm)
    out[ind, "SD"] <- colSds(sim, na.rm=na.rm)
    out[ind, paste0("q", probs)] <- colQuantiles(sim, probs=probs, na.rm=na.rm)
  }
  out[, "t-value"] <- out[, "Mean"] / out[, "SD"]
  out[, "MCSE"] <- out[, "SD"] / sqrt(out[, "n_eff"])
  if (!is.null(time)) out[, "efficiency"] <- out[, "n_eff"] / time
  class(out) <- "dc_summary"
  out
}

#' Summarize an mcdraws object
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' summary(sim)
#' par_names(sim)
#' summary(sim, c("beta", "v_sigma", "u_sigma"))
#' }
#'
#' @export
#' @method summary mcdraws
#' @param object an object of class \code{mcdraws}, typically generated by function \code{\link{MCMCsim}}.
#' @param vnames optional character vector to select a subset of parameters.
#' @param probs vector of probabilities at which to evaluate quantiles.
#' @param na.rm whether to remove NA/NaN draws in computing the summaries.
#' @param efficiency if \code{TRUE} the effective sample size per second of computation time is returned as well.
#' @param abbr if \code{TRUE} abbreviate the labels in the output.
#' @param batch.size number of parameter columns to process simultaneously for vector parameters. A larger batch size
#'  may speed things up a little, but if an out of memory error occurs it may be a good idea to use a smaller number
#'  and try again. The default is 100.
#' @param ... arguments passed to \code{\link{n_eff}}.
#' @returns A list of class \code{mcdraws_summary} summarizing \code{object}.
summary.mcdraws <- function(object, vnames=NULL, probs=c(0.05, 0.5, 0.95), na.rm=FALSE, efficiency=FALSE, abbr=FALSE, batch.size=100L, ...) {
  if (is.null(vnames)) vnames <- par_names(object)
  time <- if (efficiency) object[["_info"]]$time else NULL
  out <- list()
  for (v in vnames) {
    if (is.null(object[[v]])) {
      warn("parameter ", v, " not found")
      next
    }
    if (is.null(object[["_cluster"]]) || any("cl" == names(list(...)))) {
      out[[v]] <- summary(object[[v]], probs, na.rm, time, abbr, batch.size=batch.size, ...)
    } else {
      # if a cluster is available, it is passed to n_eff (in which usually most time is spent)
      out[[v]] <- summary(object[[v]], probs, na.rm, time, abbr, batch.size=batch.size, cl=object[["_cluster"]], ...)
    }
  }
  class(out) <- "mcdraws_summary"
  out
}

#' Display a summary of a \code{dc} object
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' print(summary(sim$u), sort="n_eff")
#' }
#'
#' @export
#' @method print dc_summary
#' @param x an object of class \code{dc_summary}.
#' @param digits number of digits to use, defaults to 3.
#' @param max.lines maximum number of lines to display.
#'  If \code{NULL}, all elements are displayed.
#' @param tail if \code{TRUE} the last instead of first at most \code{max.lines} are displayed.
#' @param sort column name on which to sort the output.
#' @param max.label.length if specified, printed row labels will be abbreviated to at most this length.
#' @param ... passed on to \code{print.default}.
print.dc_summary <- function(x, digits=3L, max.lines=1000L, tail=FALSE, sort=NULL, max.label.length=NULL, ...) {
  nlines <- if (is.null(max.lines)) nrow(x) else min(max.lines, nrow(x))
  if (tail)
    rows <- (nrow(x)):(nrow(x) - nlines + 1L)
  else
    rows <- seq_len(nlines)
  maxpr <- getOption("max.print")
  if (nlines*ncol(x) > maxpr) {
    on.exit(options(max.print=maxpr))
    options(max.print=nlines*ncol(x))
  }
  if (is.null(sort)) {
    z <- x[rows, , drop=FALSE]
  } else {
    if (!is.character(sort) || !all(sort %in% colnames(x))) stop("'sort' should be a column name of 'x'")
    z <- x[order(x[, sort])[rows], , drop=FALSE]
  }
  if (!is.null(max.label.length))
    dimnames(z)[[1L]] <- abbreviate(rownames(z), minlength=max.label.length, strict=TRUE)
  print(z, digits=digits, ...)
  suppressed <- nrow(x) - nlines
  if (suppressed > 0L) cat("...", suppressed, "elements suppressed ...\n")
}

#' Print a summary of MCMC simulation results
#'
#' Display a summary of an \code{mcdraws} object, as output by \code{\link{MCMCsim}}.
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, store.all=TRUE)
#' print(summary(sim), sort="n_eff")
#' }
#'
#' @export
#' @method print mcdraws_summary
#' @param x an object of class \code{mcdraws_summary} as output by \code{\link{summary.mcdraws}}.
#' @param digits number of digits to use, defaults to 3.
#' @param max.lines maximum number of elements per vector parameter to display.
#'  If \code{NULL}, all elements are displayed.
#' @param tail if \code{TRUE} the last instead of first \code{max.lines} of each component
#'  are displayed.
#' @param sort column name on which to sort the output.
#' @param ... passed on to \code{print.default}.
print.mcdraws_summary <- function(x, digits=3L, max.lines=10L, tail=FALSE, sort=NULL, ...) {
  for (i in seq_along(x)) {
    cat(names(x)[i], ":\n")
    print(x[[i]], digits=digits, max.lines=max.lines, tail=tail, sort=sort, ...)  # print an object of class 'dc'
    cat("\n")
  }
}

#' Get and set the variable labels of a draws component object for a vector-valued parameter
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=50, n.iter=100, n.chain=1, store.all=TRUE)
#' labels(sim$beta)
#' labels(sim$v)
#' labels(sim$beta) <- c("a", "b")
#' labels(sim$beta)
#' }
#'
## @method labels dc
#' @param object a draws component object.
#' @param value a vector of labels.
#' @param ... currently not used.
#' @returns The extractor function returns the variable labels.
#' @name labels
NULL

#' @export
#' @rdname labels
labels.dc <- function(object, ...) attr(object, "labels")

#' @export
#' @rdname labels
`labels<-` <- function(object, value) {
  if (class(object)[1L] == "dc") {
    if (length(value) != n_vars(object)) stop("wrong length of labels vector")
  }
  attr(object, "labels") <- as.character(value)
  object
}

#' Get the number of chains, samples per chain or the number of variables in a simulation object
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=50)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=5, store.all=TRUE)
#' n_chains(sim); n_chains(sim$beta)
#' n_draws(sim); n_draws(sim$beta)
#' n_vars(sim$beta); n_vars(sim$sigma_); n_vars(sim$llh_); n_vars(sim$v)
#' plot(sim, "beta")
#' n_chains(subset(sim$beta, chains=1:2))
#' n_draws(subset(sim$beta, draws=sample(1:n_draws(sim), 100)))
#' n_vars(subset(sim$u, vars=1:2))
#' }
#'
#' @param obj an mcdraws object or a draws component (dc) object.
#' @param dc a draws component object.
#' @returns The number of chains or retained samples per chain or
#'  the number of variables.
#' @name n_chains-n_draws-n_vars
NULL

#' @export
#' @rdname n_chains-n_draws-n_vars
n_chains <- function(obj) {
  switch(class(obj)[1L],
    mcdraws = obj[["_info"]][["n.chain"]],
    dc=, list = length(obj),  # NB object not necessarily (yet) of class 'dc' when called from MCMCsim
    stop("unexpected input")
  )
}

#' @export
#' @rdname n_chains-n_draws-n_vars
n_draws <- function(obj) {
  switch(class(obj)[1L],
    mcdraws = obj[["_info"]][["n.draw"]],
    dc=, list = dim(obj[[1L]])[1L],
    stop("unexpected input")
  )
}

#' @export
#' @rdname n_chains-n_draws-n_vars
n_vars <- function(dc) dim(dc[[1L]])[2L]

#' Compute MCMC diagnostic measures
#'
#' \code{R_hat} computes Gelman-Rubin convergence diagnostics based on the MCMC output
#' in a model component, and \code{n_eff} computes the effective sample sizes, .i.e.
#' estimates for the number of independent samples from the posterior distribution.
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=4, store.all=TRUE)
#' n_eff(sim$beta)
#' n_eff(sim$v_sigma)
#' n_eff(sim$v_rho)
#' R_hat(sim$beta)
#' R_hat(sim$llh_)
#' R_hat(sim$v_sigma)
#' }
#'
#' @param dc a draws component (dc) object corresponding to a model parameter.
#' @param useFFT whether to use the Fast Fourier Transform algorithm. Default is \code{TRUE} as this is typically faster.
#' @param lag.max the lag up to which autocorrelations are computed in case \code{useFFT=FALSE}.
#' @param cl a cluster for parallel computation.
#' @returns In case of \code{R_hat} the split-R-hat convergence diagnostic for each
#'  component of the vector parameter, and in case of \code{n_eff} the effective
#'  number of independent samples for each component of the vector parameter.
#' @references
#'  A. Gelman and D. B. Rubin (1992).
#'    Inference from Iterative Simulation Using Multiple Sequences.
#'    Statistical Science 7, 457-511.
#'
#'  A. Gelman, J.B. Carlin, H.S. Stern, D.B. Dunson, A. Vehtari and D.B. Rubin (2013).
#'    Bayesian Data Analysis, 3rd edition.
#'    Chapman & Hall/CRC.
#' @name MCMC-diagnostics
NULL

#' @export
#' @rdname MCMC-diagnostics
# adapted from code in R package rstan
R_hat <- function(dc) {
  n.var <- n_vars(dc)
  n.draw <- n_draws(dc)
  if (n.draw < 1L) stop("no draws in simulation object")
  n.chain <- n_chains(dc)
  if (n.chain < 2L) stop("at least 2 chains required for R_hat diagnostic")
  split <- n.draw %/% 2L  # split in two halves
  part1 <- seq_len(split)
  part2 <- (split + 1L):n.draw
  split_chain_means <- split_chain_vars <- matrix(NA_real_, 2L * n.chain, n.var)
  for (i in seq_len(n.chain)) {
    split_chain_means[i,] <- .colMeans(dc[[i]][part1, ,drop=FALSE], split, n.var)
    split_chain_vars[i,] <- colVars(dc[[i]][part1, ,drop=FALSE])
    split_chain_means[n.chain + i,] <- .colMeans(dc[[i]][part2, ,drop=FALSE], length(part2), n.var)
    split_chain_vars[n.chain + i,] <- colVars(dc[[i]][part2, ,drop=FALSE])
  }
  var_between <- split * colVars(split_chain_means)
  var_within <- .colMeans(split_chain_vars, 2L * n.chain, n.var)
  # need to set names as colVars drops them
  setNames(sqrt((var_between / var_within + split - 1) / split), labels.dc(dc))
}

#' @export
#' @rdname MCMC-diagnostics
# based on code from package rstan
n_eff <- function(dc, useFFT=TRUE, lag.max, cl=NULL) {
  n.chain <- n_chains(dc)
  n.var <- n_vars(dc)
  n.draw <- n_draws(dc)
  if (n.draw < 2L) return(rep.int(NA_real_, n.var))
  if (missing(lag.max)) {
    lag.max <- n.draw - 1L
  } else {
    if (useFFT) {
      warn("'lag.max' argument ignored because 'useFFT=TRUE'")
      lag.max <- n.draw - 1L
    } else {
      lag.max <- min(lag.max, n.draw - 1L)
    }
  }
  # function to compute autocovariance function summed over (a subset of) chains
  compute_acov_sum <- function(dc) {
    n.chain <- n_chains(dc)
    acov <- matrix(0, lag.max + 1L, n.var, dimnames=list(NULL, labels.dc(dc)))
    # TODO allow to use split chains
    for (ch in seq_len(n.chain)) {
      if (useFFT) {
        acov <- acov + ac_fft(dc[[ch]])
      } else {
        acov <- acov + apply(dc[[ch]], 2L, function(x) {
          acf(x, lag.max=lag.max, plot=FALSE, type="covariance")$acf[, , 1L]
        })
      }
    }
    acov
  }
  if (is.null(cl))
    acov <- compute_acov_sum(dc) / n.chain
  else
    acov <- Reduce(`+`, parallel::parLapply(cl, split_chains(dc, length(cl)), compute_acov_sum)) / n.chain
  chain_means <- vapply(dc, .colMeans, numeric(n.var), n.draw, n.var)  # n.var x n.chain matrix
  if (n.var == 1L) chain_means <- matrix(chain_means, ncol=n.chain)
  var_plus <- acov[1L, ]
  mean_var <- var_plus * n.draw / (n.draw - 1L)
  if (n.chain > 1L) var_plus <- var_plus + rowVars(chain_means)
  rho_hat <- 1 - (mean_var - t.default(acov[-1L, , drop=FALSE])) / var_plus
  rho_hat[is.nan(rho_hat)] <- 0
  # TODO check criterion until what lag to sum correlations (BDA3, or Geyer)
  rho_hat_sum <- apply(rho_hat, 1L, function(x) sum(x * !cumsum(x < 0)))
  ess <- n.chain * n.draw
  ess <- ess / (1 + 2 * rho_hat_sum)
  ess
}

#' Compute autocovariance or autocorrelation function via Wiener-Khinchin theorem
#'   using Fast Fourier Transform
#'
#' @keywords internal
#' @param x matrix with time (iteration number) along the rows and variables along the columns.
#' @param demean whether to subtract from each column its mean.
#' @returns A matrix of the same size as x with autocovariances at all lags from 1 to the number of rows.
# TODO check whether zero-padding to (approx.) a power of 2 is possible (should be faster); use nextn()?
ac_fft <- function(x, demean=TRUE) {
  nr <- nrow(x)
  if (demean) x <- x - rep_each(.colMeans(x, nr, ncol(x)), nr)
  Fx <- mvfft(rbind(x, 0*x)) / sqrt(2*nr)  # zero-pad and FFT
  Re(mvfft(Fx * Conj(Fx), inverse=TRUE))[seq_len(nr), , drop=FALSE] / nr
}

#' Return Metropolis-Hastings acceptance rates
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' # specify a model that requires MH sampling (in this case for a modeled
#' #   degrees of freedom parameter in the variance part of the model)
#' sampler <- create_sampler(ex$model, data=ex$dat, formula.V=~vfac(factor="fA",
#'   prior=pr_invchisq(df="modeled")))
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=4, store.all=TRUE)
#' (summary(sim))
#' acceptance_rates(sim)
#' }
#'
#' @export
#' @param obj an mcdraws object, i.e. the output of function \code{\link{MCMCsim}}.
#' @param aggregate.chains whether to return averages over chains or results per chain.
#' @returns A list of acceptance rates.
acceptance_rates <- function(obj, aggregate.chains=FALSE) {
  if (aggregate.chains)
    lapply(obj[["_accept"]], function(x) Reduce("+", x)/length(x))
  else 
    obj[["_accept"]]
}

#' Get means or standard deviations of parameters from the MCMC output in an mcdraws object
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=50)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=4)
#' get_means(sim)
#' get_means(sim, "e_")
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=4,
#'   store.mean=c("beta", "u"), store.sds=TRUE)
#' summary(sim, "beta")
#' get_means(sim, "beta")
#' get_sds(sim, "beta")
#' get_means(sim, "u")
#' get_sds(sim, "u")
#' }
# note: get_sds returns slightly different values than summary
# TODO more stable online variance computation
#'
#' @param obj an object of class \code{mcdraws}.
#' @param vnames optional character vector to select a subset of parameters.
#' @returns A list with simulation means or standard deviations.
#' @name posterior-moments
NULL

#' @export
#' @rdname posterior-moments
get_means <- function(obj, vnames=NULL) {
  # first extract mean only components and average over chains
  out <- obj[["_means"]]
  if (length(out) && !is.null(vnames)) out <- out[names(out) %in% vnames]
  nc <- n_chains(obj)
  if (length(out))
    for (v in seq_along(out)) out[[v]] <- Reduce(`+`, out[[v]]) / nc
  # append posterior means of draws components
  if (is.null(vnames)) vnames <- par_names(obj)
  vnames <- vnames[!(vnames %in% names(out)) & vnames %in% par_names(obj)]
  ni <- n_draws(obj)  # draws per chain
  c(out, lapply(obj[vnames], function(x) Reduce(`+`, lapply(x, function(ch) .colMeans(ch, ni, ncol(ch)))) / nc))
}

#' @export
#' @rdname posterior-moments
get_sds <- function(obj, vnames=NULL) {
  out <- obj[["_sds"]]
  if (length(out) && !is.null(vnames)) out <- out[names(out) %in% vnames]
  if (length(out)) {
    nc <- n_chains(obj)
    ni <- n_draws(obj)
    # include between-chain component
    for (v in names(out))
      out[[v]] <- sqrt((1/(nc*ni - 1L))*((ni - 1L)*rowSums(do.call(cbind, out[[v]])^2) + ni*(nc - 1L)*apply(do.call(cbind, obj[["_means"]][[v]]), 1L, var)))
  }
  # append posterior sds of draws components
  if (is.null(vnames)) vnames <- par_names(obj)
  vnames <- vnames[!(vnames %in% names(out)) & vnames %in% par_names(obj)]
  c(out, lapply(obj[vnames], function(x) unname(colSds(as.matrix.dc(x, colnames=FALSE)))))
}

#' Extract a list of parameter values for a single draw
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=50)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=4, store.all=TRUE)
#' get_draw(sim, iter=20, chain=3)
#' }
#'
#' @export
#' @param obj an object of class \code{mcdraws}.
#' @param iter iteration number.
#' @param chain chain number.
#' @returns A list with all parameter values of draw \code{iter} from chain \code{chain}.
get_draw <- function(obj, iter, chain) {
  p <- list()
  for (v in obj[["_info"]][["parnames"]])
    p[[v]] <- obj[[v]][[chain]][iter, ]
  if (!is.null(obj[["_info"]][["list.pars"]])) for (v in obj[["_info"]][["list.pars"]])
    p[[v]] <- obj[[v]][[chain]][[iter]]
  p
}

#' Transform one or more draws component objects into a new one by applying a function
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example(n=50)
#' sampler <- create_sampler(ex$model, data=ex$dat)
#' sim <- MCMCsim(sampler, burnin=100, n.iter=300, thin=2, n.chain=4, store.all=TRUE)
#' summary(sim$v_sigma)
#' summary(transform_dc(sim$v_sigma, fun=function(x) x^2))
#' summary(transform_dc(sim$u, sim$u_sigma, fun=function(x1, x2) abs(x1)/x2))
#' }
#'
#' @export
#' @param ... draws component object(s) of class \code{dc}.
#' @param fun a function to apply. This function should take as many arguments as there are input objects.
#'  The arguments can be arbitrarily named, but they are assumed to be in the same order as the input objects.
#'  The function should return a vector.
#' @param to.matrix if \code{TRUE} the output is in matrix format; otherwise it is a draws component object.
#' @param labels optional labels for the output object.
#' @returns Either a matrix or a draws component object.
# TODO add option chain.wise to alow application of a vectorised function to each chain at once
transform_dc <- function(..., fun, to.matrix=FALSE, labels=NULL) {
  objs <- list(...)
  nobj <- length(objs)
  nc <- n_chains(objs[[1L]])
  ni <- n_draws(objs[[1L]])
  fun <- match.fun(fun)
  if (length(formals(args(fun))) != nobj) stop("'fun' must have as many arguments as there are input objects")
  # TODO check that all objs have same number of chains, draws
  test <- do.call(fun, lapply(objs, function(x) x[[1L]][1L, ]))
  if (!is.vector(test)) stop("'fun' should return a vector")
  no <- length(test)
  if (!is.null(labels) && length(labels) != no) stop("length of 'labels' does not match dimension of 'fun'")
  if (to.matrix) {
    if (is.null(labels))
      out <- matrix(NA_real_, nc*ni, no)
    else
      out <- matrix(NA_real_, nc*ni, no, dimnames=list(NULL, labels))
    k <- 1L
    for (ch in seq_len(nc))
      for (i in seq_len(ni)) {
        out[k, ] <- do.call(fun, lapply(objs, function(x) x[[ch]][i, ]))
        k <- k + 1L
      }
  } else {
    out <- list()
    for (ch in seq_len(nc)) {
      out[[ch]] <- matrix(NA_real_, ni, no)
      for (i in seq_len(ni)) out[[ch]][i, ] <- do.call(fun, lapply(objs, function(x) x[[ch]][i, ]))
    }
    if (!is.null(labels)) attr(out, "labels") <- labels
    class(out) <- "dc"
  }
  out
}
