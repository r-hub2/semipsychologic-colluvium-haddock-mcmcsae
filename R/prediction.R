#' Generate draws from the predictive distribution
#'
#' @examples
#' \donttest{
#' n <- 250
#' dat <- data.frame(x=runif(n))
#' dat$y <- 1 + dat$x + rnorm(n)
#' sampler <- create_sampler(y ~ x, data=dat)
#' sim <- MCMCsim(sampler)
#' summary(sim)
#' # in-sample prediction
#' pred <- predict(sim, ppcheck=TRUE)
#' hist(attr(pred, "ppp"))
#' # out-of-sample prediction
#' pred <- predict(sim, newdata=data.frame(x=seq(0, 1, by=0.1)))
#' summary(pred)
#' }
#'
#' @export
## @method predict mcdraws
#' @param object an object of class \code{mcdraws}, as output by \code{\link{MCMCsim}}.
#' @param newdata data frame with auxiliary information to be used for prediction.
#  use X. instead of X to avoid argument name clash in parallel version par_predict with parLapply
#' @param X. a list of design matrices; alternatively, \code{X.} equals 'in-sample' or 'linpred'.
#'  If 'in-sample' (the default if newdata is not supplied), the design matrices for in-sample
#'  prediction are used. If 'linpred' the 'linpred_' component of \code{object} is used.
#' @param type the type of predictions. The default is \code{"data"}, meaning that
#'  new data is generated according to the predictive distribution.
#'  If \code{type="link"} only the linear predictor for the mean is generated, and
#'  in case \code{type="response"} the linear predictor is transformed to the response scale.
#'  For Gaussian models \code{type="link"} and \code{type="response"} are equivalent.
#'  For binomial models \code{type="response"} returns the simulations
#'  of the latent probabilities, and for negative binomial models the exponentiated
#'  linear predictor, which differs from the mean by a factor equal to the shape parameter.
#'  For multinomial models \code{type="link"} generates
#'  the linear predictor for all categories except the last, and \code{type="response"} transforms
#'  this vector to the probability scale, and \code{type="data"} generates the multinomial data,
#'  all in long vector format, where the output for all categories (except the last) are stacked.
#'  For multinomial models and single trials, a further option is \code{type="data_cat"},
#'  which generates the data as a categorical vector, with integer coded levels.
#' @param weights an optional formula specifying a vector of nonnegative weights.
#'  Weights are only used for data-generating prediction (\code{type = "data"}).
#'  A weight \eqn{w_i} means that the \emph{sum} of \eqn{w_i} individual
#'  predictions is generated for each unit $i$. For example, for the poststratification
#'  step of Multilevel Regression and Poststratification (MRP) the weights
#'  are the population counts and the units are the unique combinations of
#'  all auxiliary variables used in the model, typically stored in a single \code{data.frame}.
# @param var variance(s) used for out-of-sample prediction. By default 1.
# @param ny number of trials used for out-of-sample prediction in case of a binomial model. By default 1.
# use fun. instead of fun to avoid argument name clash with parLapply
#' @param fun. function applied to the vector of posterior predictions to compute one or multiple summaries
#'  or test statistics. The function can have one or two arguments. The first argument is always the vector
#'  of posterior predictions. The optional second argument must be named 'p' and represents a list of model
#'  parameters, needed only when a test statistic depends on them. The function must return an integer or
#'  numeric vector.
#' @param labels optional names for the output object. Must be a vector of the same length as the result of \code{fun.}.
#' @param ppcheck if \code{TRUE}, function \code{fun.} is also applied to the observed data and
#'  an MCMC approximation is computed of the posterior predictive probability that the test statistic for
#'  predicted data is greater than the test statistic for the observed data.
#' @param iters iterations in \code{object} to use for prediction.
#'  Default \code{NULL} means that all draws from \code{object} are used.
#' @param to.file if \code{TRUE} the predictions are streamed to file.
#' @param filename name of the file to write predictions to in case \code{to.file=TRUE}.
#' @param write.single.prec Whether to write to file in single precision. Default is \code{FALSE}.
#' @param show.progress whether to show a progress bar.
#' @param verbose whether to show informative messages.
#' @param n.cores the number of cpu cores to use. Default is one, i.e. no parallel computation.
#'  If an existing cluster \code{cl} is provided, \code{n.cores} will be set to the number
#'  of workers in that cluster.
#' @param cl an existing cluster can be passed for parallel computation. If \code{NULL} and
#'  \code{n.cores > 1}, a new cluster is created.
#' @param seed a random seed (integer). For parallel computation it is used to independently
#'  seed RNG streams for all workers.
#' @param export a character vector with names of objects to export to the workers. This may
#'  be needed for parallel execution if expressions in \code{fun.} depend on global variables.
#' @param ... currently not used.
#' @returns An object of class \code{dc}, containing draws from the posterior (or prior) predictive distribution.
#'  If \code{ppcheck=TRUE} posterior predictive p-values are returned as an additional attribute.
#'  In case \code{to.file=TRUE} the file name used is returned.
# TODO: pp check only, i.e. without storing predictions either in an object or a file
predict.mcdraws <- function(object, newdata=NULL, X.=if (is.null(newdata)) "in-sample" else NULL,
                          type=c("data", "link", "response", "data_cat"),
                          weights=NULL,
                          #var=NULL,  # DEPRECATED
                          #ny=NULL, ry=NULL,  # DEPRECATED
                          fun.=identity, labels=NULL, ppcheck=FALSE, iters=NULL,
                          to.file=FALSE, filename, write.single.prec=FALSE,
                          show.progress=TRUE, verbose=TRUE,
                          n.cores=1L, cl=NULL, seed=NULL, export=NULL, ...) {
  type <- match.arg(type)
  if (ppcheck && type != "data") stop("posterior predictive checks only possible with type = 'data'")
  model <- object[["_model"]]
  fam <- model[["family"]]
  if (fam[["family"]] == "negbinomial") {
    if (!is.null(list(...)[["ry"]])) stop("argument 'ry' is no longer supported")
  }
  par.names <- par_names(object)
  if (is.null(iters))
    iters <- seq_len(object[["_info"]]$n.draw)
  else
    if (!all(iters %in% seq_len(object[["_info"]]$n.draw))) stop("non-existing iterations selected")
  n.chain <- n_chains(object)

  fun. <- match.fun(fun.)
  arg1 <- length(formals(args(fun.)))
  if (arg1 < 1L || arg1 > 2L) stop("'fun.' must be a function of 1 or 2 arguments")
  # arg1 is flag for single argument function, i.e. function of prediction x only, otherwise function of prediction x and state p
  arg1 <- arg1 == 1L || names(formals(args(fun.)))[2L] != "p"

  if (is.null(newdata)) {
    if (is.null(X.)) stop("one of 'newdata' and 'X.' must be supplied")
    if (identical(X., "in-sample")) {
      use.linpred_ <- any("linpred_" == par.names) && model[["do.linpred"]] && is.null(model[["linpred"]])
      if (!use.linpred_ && !all(names(Filter(\(mc) mc[["type"]] != "mc_offset", model[["mod"]])) %in% par.names))
        stop("for prediction all coefficients must be stored in 'object' (use 'store.all=TRUE' in MCMCsim)")
      X. <- list()
      for (mc in model[["mod"]]) X.[[mc[["name"]]]] <- mc$make_predict(verbose=verbose)
      if (fam[["family"]] == "multinomial" && type == "data_cat")
        n <- model[["n0"]]
      else
        n <- model[["n"]]
    } else {
      if (type == "data") {
        if (fam[["family"]] == "gamma") {
          # currently disallow vector shape parameter for X.="linpred" and custom X.
          if (!fam[["alpha.scalar"]]) stop("cannot derive vector shape")
        }
        if (any(fam[["family"]] == c("gaussian", "gaussian_gamma"))) {
          if (fam[["modeled.Q"]] || fam[["Q0.type"]] != "unit") stop("for prediction based on a model with non-trivial variance structure, please use argument 'newdata'")
        }
      }
      if (identical(X., "linpred")) {
        if (all("linpred_" != par.names))
          stop("'linpred_' not found in object. Use linpred='fitted' in create_sampler.")
        n <- n_vars(object[["linpred_"]])
        use.linpred_ <- TRUE
        X. <- NULL
      } else {
        if (!is.list(X.) || length(X.) == 0L) stop("'X.' must be a non-empty list")
        if (!all(names(X.) %in% names(model[["mod"]])))
          stop("X. names not corresponding to a model component: ", paste0(setdiff(names(X.), names(model[["mod"]])), collapse=", "))
        n <- nrow(X.[[1L]])
        if (!allv(i_apply(X., NROW), n)) stop("not all matrices in 'X.' have the same number of rows")
        for (k in names(X.)) {
          X.[[k]] <- model$mod[[k]]$make_predict(Xnew=X.[[k]], verbose=verbose)
        }
        use.linpred_ <- FALSE
      }
    }
    if (!use.linpred_) {
      # check that for mec components name_X is stored
      if (is.null(X.))
        name_Xs <- unlst(lapply(model[["mod"]], `[[`, "name_X"))
      else
        name_Xs <- unlst(lapply(model[["mod"]][names(X.)], `[[`, "name_X"))
      if (!all(name_Xs %in% par.names))
        stop("(in-sample) prediction for a model with measurement error components requires ",
          "the samples of the corresponding covariates; use 'store.all=TRUE' in MCMCsim"
        )
    }
  } else {
    if (!is.null(X.)) warn("argument 'X.' is ignored")
    use.linpred_ <- FALSE
    X. <- list()
    if (verbose && nrow(newdata) > 5e4L) message("setting up model design matrices for 'newdata'")
    for (mc in model[["mod"]]) X.[[mc$name]] <- mc$make_predict(newdata, verbose=verbose)
    if (!is.null(X.) && !all(names(Filter(\(mc) mc[["type"]] != "mc_offset", model[["mod"]])) %in% par.names))
      stop("for prediction all coefficients must be stored in 'object' (use 'store.all=TRUE' in MCMCsim)")
    if (fam[["family"]] == "multinomial")
      n <- nrow(newdata) * model[["Km1"]]
    else
      n <- nrow(newdata)
  }

  # determine dimension and mode (integer/double) of test statistic
  if ((type == "data" && any(fam[["family"]] == c("binomial", "negbinomial", "poisson", "multinomial"))) || type == "data_cat") {
    # x argument of fun. is integer
    if (arg1)
      temp <- fun.(rep.int(1L, n))
    else
      temp <- fun.(rep.int(1L, n), p=get_draw(object, 1L, 1L))
  } else {
    # x argument of fun. is double
    if (arg1)
      temp <- fun.(rep.int(0.1, n))
    else
      temp <- fun.(rep.int(0.1, n), p=get_draw(object, 1L, 1L))
  }
  d <- length(temp)
  out_type <- storage.mode(temp)
  if (all(out_type != c("logical", "integer", "double"))) stop("'fun.' must return an integer or numeric vector")
  rm(temp)

  if (!is.null(labels) && length(labels) != d) stop("incompatible 'labels' vector")
  if (fam[["family"]] == "multinomial") {
    if (is.null(labels) && isTRUE(all.equal(fun., identity)) && type != "data_cat")
      labels <- paste(rep_each(model$cats[-length(model$cats)], n %/% model[["Km1"]]), rep.int(seq_len(n %/% model[["Km1"]]), model[["Km1"]]), sep="_")
    if (type == "response") {  # undo stick-breaking transformation
      n0 <- n %/% model[["Km1"]]
      ind <- seq_len(n0)
      linkinv <- function(eta) {
        p <- fam$linkinv(eta)
        pleft <- 1 - p[ind]
        for (k in seq_len(model[["Km1"]] - 1L)) {
          ind <- ind + n0
          p[ind] <- p[ind] * pleft
          pleft <- pleft - p[ind]
        }
        p
      }
    }
  } else {
    if (type == "response") linkinv <- fam[["linkinv"]]
  }

  if (type == "data" || type == "data_cat") {
    if (!is.null(weights)) {
      if (!is.null(newdata) && inherits(weights, "formula")) {
        weights <- get_var_from_formula(weights, newdata)
      } else {
        if (all(length(weights) != c(1L, n))) stop("incompatible length of weights")
      }
    }
    if (type == "data_cat")
      rpredictive <- fam$make_rpredictive_cat(if (is.null(newdata)) n else newdata, weights)
    else
      rpredictive <- fam$make_rpredictive(if (is.null(newdata)) n else newdata, weights)
  } else {
    if (!is.null(weights)) warn("argument 'weights' is ignored when type is 'link' or 'response'")
  }

  if (is.null(cl))
    n.cores <- max(as.integer(n.cores)[1L], 1L)
  else
    n.cores <- length(cl)
  if (n.cores > 1L) {
    if (n.cores > parallel::detectCores()) stop("too many cores")
    if (is.null(cl)) {
      cl <- setup_cluster(n.cores, seed)
      on.exit(stop_cluster(cl))
    } else {
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    }
    if (!is.null(export)) {
      if (!is.character(export)) stop("'export' must be a character vector")
      parallel::clusterExport(cl, export, parent.frame())
    }
    if (verbose) message(length(iters) * n.chain, " draws distributed over ", n.cores, " cores")
    sim_list <- split_iters(object, iters, parts=n.cores)
    predict_obj <- function(obj) {
      chains <- seq_len(n_chains(obj))
      iters <- seq_len(n_draws(obj))
      if (out_type == "logical" || out_type == "integer")
        out <- rep.int(list(matrix(NA_integer_, length(iters), d)), length(chains))
      else
        out <- rep.int(list(matrix(NA_real_, length(iters), d)), length(chains))
      if (ppcheck) {
        ppp <- numeric(d)
        if (arg1)
          fy <- fun.(model[["y"]])
        else
          y <- model[["y"]]
      }
      for (i in iters) {
        for (ch in chains) {
          p <- get_draw(obj, i, ch)
          if (use.linpred_)
            ystar <- p[["linpred_"]]
          else
            ystar <- model$lin_predict(p, X.)
          switch(type,
            data=, data_cat = ystar <- rpredictive(p, ystar),
            response = ystar <- linkinv(ystar)
          )
          out[[ch]][i, ] <- if (arg1) fun.(ystar) else fun.(ystar, p=p)
          if (ppcheck) ppp <- ppp + (out[[ch]][i, ] >= if (arg1) fy else fun.(y, p=p))
        }
      }
      if (ppcheck) attr(out, "ppp") <- ppp
      out
    }
    out <- combine_iters_dc(parallel::parLapply(cl, X=sim_list, fun=predict_obj))
  } else {
    if (!is.null(seed)) set.seed(seed)
    n.it <- length(iters)
    if (to.file) {
      if (missing(filename)) filename <- "MCdraws_pred.dat"
      outfile <- file(filename, "wb")
      write_header(outfile, n.it, n.chain, d, labels, write.single.prec)
      write.size <- if (write.single.prec) 4L else NA_integer_
      ppcheck <- FALSE  # TODO allow ppcheck in the case that predictions are written to file
    } else {
      if (out_type == "logical" || out_type == "integer")
        out <- rep.int(list(matrix(NA_integer_, n.it, d)), n.chain)
      else
        out <- rep.int(list(matrix(NA_real_, n.it, d)), n.chain)
    }
    if (ppcheck) {
      ppp <- numeric(d)
      if (arg1)
        fy <- fun.(model[["y"]])
      else
        y <- model[["y"]]
    }
    r <- 1L
    show.progress <- show.progress && n.it > 1L
    if (show.progress) pb <- txtProgressBar(min=1L, max=n.it, style=3L)
    for (i in iters) {
      for (ch in seq_len(n.chain)) {
        p <- get_draw(object, i, ch)
        if (use.linpred_)
          ystar <- p[["linpred_"]]
        else
          ystar <- model$lin_predict(p, X.)
        switch(type,
          data=, data_cat = ystar <- rpredictive(p, ystar),
          response = ystar <- linkinv(ystar)
        )
        if (to.file)
          writeBin(if (arg1) fun.(ystar) else fun.(ystar, p=p), con=outfile, size=write.size)
        else
          out[[ch]][r, ] <- if (arg1) fun.(ystar) else fun.(ystar, p=p)
        if (ppcheck) ppp <- ppp + (out[[ch]][r, ] >= if (arg1) fy else fun.(y, p=p))
      }
      if (show.progress) setTxtProgressBar(pb, r)
      r <- r + 1L
    }
    if (show.progress) close(pb)
    if (to.file) {
      close(outfile)
      return(filename)
    }
    if (ppcheck) attr(out, "ppp") <- ppp/(n.chain * length(iters))
  }
  if (!is.null(labels)) attr(out, "labels") <- labels
  class(out) <- "dc"
  out
}

#' Generate a data vector according to a model
#'
#' This function generates draws from the prior predictive distribution.
#' Parameter values are drawn from their priors, and consequently data is
#' generated from the sampling distribution given these parameter values.
#'
#' @examples
#' \donttest{
#' n <- 250
#' dat <- data.frame(
#'   x = rnorm(n),
#'   g = factor(sample(1:10, n, replace=TRUE)),
#'   ny = 10
#' )
#' gd <- generate_data(
#'   ~ reg(~ 1 + x, prior=pr_normal(precision=10, mean=c(0, 1)), name="beta") +
#'     gen(factor = ~ g, name="v"),
#'   family=f_binomial(n.trial = ~ ny), data=dat
#' )
#' gd
#' plot(dat$x, gd$y)
#' }
#'
#' @export
#' @param formula A model formula, see \code{\link{create_sampler}}.
#'  Any left-hand side of the formula is ignored.
#' @param data see \code{\link{create_sampler}}.
#' @param family sampling distribution family, see \code{\link{create_sampler}}.
#' @param ny NO LONGER USED; see \code{\link{create_sampler}}.
#' @param ry NO LONGER USED; see \code{\link{create_sampler}}.
#' @param r.mod NO LONGER USED; see \code{\link{create_sampler}}.
#' @param sigma.fixed DEPRECATED; see \code{\link{create_sampler}}.
#' @param sigma.mod DEPRECATED; see \code{\link{create_sampler}}.
#' @param Q0 DEPRECATED; see \code{\link{create_sampler}}.
#' @param formula.V DEPRECATED; see \code{\link{create_sampler}}.
#' @param linpred see \code{\link{create_sampler}}.
#' @returns A list with a generated data vector and a list of prior means of the
#'  parameters. The parameters are drawn from their priors.
generate_data <- function(formula, data=NULL, family="gaussian",
                          ny=NULL, ry=NULL, r.mod=NULL,  # DEPRECATED
                          sigma.fixed=NULL, sigma.mod=NULL, Q0=NULL, formula.V=NULL,  # DEPRECATED
                          linpred=NULL) {
  if (is.function(family)) {
    family <- as.character(substitute(family))
    if (startsWith(family, "f_")) family <- substring(family, 3L)
  }
  if (identical(family, "gaussian")) {
    # use a proper prior for sigma by default
    family <- f_gaussian(var.prior = pr_invchisq(df=1, scale=1))
  }
  sampler <- create_sampler(formula=formula, data=data, family=family,
    ny=ny, ry=ry, r.mod=r.mod,  # DEPRECATED
    sigma.fixed=sigma.fixed, sigma.mod=sigma.mod, Q0=Q0, formula.V=formula.V,  # DEPRECATED
    linpred=linpred, prior.only=TRUE
  )
  sim <- MCMCsim(sampler, n.iter=1L, n.chain=1L, store.all=TRUE, from.prior=TRUE, verbose=FALSE)
  pars <- lapply(sim[par_names(sim)], \(x) setNames(x[[1L]][1L, ], attr(x, "labels")))
  if (sampler$family[["family"]] == "multinomial") {
    if (allv(sampler$family[["ny0"]], 1)) {
      y <- drop(predict(sim, type="data_cat")[[1L]])
    } else {
      y <- matrix(predict(sim)[[1L]], sampler[["n0"]], sampler[["Km1"]])
      y <- cbind(y, sampler$family[["ny0"]] - dapply(y, sum, MARGIN=1L))  # NB rowSums does not preserve integer type
    }
  } else {
    y <- drop(predict(sim)[[1L]])
  }
  list(y=y, pars=pars)
}
