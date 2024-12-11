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
#'  For binomial and negative binomial models \code{type="response"} returns the simulations
#'  of the latent probabilities. For multinomial models \code{type="link"} generates
#'  the linear predictor for all categories except the last, and \code{type="response"} transforms
#'  this vector to the probability scale, and \code{type="data"} generates the multinomial data,
#'  all in long vector format, where the output for all categories (except the last) are stacked.
#'  For multinomial models and single trials, a further option is \code{type="data_cat"},
#'  which generates the data as a categorical vector, with integer coded levels.
#' @param var variance(s) used for out-of-sample prediction. By default 1.
#' @param ny number of trials used for out-of-sample prediction in case of a binomial model. By default 1.
#' @param ry fixed part of the (reciprocal) dispersion parameter in case of a negative binomial model.
# use fun. instead of fun to avoid argument name clash with parLapply
#' @param fun. function applied to the vector of posterior predictions to compute one or multiple summaries
#'  or test statistics. The function can have one or two arguments. The first argument is always the vector
#'  of posterior predictions. The optional second argument represents a list of model parameters, needed only
#'  when a test statistic depends on them. The function must return an integer or numeric vector.
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
                          var=NULL, ny=NULL, ry=NULL,
                          fun.=identity, labels=NULL, ppcheck=FALSE, iters=NULL,
                          to.file=FALSE, filename, write.single.prec=FALSE,
                          show.progress=TRUE, verbose=TRUE,
                          n.cores=1L, cl=NULL, seed=NULL, export=NULL, ...) {
  type <- match.arg(type)
  if (ppcheck && type != "data") stop("posterior predictive checks only possible with type = 'data'")
  model <- object[["_model"]]
  fam <- model$family
  par.names <- par_names(object)
  if (is.null(iters))
    iters <- seq_len(object[["_info"]]$n.draw)
  else
    if (!all(iters %in% seq_len(object[["_info"]]$n.draw))) stop("non-existing iterations selected")
  n.chain <- n_chains(object)

  fun. <- match.fun(fun.)
  arg1 <- length(formals(args(fun.)))
  if (arg1 < 1L || arg1 > 2L) stop("'fun.' must be a function of 1 or 2 arguments")
  arg1 <- arg1 == 1L  # flag for single argument function, i.e. function of prediction x only

  cholQ <- NULL
  V <- NULL
  if (is.null(newdata)) {
    if (is.null(X.)) stop("one of 'newdata' and 'X.' must be supplied")
    if (identical(X., "in-sample")) {
      use.linpred_ <- any("linpred_" == par.names) && model$do.linpred && is.null(model$linpred)
      if (!use.linpred_ && !all(names(Filter(function(mc) mc[["type"]] != "mc_offset", model[["mod"]])) %in% par.names))
        stop("for prediction all coefficients must be stored in 'object' (use 'store.all=TRUE' in MCMCsim)")
      X. <- list()
      for (mc in model$mod) X.[[mc$name]] <- mc$make_predict(verbose=verbose)
      n <- model[["n"]]
      if (!is.null(var)) {
        warn("argument 'var' ignored for in-sample prediction")
        var <- NULL
      }
      if (model$modeled.Q) {
        cholQ <- build_chol(
          switch(model$Q0.type,
            unit = rep.int(1, model$n),
            diag = model$Q0@x,
            symm = model$Q0
          )
        )
      } else {
        cholQ <- build_chol(model$Q0)
      }
      if (!is.null(ny)) {
        warn("argument 'ny' ignored for in-sample prediction")
        ny <- NULL
      }
      # BayesLogit::rpg flags ny=0, so we have set ny to a tiny value > 0 in sampler; need to undo it here as rbinom yields NA for non-integral ny
      if (fam$family == "binomial" || fam$family == "multinomial") ny <- as.integer(round(model$ny))
      if (fam$family == "multinomial" && type == "data_cat") n <- n %/% model$Km1
      if (fam$family == "negbinomial") ry <- model$ry
    } else {
      if (fam$family == "gamma") {
        # currently disallow vector shape parameter for X.="linpred" and custom X.
        if (!fam$alpha.scalar) stop("cannot derive vector shape")
      }
      if (identical(X., "linpred")) {
        if (all("linpred_" != par.names))
          stop("'linpred_' not found in object. Use linpred='fitted' in create_sampler.")
        n <- n_vars(object[["linpred_"]])
        use.linpred_ <- TRUE
        X. <- NULL
      } else {
        if (!is.list(X.) || length(X.) == 0L) stop("'X.' must be a non-empty list")
        if (!all(names(X.) %in% names(model$mod)))
          stop("X. names not corresponding to a model component: ", paste0(setdiff(names(X.), names(model$mod)), collapse=", "))
        n <- nrow(X.[[1L]])
        if (!all(i_apply(X., NROW) == n)) stop("not all matrices in 'X.' have the same number of rows")
        for (k in names(X.)) {
          X.[[k]] <- model$mod[[k]]$make_predict(Xnew=X.[[k]], verbose=verbose)
        }
        use.linpred_ <- FALSE
      }
      # by default use unit variances for new (oos) predictions, but warn if in-sample variance matrix is not unit diagonal
      if (is.null(var)) {
        if (type == "data" && model$Q0.type != "unit") warn("no 'var' specified; unit variances assumed for prediction")
        var <- 1
      }
      if (all(length(var) != c(1L, n))) stop("'var' should be either a scalar value or a vector of length 'nrow(newdata)'")
      if (model$modeled.Q && fam$family == "gaussian") stop("please use argument 'newdata' instead of 'X.' to also supply information about the factors determining the (modeled) sampling variances")
      if (any(fam$family == c("binomial", "multinomial")) && type == "data") {
        if (is.null(ny)) ny <- model$ny.input
        if (is.character(ny)) stop("please supply number of binomial trials 'ny' as a vector")
        ny <- check_ny(ny, n)
      } else if (fam$family == "negbinomial" && type == "data") {
        if (is.null(ry)) ry <- model$ry.input
        if (is.character(ry)) stop("please supply 'ry' as a vector")
        ry <- check_ry(ry, n)
      }
    }
    if (!use.linpred_) {
      # check that for mec components name_X is stored
      name_Xs <- unlst(lapply(model$mod, `[[`, "name_X"))
      if (!is.null(X.)) name_Xs <- name_Xs[names(X.)]
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
    for (mc in model$mod) X.[[mc$name]] <- mc$make_predict(newdata, verbose=verbose)
    if (!is.null(X.) && !all(names(Filter(function(mc) mc[["type"]] != "mc_offset", model[["mod"]])) %in% par.names))
      stop("for prediction all coefficients must be stored in 'object' (use 'store.all=TRUE' in MCMCsim)")
    if (fam$family == "multinomial")
      n <- nrow(newdata) * model$Km1
    else
      n <- nrow(newdata)
    if (is.null(var)) {
      if (type == "data" && model$Q0.type != "unit") warn("no 'var' specified, so unit variances are assumed for prediction")
      var <- 1
    }
    if (all(length(var) != c(1L, n))) stop("'var' should be either a scalar value or a vector of length 'nrow(newdata)'")
    if (model$modeled.Q && fam$family == "gaussian") {
      V <- list()
      for (Vmc in model$Vmod) V[[Vmc$name]] <- Vmc$make_predict_Vfactor(newdata)
    }
    if (any(fam$family == c("binomial", "multinomial")) && any(type == c("data", "data_cat"))) {
      if (is.null(ny) && fam$family == "binomial") ny <- model$ny.input
      ny <- check_ny(ny, newdata)
    } else if (fam$family == "negbinomial" && type == "data") {
      if (is.null(ry)) ry <- model$ry.input
      ry <- check_ry(ry, newdata)
    }
  }

  # determine dimension and mode (integer/double) of test statistic
  if ((type == "data" && fam$family != "gaussian") || type == "data_cat") {
    # x argument of fun. is integer
    if (arg1) {
      temp <- fun.(rep.int(1L, n))
    } else {
      temp <- fun.(rep.int(1L, n), get_draw(object, 1L, 1L))
    }
  } else {
    # x argument of fun. is double
    if (arg1) {
      temp <- fun.(rep.int(0.1, n))  
    } else {
      temp <- fun.(rep.int(0.1, n), get_draw(object, 1L, 1L))
    }
  }
  d <- length(temp)
  out_type <- storage.mode(temp)
  if (all(out_type != c("logical", "integer", "double"))) stop("'fun.' must return an integer or numeric vector")
  rm(temp)

  if (!is.null(labels) && length(labels) != d) stop("incompatible 'labels' vector")
  if (fam$family == "multinomial") {
    if (is.null(labels) && isTRUE(all.equal(fun., identity)) && type != "data_cat")
      labels <- paste(rep_each(model$cats[-length(model$cats)], n %/% model$Km1), rep.int(seq_len(n %/% model$Km1), model$Km1), sep="_")
    if (type == "response") {  # undo stick-breaking transformation
      linkinv <- function(eta) {
        p <- fam$linkinv(eta)
        ind <- seq_len(model$n0)
        pleft <- 1 - p[ind]
        for (k in 2:model$Km1) {
          ind <- ind + model$n0
          p[ind] <- p[ind] * pleft
          pleft <- pleft - p[ind]
        }
        p
      }
    }
  } else {
    if (type == "response") linkinv <- fam$linkinv
  }

  if (type %in% c("data", "data_cat")) {
    if (any(fam$family == c("gamma", "poisson"))) {
      rpred <- fam$make_rpredictive(newdata)
      rpredictive <- function(p, lp, cholQ=NULL, var=NULL, V=NULL, ny, ry)
        rpred(p, lp)
    } else
      rpredictive <- if (type == "data") model$rpredictive else model$rpredictive_cat
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
          fy <- fun.(model$y)
        else
          y <- model$y
      }
      for (i in iters) {
        for (ch in chains) {
          p <- get_draw(obj, i, ch)
          if (use.linpred_)
            ystar <- p[["linpred_"]]
          else
            ystar <- model$lin_predict(p, X.)
          switch(type,
            data=, data_cat = ystar <- rpredictive(p, ystar, cholQ, var, V, ny, ry),
            response = ystar <- linkinv(ystar)
          )
          out[[ch]][i, ] <- if (arg1) fun.(ystar) else fun.(ystar, p)
          if (ppcheck) ppp <- ppp + (out[[ch]][i, ] >= if (arg1) fy else fun.(y, p))
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
        fy <- fun.(model$y)
      else
        y <- model$y
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
          data=, data_cat = ystar <- rpredictive(p, ystar, cholQ, var, V, ny, ry),
          response = ystar <- linkinv(ystar)
        )
        if (to.file)
          writeBin(if (arg1) fun.(ystar) else fun.(ystar, p), con=outfile, size=write.size)
        else
          out[[ch]][r, ] <- if (arg1) fun.(ystar) else fun.(ystar, p)
        if (ppcheck) ppp <- ppp + (out[[ch]][r, ] >= if (arg1) fy else fun.(y, p))
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
#'   family="binomial", ny="ny", data=dat
#' )
#' gd
#' plot(dat$x, gd$y)
#' }
#'
#' @export
#' @param formula A model formula, see \code{\link{create_sampler}}.
#'  Any left-hand-side of the formula is ignored.
#' @param data see \code{\link{create_sampler}}.
#' @param family sampling distribution family, see \code{\link{create_sampler}}.
#' @param ny see \code{\link{create_sampler}}.
#' @param ry see \code{\link{create_sampler}}.
#' @param r.mod see \code{\link{create_sampler}}.
#' @param sigma.fixed see \code{\link{create_sampler}}.
#' @param sigma.mod see \code{\link{create_sampler}}.
#' @param Q0 see \code{\link{create_sampler}}.
#' @param formula.V see \code{\link{create_sampler}}.
#' @param linpred see \code{\link{create_sampler}}.
#' @returns A list with a generated data vector and a list of prior means of the
#'  parameters. The parameters are drawn from their priors.
generate_data <- function(formula, data=NULL, family="gaussian",
                          ny=NULL, ry=NULL, r.mod,
                          sigma.fixed=NULL, sigma.mod=NULL, Q0=NULL, formula.V=NULL,
                          linpred=NULL) {
  if (identical(family, "gaussian") || (is.environment(family) && family$family == "gaussian")) {
    if (is.null(sigma.fixed)) sigma.fixed <- FALSE
    if (!sigma.fixed && is.null(sigma.mod)) {
      # use a proper prior for sigma by default
      sigma.mod <- pr_invchisq(df=1, scale=1)
      sigma.mod$init(n=1L)
    }
  }
  sampler <- create_sampler(formula=formula, data=data, family=family, ny=ny, ry=ry, r.mod=r.mod,
                            sigma.fixed=sigma.fixed, sigma.mod=sigma.mod, Q0=Q0, formula.V=formula.V,
                            linpred=linpred, prior.only=TRUE)
  sim <- MCMCsim(sampler, n.iter=1L, n.chain=1L, store.all=TRUE, from.prior=TRUE, verbose=FALSE)
  pars <- lapply(sim[par_names(sim)], function(x) setNames(x[[1L]][1L, ], attr(x, "labels")))
  if (sampler$family$family == "multinomial") {
    if (is.null(ny) || all(ny == 1)) {
      y <- as.vector(predict(sim, type="data_cat")[[1L]])
    } else {
      y <- matrix(predict(sim)[[1L]], sampler$n0, sampler$Km1)
      y <- cbind(y, sampler$ny - apply(y, 1L, sum))
    }
  } else {
    y <- as.vector(predict(sim)[[1L]])
  }
  list(y=y, pars=pars)
}
