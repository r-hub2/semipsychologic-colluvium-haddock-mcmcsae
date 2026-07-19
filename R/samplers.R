#' Create a sampler object
#'
#' This function sets up a sampler object, based on the specification of a model. The object contains functions to
#' draw a set of model parameters from their prior and conditional posterior distributions, and
#' to generate starting values for the MCMC simulation. The functions share a common environment
#' containing precomputed quantities such as design matrices based on the model and the data.
#' The sampler object is the main input for the MCMC simulation function \code{\link{MCMCsim}}.
#'
#' The right hand side of the \code{formula} argument to \code{create_sampler} can be used to specify
#' additive model components. Currently four model components are supported: \code{\link{reg}(...)}
#' for regression or 'fixed' effects, \code{\link{gen}(...)} for generic random effects,
#' \code{\link{mec}(...)} for measurement in covariates effects, and \code{\link{brt}(...)}
#' for a Bayesian additive regression trees component. Note that an offset can be added
#' separately, in the usual way using \code{\link{offset}(...)}. As usual, the left hand side of
#' \code{formula} specifies the response variable. For binomial and multinomial models see
#' their respective help pages (\code{\link{f_binomial}}, \code{\link{f_multinomial}}) for
#' different ways to specify the response variable. If no response variable is specified,
#' \code{create_sampler} only sets up infrastructure to sample from the model parameters'
#' prior distributions.
#'
#' To specify family-specific parameters or options, see the help pages for the various
#' family specification functions, e.g. \code{\link{f_gaussian}} or \code{\link{f_negbinomial}}.
#'
#' @examples
#' # first generate some data
#' n <- 200
#' x <- rnorm(n)
#' y <- 0.5 + 2*x + 0.3*rnorm(n)
#' # create a sampler for a simple linear regression model
#' sampler <- create_sampler(y ~ x)
#' sim <- MCMCsim(sampler)
#' (summary(sim))
#'
#' y <- rbinom(n, 1, 1 / (1 + exp(-(0.5 + 2*x))))
#' # create a sampler for a binary logistic regression model
#' sampler <- create_sampler(y ~ x, family="binomial")
#' sim <- MCMCsim(sampler)
#' (summary(sim))
#'
#' @export
#' @param formula formula to specify the response variable and additive model components. The model components
#'  form the linear predictor part of the model. A model component on the right hand side can be either
#'  a regression term specified by \code{\link{reg}(...)}, a covariates subject to error term specified
#'  by \code{\link{mec}(...)}, or a generic random effect term specified by \code{\link{gen}(...)}.
#'  See for details the help pages for these model component creation functions.
#'  An offset can be specified as \code{offset(...)}.
#'  Other terms in the formula are collectively interpreted as ordinary regression effects,
#'  treated in the same way as a \code{reg(...)} term, but without the option to change the prior.
#' @param data data frame with n rows in which the variables specified in model components can be found.
#' @param family character string describing the data distribution. The default is 'gaussian'.
#'  Other options are 'binomial', 'multinomial', 'negbinomial' for the negative binomial distribution,
#'  'poisson', and 'gamma'. Alternatively, functions starting with 'f_' followed
#'  by the family name can be used to specify the sampling distribution and possibly
#'  further options. See \code{\link{f_gaussian}}, \code{\link{f_binomial}},
#'  \code{\link{f_multinomial}}, \code{\link{f_negbinomial}}, \code{\link{f_poisson}},
#'  \code{\link{f_gamma}} and \code{\link{f_gaussian_gamma}}.
#'  For categorical or multinomial data, use \code{family="multinomial"} or \code{family=f_multinomial()}
#'  where the second form allows to additionally specify the number of trials (if not equal to 1).
#'  A stick-breaking representation of the multinomial distribution is used for model fitting,
#'  and the logistic link function relates each category except the last to a linear predictor.
#'  The categories can be referenced in the model specification formula by 'cat_'.
#' @param ny NO LONGER USED. Please use \code{\link{f_binomial}} to specify the numbers of trials.
#' @param ry NO LONGER USED. Please use \code{\link{f_negbinomial}} to specify
#'  further options for the negative binomial sampling distribution.
#' @param r.mod NO LONGER USED. Please use \code{\link{f_negbinomial}} to specify
#'  further options for the negative binomial sampling distribution.
#' @param sigma.fixed for Gaussian models, if \code{TRUE} the residual standard deviation parameter 'sigma_' is fixed at 1. In that case
#'  argument \code{sigma.mod} is ignored. This is convenient for Fay-Herriot type models with (sampling) variances assumed to be known.
#'  Default is \code{FALSE}. DEPRECATED, please use \code{\link{f_gaussian}} to specify
#'  variance options. In particular, a fixed scalar variance parameter with value 1 can be
#'  specified with \code{var.prior=pr_fixed(value=1)}.
#' @param sigma.mod prior for the variance parameter of a gaussian sampling distribution.
#'  This can be specified by a call to one of the prior specification functions
#'  \code{\link{pr_invchisq}}, \code{\link{pr_exp}}, \code{\link{pr_gig}} or \code{\link{pr_fixed}} for
#'  inverse chi-squared, exponential, generalised inverse gaussian or degenerate prior distribution,
#'  respectively. The default is an improper prior \code{pr_invchisq(df=0, scale=1)}. A half-t prior on the
#'  standard deviation can be specified using \code{\link{pr_invchisq}} with a chi-squared distributed scale
#'  parameter. DEPRECATED, please use \code{\link{f_gaussian}} to specify variance options.
#'  In particular, to change the default prior for the scalar variance parameter use argument \code{var.prior}.
#' @param Q0 n x n data-level precision matrix for a Gaussian model. It defaults to the unit matrix.
#'  If an n-vector is provided it will be expanded to a (sparse) diagonal matrix with Q0 on its diagonal.
#'  If a name is supplied it will be looked up in \code{data} and subsequently expanded to a diagonal matrix.
#'  DEPRECATED, please use \code{\link{f_gaussian}}, in particular its \code{precision} argument,
#'  to specify unequal variances, or a non-diagonal precision matrix.
#' @param formula.V a formula specifying the terms of a variance model in the case of a Gaussian likelihood.
#'  Currently two types of terms are supported: a regression term for the log-variance
#'  specified with \code{\link{vreg}(...)}, and a term \code{\link{vfac}(...)} for multiplicative modelled factors
#'  at a certain level specified by a factor variable. By using unit-level inverse-chi-squared factors the marginal
#'  sampling distribution becomes a Student-t distribution, and by using unit-level exponential factors it becomes
#'  a Laplace or double exponential distribution. DEPRECATED, please use \code{\link{f_gaussian}},
#'  in particular its \code{var.model} argument, to specify a variance model.
#' @param logJacobian if the data are transformed the logarithm of the Jacobian can be supplied so that it
#'  is incorporated in all log-likelihood computations. This can be useful for comparing information criteria
#'  for different transformations. It should be supplied as a vector of the same size as the response variable.
#'  For example, when a log-transformation is used on response vector \code{y}, the vector \code{-log(y)}
#'  should be supplied. DEPRECATED, this argument has moved to \code{\link{f_gaussian}}.
#' @param linpred a list of matrices defining (possibly out-of-sample) linear predictors to be simulated.
#'  This allows inference on e.g. (sub)population totals or means. The list must be of the form
#'  \code{list(name_1=X_1, ...)} where the names refer to the model component names and predictions are
#'  computed by summing \code{X_i \%*\% p[[name_i]]}. Alternatively, \code{linpred="fitted"} can be used
#'  as a short-cut for simulations of the full in-sample linear predictor.
#' @param compute.weights if \code{TRUE} weights are computed for each element of \code{linpred}. Note that for
#'  a large dataset in combination with vector-valued linear predictors the weights can take up a lot of memory.
#'  By default only means are stored in the simulation carried out using \code{\link{MCMCsim}}.
#' @param block DEPRECATED, please use argument \code{control} instead, see also \code{\link{sampler_control}}.
#'  Note that this parameter is now by default set to \code{TRUE}.
#' @param prior.only whether a sampler is set up only for sampling from the prior or for sampling from both prior
#'  and posterior distributions. Default \code{FALSE}. If \code{TRUE} there is no need to specify a response in
#'  \code{formula}. This is used by \code{\link{generate_data}}, which samples from the prior predictive
#'  distribution.
#' @param control a list with further computational options. These options can
#'  be specified using function \code{\link{sampler_control}}.
#' @returns A sampler object, which is the main input for the MCMC simulation
#'  function \code{\link{MCMCsim}}. The sampler object is an environment with
#'  precomputed quantities and functions. The main functions are \code{rprior},
#'  which returns a sample from the prior distributions, \code{draw},
#'  which returns a sample from the full conditional posterior distributions,
#'  and \code{start}, which returns a list with starting values for the Gibbs
#'  sampler. If \code{prior.only} is \code{TRUE}, functions \code{draw} and
#'  \code{start} are not created.
#' @references
#'  J.H. Albert and S. Chib (1993).
#'    Bayesian analysis of binary and polychotomous response data.
#'    Journal of the American statistical Association 88(422), 669-679.
#'
#'  D. Bates, M. Maechler, B. Bolker and S.C. Walker (2015).
#'    Fitting Linear Mixed-Effects Models Using lme4.
#'    Journal of Statistical Software 67(1), 1-48.
#'
#'  S.W. Linderman, M.J. Johnson and R.P. Adams (2015).
#'    Dependent multinomial models made easy: Stick-breaking with the Polya-Gamma augmentation.
#'    Advances in Neural Information Processing Systems, 3456-3464.
#'
#'  J.S. Liu and Y.N. Wu (1999).
#'    Parameter expansion for data augmentation.
#'    Journal of the American Statistical Association, 94(448), 1264-1274.
#'
#'  P.A. Parker, S.H. Holan and R. Janicki (2024).
#'    Conjugate Modeling Approaches for Small Area Estimation with Heteroscedastic Structure.
#'    Journal of Survey Statistics and Methodology 12(4), 1061-1080.
#'
#'  N. Polson, J.G. Scott and J. Windle (2013).
#'    Bayesian Inference for Logistic Models Using Polya-Gamma Latent Variables.
#'    Journal of the American Statistical Association 108(504), 1339-1349.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
create_sampler <- function(formula, data=NULL, family="gaussian",
                           ny=NULL, ry=NULL, r.mod=NULL,  # DEPRECATED
                           sigma.fixed=NULL, sigma.mod=NULL, Q0=NULL, formula.V=NULL, logJacobian=NULL,  # DEPRECATED
                           linpred=NULL,
                           compute.weights=FALSE, block=NULL,
                           prior.only=FALSE,
                           control=sampler_control()) {

  if (!is.null(data) && !inherits(data, "data.frame")) {
    if (is.numeric(data) && length(data) == 1L) {
      # a single numeric value is interpreted as sample size
      data <- as.integer(data)
      if (data < 0L) stop("negative data size")
    } else if (!inherits(data, "data.frame")) data <- as.data.frame(data)
  }

  family <- process_family(family)

  # deprecation code
  if (!(is.null(sigma.fixed) && is.null(sigma.mod) && is.null(Q0) && is.null(formula.V) && is.null(logJacobian))) {
    if (family[["family"]] == "gaussian") {
      warn("arguments 'sigma.fixed', 'sigma.mod', 'Q0', 'formula.V' and 'logJacobian' are deprecated\n",
           "please pass variance options to argument 'family' using function f_gaussian\n")
      # TODO if non-default options have been set using f_gaussian raise an error
      family <- f_gaussian(
        var.prior = if (isTRUE(sigma.fixed)) pr_fixed(1) else if (!is.null(sigma.mod)) sigma.mod else pr_invchisq(0, 1),
        var.vec =
          if (is.null(Q0) || is_a_matrix(Q0))
            ~ 1
          else if (is.character(Q0))
            as.formula(paste0("~ I(1/", Q0, ")"))
          else if (is.vector(Q0))
            1/Q0
          else
            as.formula(paste0("~ I(1/(", deparse(substitute(Q0)), "))")),
        prec.mat = if (is_a_matrix(Q0)) Q0 else NULL,
        var.model = formula.V, logJacobian = logJacobian
      )
    } else {
      warn("arguments 'sigma.fixed', 'sigma.mod', 'Q0', 'formula.V' and 'logJacobian' are ignored for family '", family[["family"]], "'")
    }
  }
  if (family[["family"]] == "binomial" && !is.null(ny)) stop("argument 'ny' no longer supported; please specify the number of trials using function f_binomial")
  if (family[["family"]] == "negbinomial" && !is.null(ry) || !is.null(r.mod))
    stop("Arguments 'ry' and 'r.mod' are no longer supported. ",
         "Please use function f_negbinomial instead.")
  if (family[["family"]] == "poisson" && !is.null(ry)) stop("argument 'ry' is no longer supported; please use 'control' argument of f_poisson instead")

  if (missing(formula)) stop("a formula specifying response and linear predictor model components must be specified")
  if (!inherits(formula, "formula")) stop("'formula' must be a formula")
  # Poisson approximated by negative binomial with large shape parameter and (internal) offset = -log(shape)
  internal.offset <- family[["family"]] == "poisson" || (family[["family"]] == "multi" &&
    any(b_apply(family[["fam.list"]], function(f) f[["family"]] == "poisson")))
  formula <- standardise_formula(formula, data=data, internal.offset=internal.offset)
  mod <- to_mclist(formula)
  if (!length(mod)) stop("empty 'formula'")
  types <- get_types(mod)
  if (any(types %in% c("vreg", "vfac"))) stop("'vreg' and 'vfac' can only be used in variance model specification")
  if (family[["family"]] == "gamma" && !all(types %in% c("mc_offset", "reg", "gen"))) stop("only 'reg' and 'gen' (and offset/mc_offset) components supported for family 'gamma'")
  offset.only <- all(types == "mc_offset")
  has.bart <- any(types == "brt")

  control <- check_sampler_control(control)
  if (!prior.only) {
    # check Gibbs blocks for mean model specification
    if (!is.null(block)) {
      warn("'block' argument of create_sampler is deprecated; ",
           "instead, please use argument 'control' to specify the Gibbs sampling blocks; ",
           "by default a maximal blocking strategy is used so that coefficients are sampled together whenever possible")
      control$block <- block
    }
    if (is.logical(control[["block"]])) {
      if (control[["block"]]) {  # the default
        # all components of type reg, gen and mec in a single block
        control$block <- list(names(mod)[!(types %in% c("mc_offset", "brt"))])
        # single component by default not handled as a block, unless
        #   compute.weights=TRUE, or cMVN or CG sampler is used
        if (!length(control[["block"]][[1L]]) || (length(control[["block"]][[1L]]) == 1L &&
            types[!(types %in% c("mc_offset", "brt"))] != "s" && !compute.weights &&
            !control[["cMVN.or.CG"]] &&
            !(family[["link"]] == "probit" && family$control[["probit.HaarPXDA"]])))
          control$block <- list()
      } else {
        control$block <- list()
      }
    } else {
      for (bl in control[["block"]]) {
        if (!all(bl %in% names(mod)))
          stop("invalid name(s) '", paste0(setdiff(bl, names(mod)), collapse="', '"), "' in 'block'")
        if (any(types[bl] == "brt")) stop("'brt' model component cannot be part of a Gibbs block")
        if (any(types[bl] == "mc_offset")) stop("'offset/mc_offset' model component cannot be part of a Gibbs block")
      }
      rm(bl)
      if (any_duplicated(unlst(control[["block"]]))) stop("duplicate model component name(s) in 'block'")
    }
    if (family[["link"]] == "probit" && !length(control[["block"]])) family$control$probit.HaarPXDA <- FALSE
    control$single.block <- any(length(mod) == c(1L, length(unlst(control[["block"]]))))
    control$length1mod <- length(mod) == 1L
  }

  if (!prior.only && has_response(formula)) {
    y <- get_response(formula, data)
    n <- NROW(y)
    if (is.null(data))
      data <- n
    else if (n != nrow(data))
      stop("response vector has length ", n, " while data has ", nrow(data), " rows")
    if (control[["single.block"]]) {  # computation of working response may simplify
      control$recompute.e <- FALSE  # already computed in the coefficient sampling function
    }
  } else {
    if (!prior.only) {
      warn("no left hand side in 'formula': setting up prior samplers only")
      prior.only <- TRUE
    }
    if (is.null(data)) {
      vars <- all.vars(formula)
      vars <- vars[!(vars %in% c("local_", "global_"))]
      if (!length(vars)) stop("no way to determine the number of cases")
      n <- length(eval(as.name(vars[1L]), envir=environment(formula)))
      data <- n
      rm(vars)
    } else {
      n <- n_row(data)
    }
    y <- NULL
  }
  if (n < 1L) stop("empty data")

  self <- environment()
  family$sm <- self
  family$data <- data
  family$y <- y
  family$`_raw_` <- NULL
  family <- do.call(
    getFromNamespace(paste0("ff_", family[["family"]]), "mcmcsae"),
    family[-1L], envir=environment(formula)
  )
  y <- family[["y"]]

  Vmod <- NULL
  if (any(family[["family"]] == c("gaussian", "student_t", "gaussian_gamma"))) {
    if (family[["modeled.Q"]]) {
      Vmod <- family[["Vmod"]]
      if (any(names(mod) %in% names(Vmod))) stop("names of model components in 'var.model' must be distinct from names used in mean model 'formula'")
    }
  } else if (family[["family"]] == "multinomial") {
    n <- n * (family[["K"]] - 1L)
  }

  if (!prior.only) {
    if (control[["cMVN.or.CG"]]) {
      if (family[["family"]] == "gamma")
        stop("conjugate gradients and 'cMVN.sampler' algorithms not supported for 'gamma' family")
      if (!length(control[["block"]])) {
        warn("conjugate gradients and 'cMVN.sampler' algorithms currently only used for blocked Gibbs sampler")
        control$cMVN.sampler <- NULL
        control$CG <- NULL
        control$cMVN.or.CG <- FALSE
      }
    }

  }  # END if (!prior.only)

  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    mc$name <- names(mod)[k]
    mc$sc <- control
    mc$fam <- family
    mc$in.block <- any(mc[["name"]] == unlst(control[["block"]]))
    mc$prior.only <- prior.only
    mc$compute.weights <- compute.weights
    mc$data <- data
    mc <- as.list(mc)[-1L]
    mod[[mc[["name"]]]] <- do.call(
      getFromNamespace(paste0("mc_", types[k]), "mcmcsae"),
      mc, envir=environment(formula)
    )
  }
  labels_map <- list()  # information where to retrieve vector parameter labels
  for (mc in mod) {
    m <- fmatch(names(mc[["label_funs"]]), names(labels_map))
    if (!allNA(m)) stop("duplicate use of parameter name '", names(mc[["label_funs"]])[whichNA(m, invert=TRUE)], "'")
    labels_map <- c(labels_map, mc[["label_funs"]])
  }
  if (!is.null(Vmod)) {
    for (mc in Vmod) {
      m <- fmatch(names(mc[["label_funs"]]), names(labels_map))
      if (!allNA(m)) stop("duplicate use of parameter name '", names(mc[["label_funs"]])[whichNA(m, invert=TRUE)], "'")
      labels_map <- c(labels_map, mc[["label_funs"]])
    }
  }
  get_labels <- function(par.name) {
    labs <- labels_map[[par.name]]
    if (is.function(labs)) labs() else NULL
  }
  if (family[["family"]] == "poisson") {
    # add internal offset for negbinomial Poisson approximation to (first) mc_offset term
    mod[[whichv(types, "mc_offset")[1L]]]$add_internal_offset(-family[["log.shape"]])
  } else if (family[["family"]] == "multi") {
    for (fam in family[["fam.list"]]) {
      if (fam[["family"]] == "poisson") {
        mod[[whichv(types, "mc_offset")[1L]]]$add_internal_offset(-fam[["log.shape"]], sub=fam[["sub"]])
      }
    }
    rm(fam)
  }


  # compose following 3 functions:
  # draw: function to draw parameters from their full conditionals
  # rprior: function to draw parameters from their priors
  # start: function to create starting values (including for residuals/linear predictor e_)
  rprior <- function(p) {p <- list()}
  if (is.function(family$rprior)) rprior <- add(rprior, quote(p <- family$rprior(p)))

  for (k in seq_along(mod)) {
    if (is.function(mod[[k]][["rprior"]]))
      rprior <- add(rprior, bquote(p <- mod[[.(k)]]$rprior(p)))
  }
    
  if (is.null(linpred)) {
    do.linpred <- FALSE
  } else {
    if (identical(linpred, "fitted")) {
      linpred <- list()
      for (mc in mod) linpred[[mc$name]] <- mc$make_predict(verbose=FALSE)
    } else {
      if (!is.list(linpred) || length(linpred) == 0L) stop("'linpred' must be a non-empty list")
      if (!all(names(linpred) %in% names(mod)))
        stop("linpred names not corresponding to a model component: ", paste0(setdiff(names(linpred), names(mod)), collapse=", "))
      if (!all(i_apply(linpred[-1L], NROW) == NROW(linpred[[1L]]))) stop("not all matrices in 'linpred' have the same number of rows")
      for (k in names(linpred))
        linpred[[k]] <- mod[[k]]$make_predict(Xnew=linpred[[k]], verbose=FALSE)
    }
    do.linpred <- TRUE
    rprior <- add(rprior, quote(p$linpred_ <- lin_predict(p, linpred)))
  }

  lin_predict <- function(p, pred) {
    out <- pred[[1L]]$linpred(p)
    for (obj in pred[-1L]) obj$linpred_update(out, plus=TRUE, p)
    out
  }

  rprior <- add(rprior, quote(p))

  # What parameters to store by default? This information is used by MCMCsim.
  store_default <- function(prior=FALSE) {
    out <- if (prior || !control[["compute.llh"]]) NULL else "llh_"
    for (mc in mod) out <- c(out, mc[["store.default"]])
    out <- c(out, family$store_default())
    if (do.linpred) out <- c(out, "linpred_")
    out
  }
  store_mean_default <- function(prior.sampler=FALSE) {
    if (prior.sampler) {
      out <- NULL
    } else {
      out <- "e_"
      if (family[["modeled.Q"]] && any(family[["family"]] == c("gaussian", "student_t", "gaussian_gamma"))) out <- c(out, "Q_")
      if (compute.weights) out <- c(out, "weights_")
    }
    out
  }

  if (prior.only) {
    rm(k, mc, types, data)
    return(self)
  }

  if (compute.weights) {
    if (!do.linpred) stop("weights can only be computed for a linear predictor specified by argument 'linpred'")
    if (!control[["single.block"]]) stop("'compute.weights=TRUE' only supported for a fully blocked Gibbs sampler")
    if (family[["family"]] == "gamma") stop("weights computation not supported for gamma sampling distribution")
  }

  if (length(control[["block"]])) {
    mbs <- list()
    for (k in seq_along(control[["block"]]))
      mbs[[k]] <- create_mc_block(mod[control[["block"]][[k]]],
        family, control, prior.only=prior.only,
        compute.weights=compute.weights, linpred=linpred
      )
  }

  # compute residuals/linear predictor e_ for state p
  if (family[["e.is.res"]]) {  # residuals
    compute_e <- function(p) {
      e <- y - mod[[1L]]$lp(p)
      for (mc in mod[-1L]) mc$lp_update(e, plus=FALSE, p)
      e
    }
  } else {  # linear predictor
    compute_e <- function(p) {
      e <- mod[[1L]]$lp(p)
      for (mc in mod[-1L]) mc$lp_update(e, plus=TRUE, p)
      e
    }
  }
  # compute e for a collection of states, for use in waic/loo computation
  compute_e_i <- function(draws, i=seq_len(n)) {
    if (is.null(draws[["e_"]])) {
      nr <- n_chains(draws) * n_draws(draws)
      all.units <- length(i) == n
      if (family[["e.is.res"]]) {
        # residuals
        e_i <- matrix(rep_each(y[i], nr), nr, length(i))
        for (mc in mod)
          e_i <- e_i - mc$draws_linpred(draws, units = if (all.units) NULL else i, matrix=TRUE)
      } else {  # fitted values
        mc <- mod[[1L]]
        e_i <- mc$draws_linpred(draws, units = if (all.units) NULL else i, matrix=TRUE)
        for (mc in mod[-1L])
          e_i <- e_i + mc$draws_linpred(draws, units = if (all.units) NULL else i, matrix=TRUE)
      }
    } else {
      e_i <- as.matrix.dc(get_from(draws[["e_"]], vars=i))
    }
    e_i
  }

  if (control[["recompute.e"]]) {
    # compute residuals/linear predictor at each call of sampler (reduces build-up of rounding error?)
    draw <- function(p) {p$e_ <- compute_e(p)}
  } else {
    draw <- function(p) {}
  }
  start <- function(p=list()) {
    if (!is.list(p)) stop("input to 'start' function must be a list")
    if (offset.only) {
      p$e_ <- compute_e(p)
    } else {
      if (is.null(p[["e_"]])) {
        p$e_ <- Crnorm(n, sd=family[["scale.e"]])
      } else {
        if (length(p[["e_"]]) != n) stop("wrong length for 'e_' start value")
        p$e_ <- copy_obj(p[["e_"]])
      }
    }
  }
  MHpars <- NULL
  adapt <- function(ar) {}

  draw <- add(draw, quote(p <- family$draw(p)))
  if (is.function(family$start)) start <- add(start, quote(p <- family$start(p)))
  MHpars <- c(MHpars, family[["MHpars"]])
  if (is.function(family[["adapt"]])) adapt <- add(adapt, quote(family$adapt(ar)))

  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    if (mc[["type"]] == "gen") {
      if (mc[["gl"]] && mc[["usePX"]]) MHpars <- c(MHpars, mc[["name_xi"]])
      if (mc$strucA[["update.Q"]]) {
        MHpars <- c(MHpars, mc$strucA[["name_ext"]])
        if (mc$strucA$control[["adaptive"]])
          adapt <- add(adapt, bquote(mod[[.(k)]]$strucA$control$adapt(ar[[.(mc$strucA[["name_ext"]])]])))
      }
      if (!is.null(mc[["AR1.inferred"]])) {
        MHpars <- c(MHpars, mc[["name_AR1"]])
        if (mc$AR1sampler$MH[["adaptive"]])
          adapt <- add(adapt, bquote(mod[[.(k)]]$AR1sampler$MH$adapt(ar[[.(mc[["name_AR1"]])]])))
      }
      if (family[["family"]] == "gamma" && mc$control[["MHprop"]] == "LNRW") {
        MHpars <- c(MHpars, if (mc[["usePX"]]) mc[["name_sigma_raw"]] else mc[["name_sigma"]])
        adapt <- add(adapt, bquote(mod[[.(k)]]$adapt(ar)))
      }
    }
    if (is.function(mc[["start"]]))  start <- add(start, bquote(p <- mod[[.(k)]]$start(p)))
    if (is.function(mc[["draw"]])) draw <- add(draw, bquote(p <- mod[[.(k)]]$draw(p)))
  }

  if (length(control[["block"]])) {
    for (k in seq_along(mbs)) {
      draw <- add(draw, bquote(p <- mbs[[.(k)]]$draw(p)))
      start <- add(start, bquote(p <- mbs[[.(k)]]$start(p)))
    }
  }

  if (do.linpred)
    draw <- add(draw, quote(p$linpred_ <- lin_predict(p, linpred)))

  draw <- add(draw, quote(p))  # return state p

  if (!control[["recompute.e"]] && !control[["single.block"]]) {
    # adding this sometimes gives bad starting values in case of single.block (no need anyway in that case)
    start <- add(start, quote(p$e_ <- compute_e(p)))
  }
  start <- add(start, quote(p))
  if (length(body(adapt)) <= 1L) rm(adapt)

  fam_llh_i <- family[["llh_i"]]
  llh_i <- function(draws, i=seq_len(n)) {
    e_i <- compute_e_i(draws, i)
    fam_llh_i(draws, i, e_i)
  }

  if (!family[["family"]] == "multi" && !has.bart) {
    # max likelihood optimization
    # currently not available for multi-response models or models including a bart term
    # build a list of indices for vector representation of all likelihood parameters
    vec_list <- list()
    n_vec_list <- 0L
    # single block sampler assumes coefficients are at the start of the parameter vector
    for (k in seq_along(mod)) if (mod[[k]]$type != "mc_offset") {
      vec_list[[names(mod)[k]]] <- (n_vec_list + 1L):(n_vec_list + mod[[k]][["q"]])
      n_vec_list <- n_vec_list + mod[[k]][["q"]]
    }
    if (!family[["sigma.fixed"]]) {
      vec_list[["sigma_"]] <- n_vec_list + 1L
      n_vec_list <- n_vec_list + 1L
    }
    if (!is.null(Vmod)) {
      for (k in seq_along(Vmod)) {
        vec_list[[names(Vmod)[k]]] <- (n_vec_list + 1L):(n_vec_list + Vmod[[k]][["q"]])
        n_vec_list <- n_vec_list + Vmod[[k]][["q"]]
      }
    }
    switch(family[["family"]],
      negbinomial = if (!family[["shape.fixed"]]) {
        vec_list[["negbin_shape_"]] <- n_vec_list + 1L
        n_vec_list <- n_vec_list + 1L
      },
      gamma=, gaussian_gamma = if (!family[["alpha.fixed"]]) {
        vec_list[["gamma_shape_"]] <- n_vec_list + 1L
        n_vec_list <- n_vec_list + 1L
      }
    )
    pars <- names(vec_list)
    vec2list <- function(x) {
      out <- list()
      for (k in seq_along(vec_list)) out[[pars[k]]] <- x[vec_list[[k]]]
      out
    }
    list2vec <- function(p) {
      out <- rep.int(NA_real_, n_vec_list)
      for (k in seq_along(vec_list)) out[vec_list[[k]]] <- p[[pars[k]]]
      out
    }
    # log-likelihood function for optimisation: function of a vector instead of list
    llh_opt <- function(x) {
      p <- vec2list(x)
      p$e_ <- compute_e(p)
      if (any(family[["family"]] == c("gaussian", "student_t", "gaussian_gamma")))
        p$SSR_ <- dotprodC(p[["e_"]], family$Q_e(p))
      family$llh(p)
    }
  }

  # remove quantities no longer needed
  rm(k, mc, types, data)

  # return the function environment, including draw, rprior, start functions
  self
}


#' Set computational options for the sampling algorithms
#'
#' This function can be used to specify computational options
#' by passing the result to the \code{control} argument of
#' \code{\link{create_sampler}}.
#'
#' @export
#' @param add.outer.R whether to add the outer product of a constraint matrix to the
#'  conditional posterior precision matrix of coefficients sampled in a block. This is used
#'  to resolve singularity due to intrinsic GMRF components. Default is \code{TRUE}.
#'  When set to \code{NULL}, a simple heuristic is used to decide whether
#'  to add the outer product of (possibly a submatrix of) the constraint matrix.
#' @param add.eps.I whether to add a small positive multiple of the identity matrix
#'  to the conditional posterior precision matrix of coefficients sampled in a block.
#'  If needed, this can resolve singularity as an alternative to \code{add.outer.R=TRUE}.
#'  The advantage of \code{add.eps.I=TRUE} is that a sparse conditional posterior precision
#'  matrix remains sparse so that sampling is faster, at the cost of slightly deviating
#'  from the target posterior distribution, depending on the value of \code{eps}.
#'  If \code{add.eps.I=TRUE} \code{add.outer.R} will be set to FALSE.
#' @param eps a positive scalar value, used only in case \code{add.eps.I=TRUE}. This
#'  should be a small value to ensure that one is not deviating too much from the
#'  desired posterior distribution of coefficients sampled in a block. On the other
#'  hand, if it is chosen too small it may not resolve the singularity of the conditional
#'  posterior precision matrix of coefficients sampled in a block.
#' @param recompute.e when \code{FALSE}, residuals or linear predictors are only computed at the start of the simulation.
#'  This may give a modest speed-up but in some cases may be less accurate due to round-off error accumulation.
#'  Default is \code{TRUE}.
#' @param compute.llh whether to compute the log-likelihood for each MCMC draw. Default is
#'  \code{TRUE}, but setting it to \code{FALSE} may reduce computation time a little.
#' @param cMVN.sampler whether an extended linear system including dual variables is used
#'  for equality constrained multivariate normal sampling. If set to \code{TRUE} this may
#'  improve the performance of the blocked Gibbs sampler, especially in case of a large number
#'  of equality constraints, typically (intrinsic) GMRF identifiability constraints.
#' @param CG use a conjugate gradient iterative algorithm instead of Cholesky updates for sampling
#'  the model's coefficients. This must be a list with possible components \code{max.it},
#'  \code{stop.criterion}, \code{verbose}, \code{preconditioner} and \code{scale}.
#'  See the help for function \code{\link{CG_control}}, which can be used to specify these options.
#'  Conjugate gradient sampling is currently an experimental feature that can be used for
#'  blocked Gibbs sampling but with some limitations.
#' @param block if \code{TRUE}, the default, all coefficients are sampled in a single Gibbs block.
#'  If \code{FALSE}, the coefficients of each model component are sampled separately in sequence.
#'  Alternatively, a list of character vectors with names of model components can be passed to
#'  specify a grouping of model components whose coefficients should be sampled together in blocks.
#' @param auto.order.block whether Gibbs blocks should be ordered automatically in such a
#'  way that those with the most sparse design matrices come first. This way of ordering
#'  can make Cholesky updates more efficient.
#' @param chol.control options for Cholesky decomposition, see \code{\link{chol_control}}.
#' @param max.size.cps.template maximum allowed size in MB of the sparse matrix serving as a 
#'  template for the sparse symmetric crossproduct X'QX of a dgCMatrix X, where Q is a diagonal
#'  matrix subject to change.
#' @returns A list with specified computational options used by various sampling functions.
#' @references
#'  D. Bates, M. Maechler, B. Bolker and S.C. Walker (2015).
#'    Fitting Linear Mixed-Effects Models Using lme4.
#'    Journal of Statistical Software 67(1), 1-48.
#'
#'  Y. Chen, T.A. Davis, W.W. Hager and S. Rajamanickam (2008).
#'    Algorithm 887: CHOLMOD, supernodal sparse Cholesky factorization and update/downdate.
#'    ACM Transactions on Mathematical Software 35(3), 1-14.
sampler_control <- function(add.outer.R=TRUE, add.eps.I=FALSE, eps=sqrt(.Machine$double.eps),
                            recompute.e=TRUE, compute.llh=TRUE,
                            cMVN.sampler=NULL, CG=NULL,
                            block=TRUE, auto.order.block=TRUE,
                            chol.control=chol_control(),
                            max.size.cps.template=100) {
  list(add.outer.R=add.outer.R, add.eps.I=add.eps.I, eps=eps,
       recompute.e=recompute.e, compute.llh=compute.llh,
       cMVN.sampler=cMVN.sampler, CG=CG,
       block=block, auto.order.block=auto.order.block,
       chol.control=chol.control,
       max.size.cps.template=max.size.cps.template
  )
}

check_sampler_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- sampler_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  if (!is_logical_scalar(control[["recompute.e"]])) stop("'recompute.e' must be TRUE or FALSE")
  if (!is_logical_scalar(control[["compute.llh"]])) stop("'compute.llh' must be TRUE or FALSE")
  if (isTRUE(control[["add.eps.I"]])) {
    if (!(is_numeric_scalar(control[["eps"]]) && control[["eps"]] > 0)) stop("'eps' must be a single positive numerical value")
    control$add.outer.R <- FALSE
  }
  if (is.logical(control[["block"]])) {
    if (length(control[["block"]]) != 1L) stop("unexpected input for 'block'")
  } else {
    if (!is.list(control[["block"]])) stop("'block' should be either a scalar logical or a list of model component name vectors")
    if (!length(control[["block"]])) stop("'block' must contain at least one character vector")
  }
  control$chol.control <- check_chol_control(control[["chol.control"]])
  if (isTRUE(control[["CG"]])) {
    control$CG <- CG_control()
  } else if (isFALSE(control[["CG"]])) {
    control$CG <- NULL
  } else if (!is.null(control[["CG"]])) {
    control$CG <- check_CG_control(control[["CG"]])
  }
  if (isTRUE(control[["cMVN.sampler"]])) {
    control$cMVN.sampler <- cMVN_control()
  } else if (isFALSE(control[["cMVN.sampler"]])) {
    control$cMVN.sampler <- NULL
  } else if (!is.null(control[["cMVN.sampler"]])) {
    control$cMVN.sampler <- check_cMVN_control(control[["cMVN.sampler"]])
  }
  if (is.list(control[["cMVN.sampler"]]) && is.list(control[["CG"]])) stop("'cMVN.sampler' and 'CG' cannot be combined")
  control$cMVN.or.CG <- is.list(control[["cMVN.sampler"]]) || is.list(control[["CG"]])
  control
}
