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
#' separately, in the usual way using \code{\link{offset}(...)}.
#'
#' For gaussian models, \code{formula.V} can be used to specify the variance structure of the model.
#' Currently two specialised variance model components are supported, \code{\link{vreg}(...)} for regression
#' effects predicting the log-variance and \code{\link{vfac}(...)} for modelled variance factors.
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
#'  for different transformations. It should be supplied as a vector of the same size as the response variable,
#'  and is currently only supported if \code{family="gaussian"}.
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
    } else
      data <- as.data.frame(data)
  }

  if (is.function(family)) {
    family <- as.character(substitute(family))
    if (startsWith(family, "f_")) family <- substring(family, 3L)
  }
  if (is.character(family)) {
    family <- tolower(family)  # so that e.g. both "gamma" and "Gamma" allowed
    family <- match.arg(family, c("gaussian", "binomial", "negbinomial", "poisson", "multinomial", "gamma", "gaussian_gamma"))
    family <- eval(call(paste0("f_", family)))
  }

  if (!(is.null(sigma.fixed) && is.null(sigma.mod) && is.null(Q0) && is.null(formula.V) && is.null(logJacobian))) {
    if (family[["family"]] == "gaussian") {
      warn("arguments 'sigma.fixed', 'sigma.mod', 'Q0', 'formula.V' and 'logJacobian' are deprecated\n",
           "please pass variance options to argument 'family' using function f_gaussian\n")
      # TODO if non-default options have been set using f_gaussian raise an error
      family <- f_gaussian(
        var.prior = if (isTRUE(sigma.fixed)) pr_fixed(1) else if (!is.null(sigma.mod)) sigma.mod else pr_invchisq(0, 1),
        var.vec =
          if (is.null(Q0) || is_a_matrix(Q0)) {
            ~ 1
          } else if (is.character(Q0)) {
            as.formula(paste0("~ I(1/", Q0, ")"))
          } else if (is.vector(Q0)) {
            1/Q0
          } else {
            as.formula(paste0("~ I(1/(", deparse(substitute(Q0)), "))"))
          },
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
  formula <- standardize_formula(formula, data=data, internal.offset = family[["family"]] == "poisson")
  mod <- to_mclist(formula)
  if (!length(mod)) stop("empty 'formula'")
  types <- get_types(mod)
  if (any(types %in% c("vreg", "vfac"))) stop("'vreg' and 'vfac' can only be used in variance model specification")
  has.bart <- any(types == "brt")

  control <- check_sampler_control(control)

  if (has_response(formula)) {
    y <- get_response(formula, data)
    if (is.numeric(y)) {
      if (!all(is.finite(y))) stop(sum(!is.finite(y)), " missing or infinite value(s) in response variable")
    } else {
      if (anyNA(y)) stop(sum(is.na(y)), " missing value(s) in response variable")
    }
    n <- NROW(y)
    if (is.null(data)) data <- n
    y <- family$init(data, y)
  } else {
    if (!prior.only) {
      warn("no left hand side in 'formula': setting up prior samplers only")
      prior.only <- TRUE
    }
    if (is.null(data)) {
      if (!length(all.vars(formula))) stop("no way to determine the number of cases")
      n <- length(get(all.vars(formula)[1L]))
      data <- n
    } else {
      n <- n_row(data)
    }
    family$init(data)
  }
  if (n < 1L) stop("empty data")

  Vmod <- NULL
  if (any(family[["family"]] == c("gaussian", "gaussian_gamma"))) {
    if (prior.only || n == 1L) {
      scale.e <- 1
    } else {
      # scale.e used to generate default starting values for residuals (or fitted values for non-gaussian models)
      # divide by number of model components + 1: each component explains part of the total variation
      scale.e <- sd(y) / (length(mod) + 1L)
      if (scale.e == 0 || !is.finite(scale.e)) scale.e <- 1
    }
    sigma.fixed <- family[["sigma.fixed"]]
    modeled.Q <- family[["modeled.Q"]]
    Q0.type <- family[["Q0.type"]]
    Q0 <- family[["Q0"]]
    scale.sigma <- scale.e * fmean.default(sqrt(diag(Q0)), na.rm=FALSE)
    if (!sigma.fixed) sigma.mod <- family[["var.prior"]]
    if (modeled.Q) {
      Vmod <- family[["Vmod"]]
      if (any(names(mod) %in% names(Vmod))) stop("names of model components in 'var.model' must be distinct from names used in mean model 'formula'")
    }
    f_mean <- identity  # mean function acting on linear predictor
  } else {  # binomial, multinomial, negative binomial, Poisson or gamma likelihood
    sigma.fixed <- TRUE
    Q0 <- CdiagU(n)
    Q0.type <- "unit"
    modeled.Q <- family[["link"]] != "probit" && family[["family"]] != "gamma"
    scale.e <- if (family[["family"]] == "gamma") 1 else 2.5  # is 2.5 reasonable for all other cases?
    scale.sigma <- scale.e
    if (family[["family"]] == "binomial") {
      ny <- family[["ny"]]
      f_mean <- function(eta) ny / (1 + exp(-eta))  # mean function acting on linear predictor
    } else if (family[["family"]] == "multinomial") {
      cats <- family[["cats"]]
      ny <- family[["ny"]]
      Km1 <- family[["K"]] - 1L
      n0 <- n
      n <- n0 * Km1
    } else {  # negative binomial, Poisson or gamma likelihood
      if (family[["family"]] == "negbinomial") {
        if (family[["shape.fixed"]]) {
          if (!prior.only) {
            ry <- family[["shape0"]]
            ny <- y + ry  # used in Polya-Gamma full conditional
          }
          f_mean <- function(eta) family[["shape0"]] * exp(eta)  # mean function acting on linear predictor
        } else {
          if (!prior.only) y <- as.numeric(y)  # for C++ rCRT function; alternatively force to be integer
          f_mean <- function(eta) stop("negative binomial family with modelled shape parameter: mean also depends on this parameter")
        }
      } else if (family[["family"]] == "poisson") {
        ry <- family$control[["nb.shape"]]
        if (!prior.only) ny <- y + ry  # used in Polya-Gamma full conditional
        f_mean <- function(eta) exp(eta)  # mean function acting on linear predictor
      }
    }
    if (!prior.only && family[["link"]] != "probit") {
      if (control[["PG.approx"]]) {
        mPG <- as.integer(control[["PG.approx.m"]])
        if (all(length(mPG) != c(1L, n))) stop("invalid value for option 'PG.approx.m'")
        rPolyaGamma <- function(b, c) CrPGapprox(n, b, c, mPG)
      } else {
        if (!requireNamespace("BayesLogit", quietly=TRUE)) stop("please install package 'BayesLogit' and try again")
        if (family[["family"]] == "binomial") {
          ny[ny == 0] <- .Machine$double.eps  # to prevent generation of NA by BayesLogit::rpg
        }
        rpg <- BayesLogit::rpg
        rPolyaGamma <- function(b, c) rpg(n, b, c)
      }
    }
  }

  e.is.res <- family[["link"]] == "identity"
  if (e.is.res)
    y_eff <- function() copy_vector(y)
  else
    y_eff <- function() 0

  if (family[["family"]] == "gamma" && !all(types %in% c("mc_offset", "reg", "gen"))) stop("only 'reg' and 'gen' (and offset/mc_offset) components supported for family 'gamma'")
  # check Gibbs blocks for mean model specification
  if (!is.null(block)) {
    warn("'block' argument of create_sampler is deprecated; ",
      "instead, please use argument 'control' to specify the Gibbs sampling blocks; ",
      "by default a maximal blocking strategy is used so that coefficients are sampled together whenever possible")
    control$block <- block
  }
  block <- control[["block"]]
  if (is.logical(block)) {
    if (block) {
      # all components of type reg, gen and mec in a single block
      block <- list(names(mod)[!(types %in% c("mc_offset", "brt"))])
      # single component by default not handled as a block, unless
      #   compute.weights=TRUE, or expanded.cMVN or CG sampler is used
      if (!length(block[[1L]]) || (length(block[[1L]]) == 1L && !compute.weights && !is.list(control[["CG"]]) && !control[["expanded.cMVN.sampler"]]))
        block <- list()
    } else
      block <- list()
  } else {
    for (k in seq_along(block)) {
      if (!all(block[[k]] %in% names(mod)))
        stop("invalid name(s) '", paste0(setdiff(block[[k]], names(mod)), collapse="', '"), "' in 'block'")
      if (any(types[block[[k]]] == "brt")) stop("'brt' model component cannot be part of a Gibbs block")
      if (any(types[block[[k]]] == "mc_offset")) stop("'offset/mc_offset' model component cannot be part of a Gibbs block")
    }
    if (any_duplicated(unlst(block))) stop("duplicate model component name(s) in 'block'")
  }
  single.block <- any(length(mod) == c(1L, length(unlst(block))))
  offset.only <- all(types == "mc_offset")

  if (!prior.only) {
    if (control[["expanded.cMVN.sampler"]] || !is.null(control[["CG"]])) {
      if (family[["family"]] == "gamma")
        stop("conjugate gradients and 'expanded.cMVN.sampler' algorithms not supported for 'gamma' family")
      if (!length(block)) {
        warn("conjugate gradients and 'expanded.cMVN.sampler' algorithms currently only used for blocked Gibbs sampler")
        control$expanded.cMVN.sampler <- FALSE
        control$CG <- NULL
      }
    }
    if (single.block) {  # computation of working response may simplify
      control$recompute.e <- FALSE  # already computed in the coefficient sampling function
      y <- as.numeric(y)  # some C++ matrix algebra functions expect double
    }

    # function Q_e computes matrix-vector product of Q and working response e_
    if (any(family[["family"]] == c("gaussian", "gaussian_gamma"))) {
      if (modeled.Q) {
        Q_e <- switch(Q0.type,
          unit=, diag = function(p) p[["Q_"]] * p[["e_"]],
          # in case of non-diagonal Q0, full precision matrix stored as p[["QM_"]]
          symm = function(p) p[["QM_"]] %m*v% p[["e_"]]
        )
      } else {
        Q_e <- switch(Q0.type,
          unit = function(p) p[["e_"]],
          diag = function(p) Q0@x * p[["e_"]],
          symm = function(p) Q0 %m*v% p[["e_"]]
        )
      }
    } else if (family[["family"]] != "gamma") {
      # in this case Q0 is ignored
      if (family[["link"]] == "probit") {
        if (single.block)
          Q_e <- function(p) p[["z_"]]
        else
          Q_e <- function(p) p[["z_"]] - p[["e_"]]
      } else if (family[["family"]] == "negbinomial" && !family[["shape.fixed"]]) {
        if (single.block)
          Q_e <- function(p) 0.5 * (y - family$get_shape(p))
        else
          Q_e <- function(p) 0.5 * (y - family$get_shape(p)) - p[["Q_"]] * p[["e_"]]
      } else {
        if (family[["family"]] == "negbinomial" || family[["family"]] == "poisson")
          y_shifted <- 0.5 * (y - ry)
        else  # binomial or multinomial
          y_shifted <- y - 0.5 * family[["ny"]]
        if (single.block)
          Q_e <- function(p) y_shifted
        else
          Q_e <- function(p) y_shifted - p[["Q_"]] * p[["e_"]]
      }
    }
    if (!is.null(control[["CG"]]) || control[["expanded.cMVN.sampler"]]) {
      # set up a function that multiplies by L Chol factor of Q, for sampling from N(., Q)
      if (any(family[["family"]] == c("gaussian", "gaussian_gamma"))) {
        if (modeled.Q) {
          cholQ <- switch(Q0.type,
            unit=, diag = build_chol(runif(n, 0.9, 1.1)),
            symm = build_chol(crossprod_sym(Cdiag(runif(0.9, 1.1)), Q0), control=control$chol.control)
          )
        } else {
          cholQ <- switch(Q0.type,
            unit = build_chol(CdiagU(n)),
            diag = build_chol(Cdiag(Q0@x)),
            symm = build_chol(Q0, control=control$chol.control)
          )
        }
      } else {
        if (family[["link"]] == "probit") {
          cholQ <- build_chol(CdiagU(n))
        } else {
          cholQ <- build_chol(runif(n, 0.9, 1.1))
        }
      }
      # draw from MVN with variance(!) Q
      drawMVNvarQ <- function(p) {
        if (modeled.Q) {
          if (Q0.type == "symm")
            cholQ$update(p[["QM_"]])
          else
            cholQ$update(p[["Q_"]])
        }
        cholQ$Ltimes(Crnorm(n), transpose=FALSE)
      }
    }
  }  # END if (!prior.only)

  coef.names <- vector(mode="list", length(mod))
  names(coef.names) <- names(mod)
  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    mc$name <- names(mod)[k]
    mc <- as.list(mc)[-1L]
    mod[[mc[["name"]]]] <- do.call(types[k], mc, envir=environment(formula))
  }
  if (family[["family"]] == "poisson") {
    # add internal offset for negbinomial Poisson approximation to (first) mc_offset term
    mod[[whichv(types, "mc_offset")[1L]]]$add_internal_offset(-log(ry))
  }


  # compose following 3 functions:
  # draw: function to draw parameters from their full conditionals
  # rprior: function to draw parameters from their priors
  # start: function to create starting values (including for residuals/linear predictor e_)
  rprior <- function(p) {p <- list()}
  switch(family[["family"]],
    gaussian = if (!family[["sigma.fixed"]] || family[["modeled.Q"]]) rprior <- add(rprior, quote(p <- family$rprior(p))),
    negbinomial = if (!family[["shape.fixed"]]) rprior <- add(rprior, quote(p[["negbin_shape_"]] <- 1 / family$inv.shape.prior$rprior())),
    gamma=, gaussian_gamma = if (!family[["alpha.fixed"]]) rprior <- add(rprior, quote(p[["gamma_shape_"]] <- family$shape.prior$rprior()))
  )

  if (any(family[["family"]] == c("gaussian", "gaussian_gamma"))) {
    if (modeled.Q) {
      # check Gibbs blocks for variance model specification
      types <- family[["types"]]
      if (is.logical(control[["block.V"]])) {
        if (control[["block.V"]]) {
          # all components of type reg and gen in a single block
          block.V <- list(names(Vmod)[types %in% c("reg", "gen")])
          # a single component by default not handled as a block
          if (length(block.V[[1L]]) <= 1L) block.V <- NULL
        } else {
          block.V <- NULL
        }
      } else {
        block.V <- control[["block.V"]]
        for (k in seq_along(block.V)) {
          if (!all(block.V[[k]] %in% names(Vmod)))
            stop("invalid name(s) '", paste0(setdiff(block.V[[k]], names(Vmod)), collapse="', '"), "' in 'block.V'")
          if (any(types[block.V[[k]]] %in% c("vreg", "vfac"))) stop("'vreg' and 'vfac' components cannot be part of a Gibbs block")
        }
        if (any_duplicated(unlst(block.V))) stop("duplicate model component name in 'block.V'")
      }
      single.V.block <- any(length(Vmod) == c(1L, length(unlst(block.V))))

      for (k in seq_along(Vmod)) {
        mc <- Vmod[[k]]
        mc$name <- names(Vmod)[k]
        mc <- as.list(mc)[-1L]
        Vmod[[mc[["name"]]]] <- do.call(types[k], mc, envir=parent.frame())
      }
      family$set_Vmod(Vmod)
    } else {
      block.V <- NULL
    }
  }

  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    switch(mc[["type"]],
      reg = rprior <- add(rprior, bquote(p[[.(mc[["name"]])]] <- mod[[.(k)]]$rprior(p))),
      mec = rprior <- add(rprior, bquote(p <- mod[[.(k)]]$rprior(p))),
      gen = rprior <- add(rprior, bquote(p <- mod[[.(k)]]$rprior(p)))
    )
  }

  if (!has.bart) {
    # for max likelihood optimization, currently not available for models including a bart term
    # build a list of indices for vector representation of all likelihood parameters
    vec_list <- list()
    n_vec_list <- 0L
    # single block sampler assumes coefficients are at the start of the parameter vector
    for (k in seq_along(mod)) if (mod[[k]]$type != "mc_offset") {
      vec_list[[names(mod)[k]]] <- (n_vec_list + 1L):(n_vec_list + mod[[k]][["q"]])
      n_vec_list <- n_vec_list + mod[[k]][["q"]]
    }
    if (!sigma.fixed) {
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
    vec2list <- function(x) {
      pars <- names(vec_list)
      out <- list()
      for (k in seq_along(vec_list)) out[[pars[k]]] <- x[vec_list[[k]]]
      out
    }
    list2vec <- function(p) {
      pars <- names(vec_list)
      out <- rep.int(NA_real_, n_vec_list)
      for (k in seq_along(vec_list)) out[vec_list[[k]]] <- p[[pars[k]]]
      out
    }
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
  store_default <- function(prior.sampler=FALSE) {
    out <- if (prior.sampler) NULL else "llh_"
    if (!sigma.fixed) out <- c(out, "sigma_")
    if (do.linpred) out <- c(out, "linpred_")
    switch(family[["family"]],
      negbinomial = if (!family[["shape.fixed"]]) out <- c(out, "negbin_shape_"),
      gamma=, gaussian_gamma = if (!family[["alpha.fixed"]]) out <- c(out, "gamma_shape_")
    )
    for (mc in mod)
      switch(mc[["type"]],
        reg = out <- c(out, mc[["name"]]),
        mec = out <- c(out, mc[["name"]]),
        gen = {
          out <- c(out, mc[["name_sigma"]])
          if (mc[["var"]] == "unstructured") out <- c(out, mc[["name_rho"]])
          if (mc[["gl"]]) out <- c(out, mc[["name_gl"]])
          if (mc$strucA[["update.Q"]]) out <- c(out, mc$strucA[["name_ext"]])
          if (!is.null(mc[["priorA"]]) && is.list(mc$priorA[["df"]])) out <- c(out, mc[["name_df"]])
          if (!is.null(mc[["AR1.inferred"]])) out <- c(out, mc[["name_AR1"]])
        },
        brt = if (mc[["keepTrees"]]) out <- c(out, mc[["name_trees"]])
      )
    for (mc in Vmod)
      switch(mc[["type"]],
        vreg=, reg = out <- c(out, mc[["name"]]),
        vfac = if (mc$prior[["type"]] == "invchisq" && is.list(mc$prior[["df"]])) out <- c(out, mc[["name_df"]]),
        gen = out <- c(out, mc[["name_sigma"]])
      )
    out
  }
  store_mean_default <- function(prior.sampler=FALSE) {
    if (prior.sampler) {
      out <- NULL
    } else {
      out <- "e_"
      if (modeled.Q && any(family[["family"]] == c("gaussian", "gaussian_gamma"))) out <- c(out, "Q_")
      if (compute.weights) out <- c(out, "weights_")
    }
    out
  }

  if (prior.only) {
    rm(k, mc, types, data)
    return(environment())
  }


  if (compute.weights) {
    if (!do.linpred) stop("weights can only be computed for a linear predictor specified by argument 'linpred'")
    if (!single.block) stop("'compute.weights=TRUE' cannot be combined with 'block=FALSE'")
    if (family[["family"]] == "gamma") stop("weights computation not supported for gamma sampling distribution")
  }

  if (length(block)) {
    mbs <- list()
    for (k in seq_along(block)) mbs[[k]] <- create_mc_block(mod[block[[k]]])
  }

  if (any(family[["family"]] == c("gaussian", "gaussian_gamma")) && length(block.V)) {
    mbs.V <- list()
    for (k in seq_along(block.V)) mbs.V[[k]] <- create_mc_block(Vmod[block.V[[k]]])
  }

  # compute residuals/linear predictor e_ for state p
  if (e.is.res) {  # residuals
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
      if (e.is.res) {
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
        p$e_ <- Crnorm(n, sd=scale.e)
      } else {
        if (length(p[["e_"]]) != n) stop("wrong length for 'e_' start value")
        p$e_ <- copy_vector(p[["e_"]])
      }
    }
  }
  MHpars <- NULL
  adapt <- function(ar) {}

  if ((family[["family"]] == "gaussian" || family[["family"]] == "gaussian_gamma") && modeled.Q) {
    if (!single.V.block) {
      # recompute precision Q for numerical stability
      draw <- add(draw, quote(p <- family$compute_Q(p)))
    }
    for (k in seq_along(Vmod)) {
      mc <- Vmod[[k]]
      switch(mc[["type"]],
        vreg = MHpars <- c(MHpars, mc[["name"]]),
        vfac = {
          if (mc$prior[["type"]] == "invchisq" && is.list(mc$prior[["df"]])) {
            MHpars <- c(MHpars, mc[["name_df"]])
            if (mc$prior$df[["adapt"]])
              adapt <- add(adapt, bquote(Vmod[[.(k)]]$adapt(ar)))
          }
        },
        gen = {
          MHpars <- c(MHpars, if (mc[["usePX"]]) mc[["name_sigma_raw"]] else mc[["name_sigma"]])
          if (mc$control[["MHprop"]] == "LNRW")
            adapt <- add(adapt, bquote(Vmod[[.(k)]]$adapt(ar)))
        }
      )
      if (!(mc[["type"]] == "reg" && mc[["in_block"]])) {
        draw <- add(draw, bquote(p <- Vmod[[.(k)]]$draw(p)))
        start <- add(start, bquote(p <- Vmod[[.(k)]]$start(p)))
      }
    }
    if (length(block.V)) {
      for (k in seq_along(mbs.V)) {
        draw <- add(draw, bquote(p <- mbs.V[[.(k)]]$draw(p)))
        start <- add(start, bquote(p <- mbs.V[[.(k)]]$start(p)))
      }
    }
    start <- add(start, quote(p <- family$compute_Q(p)))
  }
  llh <- family$make_llh(y)
  if (family[["family"]] == "gaussian" || family[["family"]] == "gaussian_gamma") {
    draw <- draw |>
      add(quote(SSR <- dotprodC(p[["e_"]], Q_e(p)))) |>
      add(quote(p$llh_ <- llh(p, SSR)))
  } else if (family[["family"]] == "gamma") {
    draw <- add(draw, quote(p$llh_ <- llh(p)))
  }
  if (family[["family"]] == "gamma" || family[["family"]] == "gaussian_gamma") {
    if (!family[["alpha.fixed"]]) {
      if (family[["family"]] == "gamma")
        draw_shape <- family$make_draw_shape(y)
      else  # gaussian_gamma
        draw_shape <- family$make_draw_shape(family[["sigmasq"]])
      draw <- add(draw, quote(p$gamma_shape_ <- draw_shape(p)))
      # shape start value should not be too small
      start <- add(start, quote(if (is.null(p[["gamma_shape_"]])) p$gamma_shape_ <- exp(runif(1L, -3, 4))))
      MHpars <- c(MHpars, "gamma_shape_")
      if (family$control[["type"]] == "RWLN" && family$control[["adaptive"]]) {
        adapt <- add(adapt, quote(family$control$adapt(ar[["gamma_shape_"]])))
      }
    }
  } else if (family[["family"]] == "negbinomial") {
    draw <- add(draw, quote(p$llh_ <- llh(p)))
    if (!family[["shape.fixed"]]) {
      if (family$inv.shape.prior[["type"]] == "fixed") {
        start <- add(start, quote(if (is.null(p[["negbin_shape_"]])) p$negbin_shape_ <- 1/family$inv.shape.prior$rprior()))
      } else {
        start <- add(start, quote(
          if (is.null(p[["negbin_shape_"]])) {
            p$negbin_shape_ <- runif(1L, 0.1, 10)
          } else {
            if (length(p[["negbin_shape_"]]) != 1L) stop("wrong length for 'negbin_shape_' start value")
            p[["negbin_shape_"]] <- as.numeric(p[["negbin_shape_"]])
            if (is.na(p[["negbin_shape_"]]) || p[["negbin_shape_"]] <= 0) stop("'negbin_shape_' start value must be a positive number")
          }
        ))
      }
      draw_shape <- family$make_draw_shape(y)
      draw <- add(draw, quote(p$negbin_shape_ <- draw_shape(p)))
    }
    draw <- add(draw, quote(ny <- y + family$get_shape(p)))  # used in Polya-Gamma full conditional for latent precision vector
    start <- add(start, quote(ny <- y + family$get_shape(p)))
  }
  if (family[["link"]] == "probit") {  # can only be binomial family
    draw <- draw |>
      add(quote(p[["z_"]] <- CrTNprobit(p[["e_"]], y))) |>
      add(quote(p$llh_ <- llh(p)))
    start <- start |>
      add(quote(if (is.null(p[["z_"]])) p$z_ <- CrTNprobit(p[["e_"]], y))) |>
      add(quote(if (length(p[["z_"]]) != n) stop("wrong length for 'z_' start value")))
  } else if (any(family[["family"]] == c("binomial", "multinomial", "negbinomial", "poisson"))) {
    draw <- draw |>
      add(quote(p[["Q_"]] <- rPolyaGamma(ny, p[["e_"]]))) |>
      add(quote(p$llh_ <- llh(p)))
    start <- start |>
      add(quote(if (is.null(p[["Q_"]])) p$Q_ <- rPolyaGamma(ny, p[["e_"]]))) |>
      add(quote(if (length(p[["Q_"]]) != n) stop("wrong length for 'Q_' start value")))
  }

  if (!sigma.fixed) {
    draw_sigma <- function(p, SSR) {}
    # compute df.data + update SSR with contributions from reg and gen components
    df.data <- n
    for (k in seq_along(mod)) {
      mc <- mod[[k]]
      switch(mc[["type"]],
        reg=, mec = {
          if (mc$prior[["type"]] == "normal" && mc[["informative.prior"]]) {
            if (is.null(mc[["R"]]))
              df.data <- df.data + mc[["q"]]
            else
              df.data <- df.data + mc[["q"]] - ncol(mc[["R"]])
            if (mc[["zero.mean"]])
              draw_sigma <- add(draw_sigma, bquote(delta.beta <- p[[.(mc[["name"]])]]))
            else
              draw_sigma <- add(draw_sigma, bquote(delta.beta <- p[[.(mc[["name"]])]] - mod[[.(k)]]$prior[["mean"]]))
            draw_sigma <- add(draw_sigma, bquote(SSR <- SSR + dotprodC(delta.beta, mod[[.(k)]][["Q0"]] %m*v% delta.beta)))
          }
        },
        gen = {
          if (mc[["usePX"]] && mc$PX[["data.scale"]]) {
            df.data <- df.data + mc$PX[["dim"]]
            draw_sigma <- add(draw_sigma, bquote(SSR <- SSR + dotprodC(p[[.(mc[["name_xi"]])]], mod[[.(k)]][["PX_Q0"]] %m*v% p[[.(mc[["name_xi"]])]])))
          }
          if (mc[["gl"]] && mc$glp[["informative.prior"]]) {
            df.data <- df.data + mc$glp[["q"]]
            #draw_sigma <- add(draw_sigma, bquote(delta.beta <- p[[.(mc$name_gl)]] - mod[[.(k)]]$glp$b0))
            draw_sigma <- draw_sigma |>
              add(bquote(delta.beta <- p[[.(mc[["name_gl"]])]])) |>
              add(bquote(SSR <- SSR + dotprodC(delta.beta, mod[[.(k)]]$glp[["Q0"]] %m*v% delta.beta)))
          }
        }
      )
    }
    switch(sigma.mod[["type"]],
      fixed = {
        draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod[["value"]])))
      },
      invchisq = {
        if (is.list(sigma.mod[["scale"]]))
          draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod$draw(df.data, SSR, 1 / p[["sigma_"]]^2))))
        else
          draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(sigma.mod$draw(df.data, SSR))))
      },
      exp = {
        rgig <- GIGrvg::rgig
        draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(rgig(1L, 1 - 0.5*df.data, SSR, 2/sigma.mod[["scale"]]))))
      },
      gig = {
        rgig <- GIGrvg::rgig
        draw_sigma <- add(draw_sigma, quote(p$sigma_ <- sqrt(rgig(1L, sigma.mod[["p"]] - 0.5*df.data, sigma.mod[["b"]] + SSR, sigma.mod[["a"]]))))
      }
    )
    draw_sigma <- add(draw_sigma, quote(p))

    draw <- add(draw, quote(p <- draw_sigma(p, SSR)))

    if (sigma.mod[["type"]] == "fixed") {
      start <- add(start, quote(if (is.null(p[["sigma_"]])) p$sigma_ <- sqrt(sigma.mod$rprior())))
    } else {
      start <- add(start, quote(
        if (is.null(p[["sigma_"]]))
          p$sigma_ <- runif(1L, 0.1 * scale.sigma, scale.sigma)
        else if (length(p[["sigma_"]]) != 1L)
          stop("wrong length for 'sigma_' start value")
      ))
    }
  }  # END if (!sigma.fixed)

  for (k in seq_along(mod)) {
    mc <- mod[[k]]
    if (mc[["type"]] == "gen") {
      if (mc[["gl"]] && mc[["usePX"]]) MHpars <- c(MHpars, mc[["name_xi"]])
      if (mc$strucA[["update.Q"]]) {
        MHpars <- c(MHpars, mc$strucA[["name_ext"]])
        if (mc$strucA$control[["adaptive"]]) {
          adapt <- add(adapt, bquote(mod[[.(k)]]$strucA$control$adapt(ar[[.(mc$strucA[["name_ext"]])]])))
        }
      }
      if (!is.null(mc[["AR1.inferred"]])) {
        MHpars <- c(MHpars, mc[["name_AR1"]])
        if (mc$AR1sampler$MH[["adaptive"]]) {
          adapt <- add(adapt, bquote(mod[[.(k)]]$AR1sampler$MH$adapt(ar[[.(mc[["name_AR1"]])]])))
        }
      }
      if (family[["family"]] == "gamma" && mc$control[["MHprop"]] == "LNRW") {
        MHpars <- c(MHpars, if (mc[["usePX"]]) mc[["name_sigma_raw"]] else mc[["name_sigma"]])
        adapt <- add(adapt, bquote(mod[[.(k)]]$adapt(ar)))
      }
    }
    if (mc[["type"]] != "mc_offset" && !(mc[["type"]] == "reg" && mc[["in_block"]])) {
      start <- add(start, bquote(p <- mod[[.(k)]]$start(p)))
      draw <- add(draw, bquote(p <- mod[[.(k)]]$draw(p)))
    }
  }

  if (length(block)) {
    for (k in seq_along(mbs)) {
      draw <- add(draw, bquote(p <- mbs[[.(k)]]$draw(p)))
      start <- add(start, bquote(p <- mbs[[.(k)]]$start(p)))
    }
  }

  if (do.linpred)
    draw <- add(draw, quote(p$linpred_ <- lin_predict(p, linpred)))

  draw <- add(draw, quote(p))  # return state p

  if (!control[["recompute.e"]] && !single.block) {
    # adding this sometimes gives bad starting values in case of single.block (no need anyway in that case)
    start <- add(start, quote(p$e_ <- compute_e(p)))
  }
  start <- add(start, quote(p))
  if (length(body(adapt)) <= 1L) rm(adapt)

  fam_llh_i <- family$make_llh_i(y)
  llh_i <- function(draws, i=seq_len(n)) {
    e_i <- compute_e_i(draws, i)
    fam_llh_i(draws, i, e_i)
  }

  # log-likelihood function for optimization: function of a vector instead of list
  if (!has.bart) {
    llh_opt <- function(x) {
      p <- vec2list(x)
      p$e_ <- compute_e(p)
      if (family[["family"]] == "gaussian") {
        SSR <- dotprodC(p[["e_"]], Q_e(p))
        llh(p, SSR)
      } else {
        llh(p)
      }
    }
  }

  # remove quantities no longer needed
  rm(k, mc, types, data)

  # return the function environment, including draw, rprior, start functions
  environment()
}


#' Set computational options for the sampling algorithms
#'
#' @export
#' @param add.outer.R whether to add the outer product of a constraint matrix to the
#'  conditional posterior precision matrix of coefficients sampled in a block. This is used
#'  to resolve singularity due to intrinsic GMRF components.
#'  By default, \code{add.outer.R=NULL}, a simple heuristic is used to decide whether
#'  to add the outer product of possibly a submatrix of the constraint matrix.
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
#' @param expanded.cMVN.sampler whether an expanded linear system including dual variables is used
#'  for equality constrained multivariate normal sampling. If set to \code{TRUE} this may
#'  improve the performance of the blocked Gibbs sampler, especially in case of a large number of equality
#'  constraints, typically GMRF identifiability constraints.
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
#' @param block.V if \code{TRUE}, the default, all coefficients of \code{reg} and \code{gen} components
#'  in a variance model formula are sampled in a single block. Alternatively, a list of
#'  character vectors with names of model components whose coefficients should be sampled together in blocks.
#' @param auto.order.block whether Gibbs blocks should be ordered automatically in such a
#'  way that those with the most sparse design matrices come first. This way of ordering
#'  can make Cholesky updates more efficient.
#' @param chol.control options for Cholesky decomposition, see \code{\link{chol_control}}.
#' @param max.size.cps.template maximum allowed size in MB of the sparse matrix serving as a 
#'  template for the sparse symmetric crossproduct X'QX of a dgCMatrix X, where Q is a diagonal
#'  matrix subject to change.
#' @param PG.approx whether Polya-Gamma draws for logistic binomial models are
#'  approximated by a hybrid gamma convolution approach. If not, \code{BayesLogit::rpg}
#'  is used, which is exact for some values of the shape parameter.
#' @param PG.approx.m if \code{PG.approx=TRUE}, the number of explicit gamma draws in the
#'  sum-of-gammas representation of the Polya-Gamma distribution. The remainder (infinite)
#'  convolution is approximated by a single moment-matching gamma draw. Special values are:
#'  \code{-2L} for a default choice depending on the value of the shape parameter
#'  balancing performance and accuracy, \code{-1L} for a moment-matching normal approximation,
#'  and \code{0L} for a moment-matching gamma approximation.
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
                            recompute.e=TRUE,
                            expanded.cMVN.sampler=FALSE, CG=NULL,
                            block=TRUE, block.V=TRUE, auto.order.block=TRUE,
                            chol.control=chol_control(),
                            max.size.cps.template=100,
                            PG.approx=TRUE, PG.approx.m=-2L) {
  list(add.outer.R=add.outer.R, add.eps.I=add.eps.I, eps=eps,
       recompute.e=recompute.e,
       expanded.cMVN.sampler=expanded.cMVN.sampler, CG=CG,
       block=block, block.V=block.V, auto.order.block=auto.order.block,
       chol.control = chol.control,
       max.size.cps.template=max.size.cps.template,
       PG.approx=PG.approx, PG.approx.m=PG.approx.m
  )
}

check_sampler_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- sampler_control()
  w <- which(!(names(control) %in% names(defaults)))
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
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
  if (is.logical(control[["block.V"]])) {
    if (length(control[["block.V"]]) != 1L) stop("unexpected input for 'block.V'")
  } else {
    if (!is.list(control[["block.V"]])) stop("'block.V' should be either a scalar logical or a list of variance model component name vectors")
    if (!length(control[["block.V"]])) stop("'block.V' must contain at least one character vector")
  }
  control$chol.control <- check_chol_control(control[["chol.control"]])
  if (isTRUE(control[["CG"]])) {
    control$CG <- CG_control()
  } else if (isFALSE(control[["CG"]])) {
    control$CG <- NULL
  } else if (!is.null(control[["CG"]])) {
    control$CG <- check_CG_control(control[["CG"]])
  }
  if (control[["expanded.cMVN.sampler"]] && is.list(control[["CG"]])) stop("'expanded.cMVN.sampler' and 'CG' cannot currently be combined")
  control
}
