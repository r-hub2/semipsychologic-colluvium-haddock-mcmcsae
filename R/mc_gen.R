#' Create a model component object for a generic random effects component in the linear predictor
#'
#' This function is intended to be used on the right hand side of the \code{formula} argument to
#' \code{\link{create_sampler}} or \code{\link{generate_data}}.
#'
#' @export
#' @param formula a model formula specifying the effects that vary over the levels of
#'  the factor variable(s) specified by argument \code{factor}. Defaults to \code{~1},
#'  corresponding to random intercepts. If \code{X} is specified \code{formula} is ignored.
#'  Variable names are looked up in the data frame passed as \code{data} argument to
#'  \code{\link{create_sampler}} or \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param factor a formula with factors by which the effects specified in the \code{formula}
#'  argument vary. Often only one such factor is needed but multiple factors are allowed so that
#'  interaction terms can be modeled conveniently. The formula must take the form
#'  \code{ ~ f1(fac1, ...) * f2(fac2, ...) ...}, where
#'  \code{fac1, fac2} are factor variables and \code{f1, f2} determine the
#'  correlation structure assumed between levels of each factor, and the \code{...} indicate
#'  that for some correlation types further arguments can be passed. Correlation structures
#'  currently supported include \code{iid} for independent identically distributed effects,
#'  \code{RW1} and \code{RW2} for random walks of first or second order over the factor levels,
#'  \code{AR1} for first-order autoregressive effects, \code{season} for seasonal effects,
#'  \code{spatial} for spatial (CAR) effects and \code{custom} for supplying a custom
#'  precision matrix corresponding to the levels of the factor. For further details about
#'  the correlation structures, and further arguments that can be passed, see
#'  \code{\link{correlation}}. Argument \code{factor} is ignored if \code{X} is specified.
#'  The factor variables are looked up in the data frame passed as \code{data} argument to
#'  \code{\link{create_sampler}} or \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param remove.redundant whether redundant columns should be removed from the model matrix
#'  associated with \code{formula}. Default is \code{FALSE}.
#' @param drop.empty.levels whether to remove factor levels without observations.
#' @param X A (possibly sparse) design matrix. If \code{X} is specified, \code{formula} and \code{factor}
#'  are only used to derive the random effects' structured precision matrix.
#' @param var the (co)variance structure among the varying effects defined by \code{formula}
#'  over the levels of the factors defined by \code{factor}.
#'  The default is \code{"unstructured"}, meaning that a full covariance matrix
#'  parametrization is used. For uncorrelated effects with unequal variances use
#'  \code{var="diagonal"}. For uncorrelated effects with equal variances use \code{var="scalar"}.
#'  In the case of a single varying effect there is no difference between these choices.
#' @param prior the prior specification for the variance parameters of the random effects.
#'  These can currently be specified by a call to \code{\link{pr_invwishart}} in case
#'  \code{var="unstructured"} or by a call to \code{\link{pr_invchisq}} otherwise.
#'  See the documentation of those prior specification functions for more details.
#' @param Q0 precision matrix associated with \code{formula}. This can only be used in
#'  combination with \code{var="scalar"}.
#' @param PX whether parameter expansion should be used. Default is \code{TRUE}, which
#'  applies parameter expansion with default options. The only exception is that for
#'  gamma sampling distributions the default is \code{FALSE}, i.e. no parameter expansion.
#'  Alternative options can be specified
#'  by supplying a list with one or more of the following components:
#'  \describe{
#'    \item{prior}{prior for the multiplicative expansion parameter. Defaults to a
#'      normal prior with mean 0 and standard deviation 1, unless the sampling
#'      distribution is gamma in which case the default is a Multivariate Log
#'      inverse Gamma prior. The default parameters can be changed using functions
#'      \code{\link{pr_normal}} or \code{\link{pr_MLiG}}.}
#     \item{mu0}{location (vector) parameter for parameter expansion. By default \code{0}.}
#     \item{Q0}{precision (matrix) parameter for parameter expansion. Default is the identity matrix.}
#'    \item{vector}{whether a redundant multiplicative expansion parameter is used for each varying effect
#'      specified by \code{formula}. The default is \code{TRUE} except when \code{var="scalar"}.
#'      If \code{FALSE} a single redundant multiplicative parameter is used.}
#'    \item{data.scale}{whether the data level scale is used as a variance factor for the expansion
#'      parameters. Default is \code{TRUE}.}
#     \item{sparse}{UNDOCUMENTED}
#'  }
#' @param priorA prior distribution for scale factors at the variance scale associated with \code{QA}.
#'  In case of IGMRF models the scale factors correspond to the innovations.
#'  The default \code{NULL} means not to use any local scale factors.
#'  A prior can currently be specified using \code{\link{pr_invchisq}} or \code{\link{pr_exp}}.
#' @param strucA this option can be used to modify the default structure encoded by
#'  \code{factor} to a 'bym2' or 'leroux' structure. See \code{\link{GMRF_structure}}
#'  for details.
#' @param R0 an optional equality restriction matrix acting on the coefficients defined by \code{formula}, for each
#'  level defined by \code{factor}. If c is the number of restrictions, \code{R0} is a
#'  q0 x c matrix where q0 is the number of columns of the design matrix derived
#'  from \code{formula}. Together with \code{RA} it defines the set of equality constraints
#'  to be imposed on the vector of coefficients. Only allowed in combination with \code{var="scalar"}.
#' @param RA an optional equality restriction matrix acting on the coefficients defined by \code{factor},
#'  for each effect defined by \code{formula}. If c is the number of restrictions, \code{RA} is a
#'  l x c matrix where l is the number of levels defined by \code{factor}.
#'  Together with \code{R0} this defines the set of equality constraints to be imposed on the vector
#'  of coefficients.
#'  If \code{constr=TRUE}, additional constraints are imposed, corresponding to the
#'  null-vectors of the singular precision matrix in case of an intrinsic Gaussian Markov Random Field.
#' @param constr whether constraints corresponding to the null-vectors of the precision matrix
#'  are to be imposed on the vector of coefficients. By default this is \code{TRUE} for
#'  improper or intrinsic GMRF model components, i.e. components with a singular precision matrix
#'  such as random walks or CAR spatial components.
#' @param S0 an optional inequality restriction matrix acting on the coefficients defined by \code{formula}, for each
#'  level defined by \code{factor}. If c is the number of restrictions, \code{S0} is a
#'  q0 x c matrix where q0 is the number of columns of the design matrix derived
#'  from \code{formula}. Together with \code{SA} it defines the set of inequality constraints
#'  to be imposed on the vector of coefficients.
#   TODO does this work, also for var != "scalar"?
#' @param SA an optional inequality restriction matrix acting on the coefficients defined by \code{factor},
#'  for each effect defined by \code{formula}. If c is the number of restrictions, \code{SA} is a
#'  l x c matrix where l is the number of levels defined by \code{factor}.
#'  Together with \code{S0} this defines the set of constraints to be imposed on the vector
#'  of coefficients.
# TODO R,S instead of R0,RA and S0,SA, and rhs r,s
#' @param formula.gl a formula of the form \code{~ glreg(...)} for group-level predictors
#'  around which the random effect component is hierarchically centered.
#'  See \code{\link{glreg}} for details.
#' @param a only used in case the effects are MLiG distributed, as
#'  assumed in case of a gamma sampling distribution, or for
#'  gaussian variance modelling. In those cases \code{a} controls how close
#'  the effects' prior is to a normal prior, see \code{\link{pr_MLiG}}.
#' @param name the name of the model component. This name is used in the output of the MCMC simulation
#'  function \code{\link{MCMCsim}}. By default the name will be 'gen' with the number of the model term attached.
#' @param sparse whether the model matrix associated with \code{formula} should be sparse.
#'  The default is based on a simple heuristic based on storage size.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @param control a list with further computational options. These options can
#'  be specified using function \code{\link{gen_control}}.
#' @returns An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Intended
#'  for internal use by other package functions.
#' @references
#'  J. Besag and C. Kooperberg (1995).
#'    On Conditional and Intrinsic Autoregression.
#'    Biometrika 82(4), 733-746.
#'
#'  C.M. Carvalho, N.G. Polson and J.G. Scott (2010).
#'    The horseshoe estimator for sparse signals.
#'    Biometrika 97(2), 465-480.
#'
#'  L. Fahrmeir, T. Kneib and S. Lang (2004).
#'    Penalized Structured Additive Regression for Space-Time Data:
#'    a Bayesian Perspective.
#'    Statistica Sinica 14, 731-761.
#'
#'  A. Gelman (2006).
#'    Prior distributions for variance parameters in hierarchical models.
#'    Bayesian Analysis 1(3), 515-533.
#'
#'  A. Gelman, D.A. Van Dyk, Z. Huang and W.J. Boscardin (2008).
#'    Using Redundant Parameterizations to Fit Hierarchical Models.
#'    Journal of Computational and Graphical Statistics 17(1), 95-122.
#'
#'  T. Park and G. Casella (2008).
#'    The Bayesian Lasso.
#'    Journal of the American Statistical Association 103(482), 681-686.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
gen <- function(formula = ~ 1, factor=NULL,
                remove.redundant=FALSE, drop.empty.levels=FALSE, X=NULL,
                var=NULL, prior=NULL, Q0=NULL, PX=NULL,
                priorA=NULL,
                strucA=GMRF_structure(),
                R0=NULL, RA=NULL, constr=NULL, S0=NULL, SA=NULL,
                formula.gl=NULL, a=1000,
                name="", sparse=NULL, control=gen_control(), debug=FALSE) {

  e <- sys.frame(-2L)
  type <- "gen"  # for generic (random effects)
  if (name == "") stop("missing model component name")

  if (e$family$family == "gamma") {
    modus <- "gamma"       # model for log(mean) of gamma
  } else if (any(name == names(e$Vmod))) {
    if (e$family$family == "gaussian_gamma")
      modus <- "vargamma"  # model for log(var) of gaussian and log(mean) of gamma
    else
      modus <- "var"       # model for log(var) of gaussian
  } else {
    modus <- "regular"     # model for mean of gaussian/binomial/...
  }

  # check supplied var component; if NULL leave it to be filled in later
  if (!is.null(var)) var <- match.arg(var, c("unstructured", "diagonal", "scalar"))
  if (!is.null(factor) && !inherits(factor, "formula")) stop("element 'factor' of a model component must be a formula")
  if (!is.environment(strucA)) stop("'strucA' must be an environment created by GMRF_structure")

  if (e$family$family == "multinomial") {
    edat.formula <- new.env(parent = environment(formula))
    environment(formula) <- edat.formula
    edat.factor <- new.env(parent = environment(factor))
    environment(factor) <- edat.factor
  }

  varnames <- factor.cols.removed <- NULL
  if (is.null(X)) {
    if (e$family$family == "multinomial") {
      for (k in seq_len(e[["Km1"]])) {
        edat.factor$cat_ <- edat.formula$cat_ <- factor(rep.int(e$cats[k], e[["n0"]]), levels=e$cats[-length(e$cats)])
        if (k == 1L) {
          X0 <- model_matrix(formula, e[["data"]], sparse=sparse)
          XA <- compute_XA(factor, e[["data"]])
        } else {
          X0 <- rbind(X0, model_matrix(formula, e[["data"]], sparse=sparse))
          if (!is.null(XA)) XA <- rbind(XA, compute_XA(factor, e[["data"]]))
        }
      }
      edat.factor$cat_ <- edat.formula$cat_ <- NULL
    } else {
      X0 <- model_matrix(formula, e[["data"]], sparse=sparse)
      XA <- compute_XA(factor, e[["data"]])
    }
    if (remove.redundant) X0 <- remove_redundancy(X0)
    varnames <- colnames(X0)
    if (!is.null(XA)) {
      if (drop.empty.levels) {
        factor.cols.removed <- which(zero_col(XA))
        if (length(factor.cols.removed)) XA <- XA[, -factor.cols.removed, drop=FALSE]
      }
      if (strucA$type == "bym2") XA <- cbind(XA, zeroMatrix(e[["n"]], ncol(XA)))
      X <- combine_X0_XA(X0, XA)
    } else {
      X <- X0
    }
    rm(X0, XA)
  } else {
    if (is.null(colnames(X))) colnames(X) <- seq_len(ncol(X))
  }
  if (nrow(X) != e[["n"]]) stop("design matrix with incompatible number of rows")
  e$coef.names[[name]] <- colnames(X)
  X <- economizeMatrix(X, sparse=sparse, strip.names=TRUE, check=TRUE)
  q <- ncol(X)
  in_block <- any(name == unlst(e[["block"]])) || any(name == unlst(e[["block.V"]]))

  self <- environment()

  if (e$family$family == "multinomial") {
    edat.factor$cat_ <- e$cats[-length(e$cats)]  # K-1 categories; TODO for general GMRF prediction need all K categories?
  }
  info <- get_factor_info(factor, e[["data"]])
  is.proper <- all(info$types %in% c("iid", "AR1")) || strucA$type == "leroux"
  # by default, IGMRF constraints are imposed, unless the GMRF is proper
  # NB constr currently only refers to QA, not Q0
  if (is.null(constr)) constr <- !is.proper
  # fastGMRFprior currently only for IGMRF without user-defined constraints
  #   and currently only used for IGMRFs
  #   and not possible for spatial as typically n_edges > n_vertices
  #   and only if any custom factor has (incidence) matrix D specified
  fastGMRFprior <- !is.null(info) && is.null(priorA) && modus == "regular" &&
    is.null(RA) && is.null(R0) && is.null(SA) && is.null(S0) &&
    strucA$type == "default" &&
    !all(info$types %in% c("iid", "AR1")) && all(info$types != "spatial") &&
    all(b_apply(info$factors[info$types == "custom"], function(x) !is.null(x[["D"]]))) &&
    all(!b_apply(info$factors, function(x) isTRUE(x[["circular"]])))
  GMRFmats <- compute_GMRF_matrices(info, e[["data"]],
    D=fastGMRFprior || !is.null(priorA) || !is.null(e$control$CG) || e$control$expanded.cMVN.sampler,
    Q=!(fastGMRFprior || !is.null(priorA)),
    R=constr, sparse=if (in_block) TRUE else NULL,
    cols2remove=factor.cols.removed, scale.precision=strucA$scale.precision, drop.zeros=TRUE
  )

  AR1.inferred <- which(info$types == "AR1" & !b_apply(info$extra, is.null))
  if (length(AR1.inferred) == 1L) {
    if (!is.null(formula.gl)) stop("AR1 with non-trivial parameter prior not supported in combination with group-level effects")
    if (strucA$type != "default") stop("AR1 with non-trivial parameter prior not supported in combination with GMRF extension")
    if (!is.null(info$factors[[AR1.inferred]][["w"]])) stop("AR1 with non-trivial parameter prior not supported in combination with AR1 weights")
  } else if (length(AR1.inferred) == 2L) {
    stop("only one AR1 factor with a parameter to be inferred allowed")
  } else {
    AR1.inferred <- NULL
  }

  if (fastGMRFprior || !is.null(priorA) || !is.null(e$control$CG) || e$control$expanded.cMVN.sampler) {
    if (is.null(AR1.inferred)) {
      DA <- GMRFmats[["D"]]  # lD x l incidence matrix DA
      l <- ncol(DA)
      lD <- nrow(DA)
    } else {
      DA.template <- DA_AR1_template(info, GMRFmats[["D"]], AR1.inferred)
      l <- ncol(DA.template$DA0.5)
      lD <- nrow(DA.template$DA0.5)
    }
    if (is.null(priorA)) {
      if (is.null(AR1.inferred)) {
        QA <- economizeMatrix(crossprod(DA), sparse=if (in_block) TRUE else NULL, symmetric=TRUE, check=TRUE)
      } else {
        QA.template <- QA_AR1_template(info,
          economizeMatrix(crossprod(DA.template$DA0.5), sparse=TRUE, symmetric=TRUE, check=TRUE),
          AR1.inferred
        )
      }
    }
  } else {  # compute precision matrix QA only
    if (is.null(AR1.inferred)) {
      QA <- economizeMatrix(GMRFmats[["Q"]], sparse=if (in_block) TRUE else NULL, symmetric=TRUE, check=TRUE)
      l <- nrow(QA)
    } else {
      QA.template <- QA_AR1_template(info, GMRFmats[["Q"]], AR1.inferred)
      l <- ncol(QA.template$QA0.5)
    }
  }
  if (strucA$type == "bym2") l <- 2L * l
  if (q %% l != 0L) stop("incompatible dimensions of design and precision matrices")
  q0 <- q %/% l  # --> q = q0 * l

  if (!is.null(priorA)) {  # scaled precision matrix QA = DA' QD DA with QD modeled
    name_omega <- paste0(name, "_omega")
    if (strucA$type != "default") stop("not implemented: GMRF extension in combination with non-normal random effects")
    switch(priorA$type,
      invchisq = {
        df.data.omega <- q0  # TODO is this always the correct value?
      },
      exp = {},
      stop("priorA argument expects a prior specified with pr_invchisq or pr_exp")
    )
    priorA$init(lD, !e$prior.only)
  }

  if (q0 == 1L) {
    var <- "scalar"  # single varying effect
  } else if (is.null(var)) {
    if (l == 1L)
      var <- "diagonal"  # single group, default ARD prior
    else if (q0 <= 500L)
      var <- "unstructured"
    else {
      warn("model component '", name, "': default variance structure changed
        to 'diagonal' because of the large number ", q0, " of varying effects.
        If you instead intend to model a full covariance matrix among all varying effects,
        use var='unstructured', though it may be very slow.")
      var <- "diagonal"
    }
  }

  if (is.null(PX)) PX <- modus == "regular"
  usePX <- !isFALSE(PX)
  if (usePX) {
    if (modus == "regular")
      PX.defaults <- list(prior=pr_normal(mean=0, precision=1), data.scale = !e$sigma.fixed, vector=var != "scalar", sparse=NULL)
    else
      PX.defaults <- list(prior=pr_MLiG(mean=0, precision=1), data.scale=FALSE, vector=FALSE, sparse=FALSE)
    if (is.list(PX)) {
      if (!all(names(PX) %in% names(PX.defaults))) stop("invalid 'PX' options list")
      PX <- modifyList(PX.defaults, PX)
      if (modus != "regular" && PX$data.scale) {
        warn("PX data.scale has been set to FALSE")
        PX$data.scale <- FALSE
      }
    } else {
      if (is.logical(PX) && length(PX) == 1L)
        PX <- PX.defaults
      else
        stop("wrong input for 'PX'")
    }
    rm(PX.defaults)

    if (is.null(PX$sparse)) PX$sparse <- q0 > 10L
    if (q0 == 1L) PX$vector <- FALSE
    if (!PX$vector) PX$sparse <- FALSE
    PX$dim <- if (PX$vector) q0 else 1L
    switch(PX$prior$type,
      normal = {
        if (modus != "regular") stop("please use 'pr_MLiG' to specify a (conjugate) prior for gamma family or variance modelling")
        PX$prior$init(n=PX$dim, sparse=PX$sparse, sigma=PX$data.scale)
        if (PX$vector)
          base_tcrossprod <- base::tcrossprod  # faster access (Matrix promotes tcrossprod to S4 generic)
        PX_Q0 <- PX$prior$precision
        if (length(PX$prior$mean) == 1L)
          PX_Qmu0 <- PX_Q0 %m*v% rep.int(PX$prior$mean, PX$prior$n)
        else
          PX_Qmu0 <- PX_Q0 %m*v% PX$prior$mean
      },
      MLiG = {
        if (modus == "regular") stop("MLiG prior only available in combination with gamma family or variance modelling")
        PX$prior$init(n=PX$dim)
      }
    )
    name_xi <- paste0(name, "_xi_")  # trailing '_' --> not stored by MCMCsim even if store.all=TRUE
  }  # END if (usePX)

  # construct equality restriction matrix
  # TODO if R is supplied, it is assumed to be the full constraint matrix...
  #      in that case the derivation from R0 and RA can be skipped
  if (!is.null(R0)) {
    if (var != "scalar") stop("R0 constraint matrix only allowed if var='scalar'")
    if (nrow(R0) != q0) stop("incompatible matrix R0")
  }

  # if a matrix 'RA' is provided by the user, use these constraints in addition to possible GMRF constraints
  if (is.null(RA)) {
    RA <- GMRFmats[["R"]]
  } else {
    if (nrow(RA) != l) stop("incompatible matrix RA")
    if (!is.null(GMRFmats[["R"]])) {
      RA <- economizeMatrix(cbind(RA, GMRFmats[["R"]]), allow.tabMatrix=FALSE, check=TRUE)
      RA <- remove_redundancy(RA)  # combination of user-supplied and IGMRF constraints may contain redundant columns
    }
  }

  # RA is needed here for scale.precision, and it is expanded in case of bym2
  strucA <- create_GMRF_structure(strucA, self, e$prior.only)
  
  R <- switch(paste0(is.null(RA), is.null(R0)),
    TRUETRUE = NULL,
    TRUEFALSE = kronecker(CdiagU(l), R0),
    FALSETRUE = kronecker(RA, CdiagU(q0)),
    FALSEFALSE = cbind(kronecker(CdiagU(l), R0), kronecker(RA, CdiagU(q0)))
  )
  if (!is.null(R)) {
    # tabMatrix doesn't work yet because of lack of solve method
    R <- economizeMatrix(R, allow.tabMatrix=FALSE, check=TRUE)
  }
  if (!is.null(R0) && !is.null(RA)) R <- remove_redundancy(R)
  rm(R0, GMRFmats)

  # construct inequality restriction matrix (TODO S instead of S0, SA)
  if (!is.null(S0) && nrow(S0) != q0) stop("incompatible matrix S0")
  if (!is.null(SA) && nrow(SA) != l) stop("incompatible matrix SA")
  S <- switch(paste0(is.null(SA), is.null(S0)),
    TRUETRUE = NULL,
    TRUEFALSE = kronecker(CdiagU(l), S0),
    FALSETRUE = kronecker(SA, CdiagU(q0)),
    FALSEFALSE = cbind(kronecker(CdiagU(l), S0), kronecker(SA, CdiagU(q0)))
  )
  rm(S0, SA)

  if (is.null(prior)) {
    prior <- switch(var,
      unstructured = pr_invwishart(df=q0+1, scale=diag(q0)),
      diagonal=, scalar = pr_invchisq(if (usePX) 1 else 1e-3, 1)
    )
  }
  if (all(prior$type != if (var == "unstructured") "invwishart" else c("invchisq", "exp"))) stop("unsupported prior")
  if (is.list(prior$df)) stop("not supported: modeled df for random effect variance prior")
  switch(prior$type,
    invwishart = prior$init(q0),
    invchisq = {
      if (var == "scalar") {
        prior$init(1L, !e$prior.only)
        df.data <- if (is.null(R)) q else q - ncol(R)
      } else {  # var "diagonal"
        prior$init(q0, !e$prior.only)
        # df based on number of unconstrained coefficients
        df.data <- if (is.null(RA)) l else l - ncol(RA)
      }
    },
    exp = {
      if (var == "scalar") {
        prior$init(1L, !e$prior.only)
        ppar <- if (is.null(R)) q else q - ncol(R)
      } else {  # var "diagonal"
        prior$init(q0, !e$prior.only)
        ppar <- if (is.null(RA)) l else l - ncol(RA)
      }
      ppar <- 1 - ppar/2
    }
  )

  if (var == "scalar") {
    # also check given precision matrix Q0, if any
    if (!is.null(Q0)) {
      Q0 <- economizeMatrix(Q0, symmetric=TRUE, vec.diag=TRUE, check=TRUE)
      if (is.vector(Q0)) {
        if (all(length(Q0) != c(1L, q0))) stop("incompatible 'Q0'")
      } else {
        if (!identical(dim(Q0), c(q0, q0))) stop("incompatible 'Q0'")
      }
      if (is.vector(Q0) && all(Q0 == 1)) Q0 <- NULL  # identity Q0, can be ignored
    }
  } else {
    if (var == "unstructured") {
      # TODO add rprior and draw methods to pr_invwishart
      df <- prior$df + l
      if (!is.null(RA)) df <- df - ncol(RA)
    }
    if (!is.null(Q0)) warn("'Q0' ignored")
  }
  rm(RA)

  gl <- !is.null(formula.gl)
  if (gl) {  # group-level covariates
    if (!inherits(formula.gl, "formula")) stop("'formula.gl' must be a formula")
    vars <- as.list(attr(terms(formula.gl), "variables"))[-1L]
    if (length(vars) != 1L || as.character(vars[[1L]][[1L]]) != "glreg") stop("'formula.gl' should be a formula with a single 'glreg' component")
    name_gl <- paste0(name, "_gl")
    glp <- vars[[1L]]
    glp$name <- name_gl
    rm(vars)
  }

  if (e$modeled.Q) {
    if (gl) {
      # ensure that XX is always sparse because we need a sparse block diagonal template in this case
      X <- economizeMatrix(X, sparse=TRUE)  # --> crossprod_sym(X, Q) also sparse
    }
    XX <- crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]]))
  } else {
    XX <- economizeMatrix(crossprod_sym(X, e[["Q0"]]), symmetric=TRUE, drop.zeros=TRUE)
  }

  # for both memory and speed efficiency define unit_Q case (random intercept and random slope with scalar var components, no constraints)
  unit_Q <- is.null(priorA) && is.null(AR1.inferred) && is_unit_diag(QA) && (var == "scalar" && is.null(Q0)) && is.null(R)
  if (modus != "regular") {
    if (!(unit_Q && q0 == 1L && !gl && strucA$type == "default" && is.null(S))) stop("gamma family and variance modelling can currently be combined with basic generic model components only")
  }

  generate_Qv <-
    if (var == "unstructured") {
      function() rWishart(1L, df=prior$df, Sigma=diag(q0))[,,1L]
    } else if (var == "scalar" && !is.null(Q0)) {
      function() scale_mat(Q0, runif(if (usePX && PX$vector) q0 else 1L, 0.25, 0.75))
    } else if (var == "diagonal" || (usePX && PX$vector)) {
      function() runif(q0, 0.25, 0.75)
    } else {  # scalar var, scalar or no PX
      function() runif(1L, 0.25, 0.75)
    }

  if (gl) {  # group-level covariates
    glp <- eval(glp)
    sparse_template(self, update.XX=e$modeled.Q)
    if (!in_block && !e$prior.only) glp$R <- NULL
  } else if (modus == "regular") {
    sparse_template(self, update.XX=e$modeled.Q)
  } else {
    if (in_block) {
      Q <- Cdiag(rep.int(1/a, q))  # to construct QT template in mc_block
    } else {
      if (modus == "vargamma")
        cholHH <- build_chol(2 * XX  + (1/a) * runif(1L, 0.9, 1.1) * CdiagU(q))
      else
        cholHH <- build_chol(XX + (1/a) * runif(1L, 0.9, 1.1) * CdiagU(q))
    }
  }

  if (!is.null(AR1.inferred)) AR1sampler <- AR1_sampler(self)

  name_sigma <- paste0(name, "_sigma")
  if (var == "unstructured") name_rho <- paste0(name, "_rho")

  if (q0 > 1L) {  # add labels for sigma, rho, xi to e$coef.names
    if (var != "scalar" || (usePX && PX$vector)) {
      e$coef.names[[name_sigma]] <- varnames
      if (var == "unstructured")
        e$coef.names[[name_rho]] <- outer(varnames, varnames, FUN=paste, sep=":")[upper.tri(diag(q0))]
      if (usePX && PX$vector)
        e$coef.names[[name_xi]] <- varnames
    }
    if (gl) {
      if (is.null(e$coef.names[[name_gl]]))
        e$coef.names[[name_gl]] <- varnames
      else
        e$coef.names[[name_gl]] <- as.vector(outer(varnames, e$coef.names[[name_gl]], FUN=paste, sep=":"))
    }
  }
  rm(varnames)

  if (modus == "var" || modus == "vargamma") {
    if (is_ind_matrix(X) && q < e[["n"]])
      compute_Qfactor <- function(p) X %m*v% exp(-p[[name]])
    else
      compute_Qfactor <- function(p) exp(X %m*v% (-p[[name]]))
  }
  lp <- function(p) X %m*v% p[[name]]
  lp_update <- function(x, plus=TRUE, p) mv_update(x, plus, X, p[[name]])
  # draws_linpred method acts on (subset of) mcdraws object, used in fitted() and pointwise log-likelihood llh_i functions
  draws_linpred <- function(obj, units=NULL, chains=NULL, draws=NULL, matrix=FALSE) {
    if (is.null(units)) Xi <- X else Xi <- X[units, , drop=FALSE]
    if (matrix) {
      out <- tcrossprod(as.matrix.dc(get_from(obj[[name]], chains=chains, draws=draws), colnames=FALSE), Xi)
    } else {
      out <- list()
      for (ch in seq_along(chains))
        out[[ch]] <- tcrossprod(get_from(obj[[name]], chains=chains[ch], draws=draws)[[1L]], Xi)
    }
    out
  }

  make_predict <- function(newdata=NULL, Xnew=NULL, verbose=TRUE) {
    if (modus == "vargamma") stop("prediction for 'gaussian_gamma' family not supported")
    linpred <- function(p) Xnew %m*v% p[[name]]
    linpred_update <- function(x, plus=TRUE, p) mv_update(x, plus, Xnew, p[[name]])
    if (is.null(newdata)) {
      if (is.null(Xnew)) {
        # in-sample prediction
        Xnew <- X
      } else {
        Xnew <- economizeMatrix(Xnew, strip.names=TRUE, check=TRUE)
        if (ncol(Xnew) != q) stop("wrong number of columns for Xnew matrix of component ", name)
      }
    } else {
      if (e$family$family == "multinomial") {
        for (k in seq_len(e[["Km1"]])) {
          edat.factor$cat_ <- edat.formula$cat_ <- factor(rep.int(e$cats[k], nrow(newdata)), levels=e$cats[-length(e$cats)])
          if (k == 1L) {
            X0 <- model_matrix(formula, data=newdata, sparse=sparse)
            XA <- compute_XA(factor, newdata)
          } else {
            X0 <- rbind(X0, model_matrix(formula, data=newdata, sparse=sparse))
            if (!is.null(XA)) XA <- rbind(XA, compute_XA(factor, newdata))
          }
        }
        edat.factor$cat_ <- edat.formula$cat_ <- NULL
      } else {
        X0 <- model_matrix(formula, newdata, sparse=sparse)
        XA <- compute_XA(factor, newdata)
      }
      if (!is.null(XA))
        Xnew <- combine_X0_XA(X0, XA)
      else
        Xnew <- X0
      rm(X0, XA)
      if (unit_Q) {
        # in this case we allow oos random effects and mixed cases by matching training and test levels
        m <- match(colnames(Xnew), e$coef.names[[name]])
        qnew <- ncol(Xnew)
        Xnew <- economizeMatrix(Xnew, sparse=sparse, strip.names=TRUE, check=TRUE)
        if (all(is.na(m))) {
          # all random effects in Xnew are new
          if (verbose) message("model component '", name, "': all categories in 'newdata' are out-of-sample categories")
          linpred <- function(p) Xnew %m*v% Crnorm(qnew, sd = p[[name_sigma]])
          linpred_update <- function(x, plus=TRUE, p) mv_update(x, plus, Xnew, Crnorm(qnew, sd = p[[name_sigma]]))
        } else if (anyNA(m)) {
          # both existing and new levels in newdata
          q_unmatched <- sum(is.na(m))
          if (verbose) message("model component '", name, "': ", q_unmatched, " out of ", qnew, " categories in 'newdata' are out-of-sample categories")
          I_matched <- which(!is.na(m))
          I_unmatched <- setdiff(seq_len(qnew), I_matched)
          I_coef <- m[!is.na(m)]
          # TODO next function in C++ (no need to initialize)
          linpred <- function(p) {
            v <- numeric(qnew)
            v[I_matched] <- p[[name]][I_coef]  # TODO distinguish case where I_coef = 1:q
            v[I_unmatched] <- Crnorm(q_unmatched, sd=p[[name_sigma]])
            Xnew %m*v% v
          }
          linpred_update <- function(x, plus=TRUE, p) {
            v <- numeric(qnew)
            v[I_matched] <- p[[name]][I_coef]  # TODO distinguish case where I_coef = 1:q
            v[I_unmatched] <- Crnorm(q_unmatched, sd=p[[name_sigma]])
            mv_update(x, plus, Xnew, v)
          }
        } else {
          # all columns of Xnew correspond to existing levels
          if (!identical(m, seq_len(q))) {
            I_coef <- m
            linpred <- function(p) Xnew %m*v% p[[name]][I_coef]
            linpred_update <- function(x, plus=TRUE, p) mv_update(x, plus, Xnew, p[[name]][I_coef])
          }
        }
        rm(m)
      } else {
        # generic case: here we expect the same levels in training and test sets
        # for prediction do not (automatically) remove redundant columns!
        if (remove.redundant || drop.empty.levels) {
          Xnew <- Xnew[, e$coef.names[[name]], drop=FALSE]
        }
        if (ncol(Xnew) != q) stop("'newdata' yields ", ncol(Xnew), " predictor column(s) for model term '", name, "' versus ", q, " originally")
        Xnew <- economizeMatrix(Xnew, sparse=sparse, strip.names=TRUE, check=TRUE)
      }
    }
    rm(newdata, verbose)
    environment()
  }

  # BEGIN rprior function
  # draw from prior; draws are independent and stored in p
  if (usePX) {
    rprior <- function(p) {
      xi <- PX$prior$rprior(p)
      p[[name_xi]] <- xi
      inv_xi <- 1/xi
    }
  } else {
    rprior <- function(p) {inv_xi <- 1}
  }

  if (var == "unstructured") {
    if (is.list(prior$scale)) {
      psi0 <- prior$scale$df / prior$scale$scale
      if (prior$scale$common)
        rprior <- add(rprior, quote(psi0 <- diag(rep.int(rchisq_scaled(1L, prior$scale$df, psi=psi0), q0))))
      else
        rprior <- add(rprior, quote(psi0 <- diag(rchisq_scaled(q0, prior$scale$df, psi=psi0))))
    } else {
      psi0 <- prior$scale  # matrix
    }
    rprior <- add(rprior, quote(Qraw <- rWishart(1L, df=prior[["df"]], Sigma=inverseSPD(psi0))[,,1L]))
    if (usePX && PX$vector)
      rprior <- add(rprior, quote(inv_xi_qform <- base_tcrossprod(inv_xi)))
    else
      rprior <- add(rprior, quote(inv_xi_qform <- inv_xi^2))
    rprior <- add(rprior, quote(Qv <- Qraw * inv_xi_qform))  # use Qv to distinguish from data-level precision Q
    # convert precision matrix to standard errors and correlations
    names_se_cor <- c(name_sigma, name_rho)
    rprior <- add(rprior, quote(p[names_se_cor] <- prec2se_cor(Qv)))
  } else {
    rprior <- rprior |>
      add(quote(Qraw <- 1 / prior$rprior())) |>
      add(quote(Qv <- Qraw * inv_xi^2)) |>
      add(bquote(p[[.(name_sigma)]] <- sqrt(1/Qv)))
  }
  if (strucA$update.Q) {
    rprior <- add(rprior, bquote(p[[.(strucA$name_ext)]] <- strucA$rprior()))
    rprior <- add(rprior, bquote(QA <- strucA$update_Q(QA, p[[.(strucA$name_ext)]])))
  }
  if (!is.null(AR1.inferred)) {
    name_AR1 <- paste0(name, "_AR1")
    rprior <- add(rprior, bquote(p[[.(name_AR1)]] <- info$extra[[AR1.inferred]]$prior$rprior()))
  }
  if (gl) {  # draw group-level predictor coefficients
    rprior <- add(rprior, bquote(p[[.(name_gl)]] <- glp$prior$rprior(p)))
  }

  if (!is.null(priorA)) {
    switch(priorA$type,
      invchisq = {
        if (is.list(priorA$df)) {
          name_df <- paste0(name, "_df")
          rprior <- rprior |>
            add(bquote(p[[.(name_df)]] <- priorA$rprior_df())) |>
            add(bquote(p[[.(name_omega)]] <- priorA$rprior(p[[.(name_df)]])))
        } else {
          rprior <- add(rprior, bquote(p[[.(name_omega)]] <- priorA$rprior()))
        }
      },
      exp = {
        rprior <- add(rprior, bquote(p[[.(name_omega)]] <- priorA$rprior()))
      }
    )
    rprior <- add(rprior, bquote(QA <- crossprod_sym(DA, 1/p[[.(name_omega)]])))
  }

  # draw coefficient from prior
  if (unit_Q) {  # slightly faster alternative for random intercept models (or random slope with scalar variance)
    if (modus == "regular") {
      rprior <- add(rprior, bquote(p[[.(name)]] <- Crnorm(.(q), sd=p[[.(name_sigma)]])))
    } else {
      rprior <- add(rprior, bquote(p[[.(name)]] <- sqrt(a) * p[[.(name_sigma)]] * rMLiG(.(q), a, a)))
    }
  } else {
    if (fastGMRFprior) {
      rGMRFprior <- NULL  # to be created at the first call of rprior
      rprior <- rprior |>
        add(quote(if (is.null(rGMRFprior)) setup_priorGMRFsampler(self, Qv))) |>
        add(bquote(p[[.(name)]] <- rGMRFprior(Qv)))
    } else {
      if (q0 == 1L)
        rprior <- add(rprior, quote(Q <- Qv * QA))
      else
        rprior <- add(rprior, quote(Q <- kron_prod(QA, Qv)))
      rMVNprior <- NULL  # to be created at the first call of rprior
      rprior <- rprior |>
        add(quote(if (is.null(rMVNprior)) setup_priorMVNsampler(self, Q))) |>
        add(bquote(p[[.(name)]] <- rMVNprior(p, Q)))
    }
  }
  if (gl) {
    # add U alpha, where U = glp$X x I_q0 and alpha are the gl effects
    if (glp$p0 == 1L)
      rprior <- add(rprior, bquote(p[[.(name)]] <- p[[.(name)]] + glp[["X"]] %m*m% matrix(p[[.(name_gl)]], 1L)))
    else
      rprior <- add(rprior, quote(stop("TBI: include multiple group-level effect centering in generation of random effects")))
  }
  rprior <- add(rprior, quote(p))
  # END rprior function

  if (e$prior.only) return(self)  # no need for posterior draw function in this case

  if (!is.null(priorA) || is.list(prior$scale) || strucA$update.Q) {
    # store Qraw in state p --> no need to recompute at next MCMC iteration
    name_Qraw <- paste0(name, "_Qraw_")
  }

  # BEGIN draw function
  draw <- if (debug) function(p) {browser()} else function(p) {}
  if (in_block || !is.null(AR1.inferred)) {
    # store QA x Qv for use in block sampler, and/or AR1 parameter sampler
    name_Q <- paste0(name, "_Q_")  # trailing "_" --> only temporary storage
  }
  if (!is.null(e$control$CG) || e$control$expanded.cMVN.sampler ||
      (!is.null(AR1.inferred) && AR1sampler$MH$type != "TN")) {
    name_Qv <- paste0(name, "_Qv_")
  }
  if (!is.null(AR1.inferred)) {
    draw <- add(draw, bquote(phi <- p[[.(name_AR1)]]))
    draw <- add(draw, bquote(p[[.(name_AR1)]] <- AR1sampler$draw(phi, p)))
    if (fastGMRFprior || !is.null(priorA) || !is.null(e$control$CG) || e$control$expanded.cMVN.sampler) {
      draw <- add(draw, bquote(DA <- DA.template$update(p[[.(name_AR1)]])))
      if (is.null(priorA))
        draw <- add(draw, bquote(QA <- QA.template$update(p[[.(name_AR1)]])))
    } else {
      draw <- add(draw, bquote(QA <- QA.template$update(p[[.(name_AR1)]])))
    }
  }
  if (modus == "var" || modus == "vargamma") {
    if (!e$single.V.block) {
      if (is_ind_matrix(X) && q < e[["n"]])
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * (X %m*v% exp(p[[.(name)]]))))
      else
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(X %m*v% p[[.(name)]])))
    }
  } else {
    if (e$single.block && length(e$mod) == 1L) {
      # optimization in case of a single regression term, only in case of single mc (due to PX)
      draw <- add(draw, quote(p$e_ <- e$y_eff()))
    } else {
      if (e$e.is.res)
        draw <- add(draw, bquote(mv_update(p[["e_"]], plus=TRUE, X, p[[.(name)]])))
      else
        draw <- add(draw, bquote(mv_update(p[["e_"]], plus=FALSE, X, p[[.(name)]])))
    }
  }

  if (!e$e.is.res && e$single.block && length(e$mod) > 1L && modus == "regular") {
    # need to correct Q_e function; this should not be necessary if we draw all xi's in a single block!!
    if (e$family$link == "probit")
      draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p) - p[["e_"]])))
    else
      draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p) - p[["Q_"]] * p[["e_"]])))
  } else if (modus == "regular") {
    draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p))))
  } else {
    if (modus == "var" || modus == "vargamma") {  # variance modelling
      if (e$single.V.block)
        draw <- add(draw, bquote(vkappa <- .(if (e$sigma.fixed) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2))
      else
        draw <- add(draw, bquote(vkappa <- .(if (e$sigma.fixed) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2 * p[["Q_"]]))
    }
    if (modus == "gamma") {
      if (e$family$alpha.fixed) {
        alpha <- e$family$get_shape()
        if (e$single.block)
          kappa <- alpha * e[["y"]]
        else {
          kappa0 <- alpha * e[["y"]]
          draw <- add(draw, quote(kappa <- kappa0 * exp(-p[["e_"]])))
        }
      } else {
        draw <- add(draw, quote(alpha <- e$family$get_shape(p)))
        if (e$single.block)
          draw <- add(draw, quote(kappa <- alpha * e[["y"]]))
        else
          draw <- add(draw, quote(kappa <- alpha * e[["y"]] * exp(-p[["e_"]])))
      }
    } else if (modus == "vargamma") {
      if (e$family$alpha.fixed) {
        alpha <- e$family$get_shape()
        if (e$single.V.block)
          kappa <- alpha * e$family[["sigmasq"]]
        else {
          kappa0 <- alpha * e$family[["sigmasq"]]
          draw <- add(draw, quote(kappa <- kappa0 * p[["Q_"]]))
        }
      } else {
        draw <- add(draw, quote(alpha <- e$family$get_shape(p)))
        if (e$single.V.block) {
          draw <- add(draw, quote(kappa <- alpha * e$family[["sigmasq"]]))
        } else {
          draw <- add(draw, quote(kappa <- alpha * e$family[["sigmasq"]] * p[["Q_"]]))
        }
      }
    }
  }
  if (e$modeled.Q && modus == "regular") {
    if (e$Q0.type == "symm")
      draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
    else
      draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
  }

  if (usePX)
    draw <- add(draw, bquote(inv_xi <- 1 / p[[.(name_xi)]]))
  else
    draw <- add(draw, quote(inv_xi <- 1))
  draw <- add(draw, bquote(coef_raw <- inv_xi * p[[.(name)]]))
  if (q0 == 1L)
    draw <- add(draw, quote(M_coef_raw <- coef_raw))
  else
    draw <- add(draw, bquote(M_coef_raw <- matrix(coef_raw, nrow=.(l), byrow=TRUE)))

  draw <- add(draw, bquote(tau <- .(if (e$sigma.fixed) 1 else quote(1 / p[["sigma_"]]^2))))

  if (usePX) {
    if (PX$vector) {
      if (PX$sparse) {  # sparse M_ind, Dv
        M_ind <- kronecker(rep.int(1, l), CdiagU(q0))  # dgC
        mat_sum_xi <- make_mat_sum(M0=PX_Q0, M1=crossprod_sym(M_ind, XX))
        chol_xi <- build_chol(mat_sum_xi(crossprod_sym(M_ind, XX)))
        draw <- draw |>
          add(quote(Dv <- M_ind)) |>
          add(quote(attr(Dv, "x") <- as.numeric(M_coef_raw)))
        if (PX$data.scale)
          draw <- add(draw, quote(chol_xi$update(mat_sum_xi(crossprod_sym(Dv, XX)))))
        else
          draw <- add(draw, quote(chol_xi$update(mat_sum_xi(crossprod_sym(Dv, XX), w1=tau))))
      } else {  # matrix M_ind, Dv
        M_ind <- kronecker(rep.int(1, l), diag(q0))
        chol_xi <- build_chol(crossprod_sym(M_ind, XX) + PX_Q0)
        draw <- add(draw, quote(Dv <- coef_raw * M_ind))
        if (PX$data.scale)
          draw <- add(draw, quote(chol_xi$update(crossprod_sym(Dv, XX) + PX_Q0)))
        else
          draw <- add(draw, quote(chol_xi$update(crossprod_sym(Dv, XX) * tau + PX_Q0)))
      }
      if (PX$data.scale)
        draw <- add(draw, bquote(xi <- drawMVN_cholQ(chol_xi, crossprod_mv(Dv, Xy) + PX_Qmu0, sd=.(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])))))
      else
        draw <- add(draw, quote(xi <- drawMVN_cholQ(chol_xi, crossprod_mv(Dv, Xy) * tau + PX_Qmu0)))
    } else {  # scalar xi
      if (modus == "regular") {
        if (PX$data.scale) {
          draw <- draw |>
            add(quote(V <- 1 / (dotprodC(coef_raw, XX %m*v% coef_raw) + PX_Q0))) |>
            add(quote(E <- V * (dotprodC(coef_raw, Xy) + PX_Qmu0))) |>
            add(bquote(xi <- rnorm(1L, mean=E, sd=.(if (e$sigma.fixed) quote(sqrt(V)) else quote(p[["sigma_"]] * sqrt(V))))))
        } else {
          draw <- draw |>
            add(quote(V <- 1 / (dotprodC(coef_raw, XX %m*v% coef_raw) * tau + PX_Q0))) |>
            add(quote(E <- V * (dotprodC(coef_raw, Xy) * tau + PX_Qmu0))) |>
            add(quote(xi <- rnorm(1L, mean=E, sd=sqrt(V))))
        }
      } else {
        if (modus == "var" || modus == "vargamma")
          draw <- add(draw, bquote(Hz.xi <- dotprodC(coef_raw, crossprod_mv(X, rMLiG(.(e[["n"]]), 0.5, vkappa)))))
        if (modus == "gamma")
          draw <- add(draw, bquote(Hz.xi <- dotprodC(coef_raw, crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa)))))
        else if (modus == "vargamma")
          draw <- add(draw, bquote(Hz.xi <- Hz.xi + dotprodC(coef_raw, crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa)))))
        log.kappa.xi <- log(PX$prior$a) + sqrt(PX$prior$precision / PX$prior$a) * PX$prior$mean
        draw <- add(draw, bquote(Hz.xi <- Hz.xi + sqrt(PX$prior$precision / PX$prior$a) * rMLiG(1L, PX$prior$a, log.kappa=log.kappa.xi)))
        if (modus == "vargamma")
          draw <- add(draw, quote(xi <- Hz.xi / (2 * dotprodC(coef_raw, XX %m*v% coef_raw) + PX$prior$precision / PX$prior$a)))
        else
          draw <- add(draw, quote(xi <- Hz.xi / (dotprodC(coef_raw, XX %m*v% coef_raw) + PX$prior$precision / PX$prior$a)))
      }
    }
    if (!is.null(S)) {
      stop("TBI: inequality constrained coefficients with parameter expansion")  # --> f.c. of xi also TMVN(?)
      draw <- add(draw, bquote(p[[.(name)]] <- xi * coef_raw))  # need updated coefficient as input for TMVN sampler
    }
  } else {
    xi <- 1  # no parameter expansion
  }
  if (modus != "regular") {
    # for use in computation of acceptance ratio
    draw <- add(draw, bquote(sigma2_raw <- (inv_xi * p[[.(name_sigma)]])^2))
  }

  if (e$modeled.Q && modus == "regular") rm(XX)  # XX recomputed at each iteration

  if (gl) {  # group-level predictors
    if (usePX) {
      # MH step
      if (PX$vector)
        draw <- add(draw, bquote(logr <- .(glp$p0) * sum(log(abs(xi/p[[.(name_xi)]])))))
      else
        draw <- add(draw, bquote(logr <- .(q0 * glp$p0) * log(abs(xi/p[[.(name_xi)]]))))
      draw <- draw |>
        add(bquote(delta <- (xi * inv_xi) * p[[.(name_gl)]])) |>
        # in case !sigma.fixed, we need sigma^-2 factor in 2nd term on next line(?)
        add(quote(logr <- logr - 0.5 * dotprodC(delta, glp[["Q0"]] %m*v% delta))) |>
        add(bquote(delta <- p[[.(name_gl)]])) |>
        # in case !sigma.fixed, we need sigma^-2 factor in 2nd term on next line(?)
        add(quote(logr <- logr + 0.5 * dotprodC(delta, glp[["Q0"]] %m*v% delta))) |>
        add(bquote(if (logr < log(runif(1L))) xi <- p[[.(name_xi)]]))  # NB reject, otherwise accept
    }
    # NB use old value of xi to get 'raw' group-level coefficients
    if (q0 == 1L) {
      draw <- add(draw, bquote(M_coef_raw <- M_coef_raw - glp[["X"]] %m*v% (inv_xi * p[[.(name_gl)]])))
    } else {
      draw <- add(draw, bquote(M_coef_raw <- M_coef_raw - glp[["X"]] %m*m% matrix(inv_xi * p[[.(name_gl)]], nrow=.(glp$p0), byrow=TRUE)))
    }
  }
  if (usePX) {
    draw <- draw |>
      add(bquote(p[[.(name_xi)]] <- xi)) |>
      add(quote(inv_xi <- 1 / xi))
  }

  if (!is.null(priorA)) {
    if (q0 == 1L)  # DAM is lD vector
      draw <- add(draw, quote(DAM <- DA %m*v% M_coef_raw))
    else  # DAM is of type matrix (lD x q0) as M_coef_raw is
      draw <- add(draw, quote(DAM <- DA %m*m% M_coef_raw))
    switch(var,
      unstructured = {
        draw <- add(draw, bquote(SSR <- .rowSums((DAM %m*m% p[[.(name_Qraw)]]) * DAM, .(lD), .(q0))))
      },
      diagonal = {
        #draw <- add(draw, bquote(SSR <- .rowSums((DAM %*% Cdiag(p[[.(name_temp)]]$Qraw)) * DAM, .(lD), .(q0))))
        draw <- add(draw, bquote(SSR <- .rowSums(rep_each(p[[.(name_Qraw)]], .(lD)) * DAM^2, .(lD), .(q0))))
      },
      scalar = {
        if (q0 == 1L) {
          if (is.null(Q0)) {
            draw <- add(draw, bquote(SSR <- p[[.(name_Qraw)]] * DAM^2))
          } else {
            draw <- add(draw, bquote(SSR <- Q0 * p[[.(name_Qraw)]] * DAM^2))
          }
        } else {
          if (is.null(Q0)) {
            draw <- add(draw, bquote(SSR <- p[[.(name_Qraw)]] * .rowSums(DAM * DAM, .(lD), .(q0))))
          } else {
            draw <- add(draw, bquote(SSR <- p[[.(name_Qraw)]] * .rowSums((DAM %m*m% Q0) * DAM, .(lD), .(q0))))
          }
        }
      }
    )
    switch(priorA$type,
      invchisq = {
        if (is.list(priorA$df)) {
          draw <- add(draw, bquote(p[[.(name_df)]] <- priorA$draw_df(p[[.(name_df)]], 1 / p[[.(name_omega)]])))
          if (is.list(priorA$scale))
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(p[[.(name_df)]], df.data.omega, SSR, 1 / p[[.(name_omega)]])))
          else
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(p[[.(name_df)]], df.data.omega, SSR)))
        } else {
          if (is.list(priorA$scale))
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(df.data.omega, SSR, 1 / p[[.(name_omega)]])))
          else
            draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(df.data.omega, SSR)))
        }
      },
      exp = {
        aparA <- 2 / priorA$scale
        draw <- add(draw, bquote(p[[.(name_omega)]] <- priorA$draw(1 - q0/2, aparA, SSR)))
      }
    )
    # then update QA
    draw <- add(draw, bquote(QA <- crossprod_sym(DA, 1 / p[[.(name_omega)]])))
  }

  if (strucA$update.Q) {
    draw <- draw |>
      add(quote(p <- strucA$draw(p, M_coef_raw))) |>
      add(bquote(QA <- strucA$update_Q(QA, p[[.(strucA$name_ext)]])))
  }

  # draw V
  if (modus != "regular") {
    if (usePX) {
      name_sigma_raw <- paste0(name, "_sigma_raw_")  # MH parameter
    }
    if (control$MHprop == "LNRW") {
      # log-normal proposal
      sd.sigma.prop <- 0.2  # start value of RW stdev
      if (usePX) {
        adapt <- function(ar) {
          if (ar[[name_sigma_raw]] < .2)
            sd.sigma.prop <<- sd.sigma.prop * runif(1L, 0.6, 0.9)
          else if (ar[[name_sigma_raw]] > .7)
            sd.sigma.prop <<- sd.sigma.prop * runif(1L, 1.1, 1.5)
        }
      } else {
        adapt <- function(ar) {
          if (ar[[name_sigma]] < .2)
            sd.sigma.prop <<- sd.sigma.prop * runif(1L, 0.6, 0.9)
          else if (ar[[name_sigma]] > .7)
            sd.sigma.prop <<- sd.sigma.prop * runif(1L, 1.1, 1.5)
        }
      }
      draw <- add(draw, quote(sigma2_raw.star <- sigma2_raw * exp(rnorm(1L, sd=sd.sigma.prop))))
    } else {  # GiG proposal
      if (prior$type == "invchisq")
        draw <- add(draw, quote(sigma2_raw.star <- 1/rchisq_scaled(1L, q + prior$df, psi = dotprodC(coef_raw, coef_raw) + prior$df * prior$scale)))
      else  # exponential prior
        draw <- add(draw, quote(sigma2_raw.star <- Crgig(1L, p = 1 - 0.5*q, a=2/prior$scale, b=dotprodC(coef_raw, coef_raw))))
    }
    switch(prior$type,
      invchisq = {
        draw <- add(draw, bquote(
          log.ar <- 
            .(if (control$MHprop == "LNRW") 0.5 * (q + prior$df) else 0.5 * (q + prior$df + 2)) * log(sigma2_raw/sigma2_raw.star) +
            0.5 * prior$df * prior$scale * (1/sigma2_raw - 1/sigma2_raw.star)
        ))
      },
      exp = {
        draw <- add(draw, bquote(
          log.ar <- 
            .(if (control$MHprop == "LNRW") 0.5*(q - 2) else 0.5*q) * log(sigma2_raw/sigma2_raw.star) +
            (sigma2_raw - sigma2_raw.star) / prior$scale
        ))
      },
      stop("unsupported prior for gen() variance component and gamma sampling distribution")
    )
    draw <- add(draw, quote(
      log.ar <- log.ar +
        sqrt(a) * sum(coef_raw) * (1/sqrt(sigma2_raw) - 1/sqrt(sigma2_raw.star)) +
        a * sum(exp(-coef_raw/(sqrt(a * sigma2_raw))) - exp(-coef_raw/(sqrt(a * sigma2_raw.star))))
    ))
    if (usePX) {
      # store raw sigma, as it gets accepted/rejected, where sigma always changes due to PX
      draw <- draw |>
        add(bquote(p[[.(name_sigma_raw)]] <- sqrt(if (log(runif(1L)) < log.ar) sigma2_raw.star else sigma2_raw))) |>
        add(bquote(p[[.(name_sigma)]] <- abs(xi) * p[[.(name_sigma_raw)]]))
    } else {
      draw <- add(draw, bquote(p[[.(name_sigma)]] <- abs(xi) * sqrt(if (log(runif(1L)) < log.ar) sigma2_raw.star else sigma2_raw)))
    }
  } else if (var == "unstructured") {
    if (is.list(prior$scale)) {
      # Qraw stored in p; NB Qv has changed because of updated xi (PX), but Qraw has not
      if (prior$scale$common) {
        draw <- draw |>
          add(bquote(psi0 <- psi0 + sum(diagC(p[[.(name_Qraw)]])))) |>
          add(quote(psi0 <- rchisq_scaled(1L, q0 * prior$df + prior$scale$df, psi = psi0)))
      } else {
        draw <- draw |>
          add(bquote(psi0 <- psi0 + diagC(p[[.(name_Qraw)]]))) |>
          add(bquote(psi0 <- rchisq_scaled(.(q0), prior$df + prior$scale$df, psi=psi0)))
      }
      draw <- add(draw, quote(Qraw <- rWishart(1L, df=df, Sigma = inverseSPD(add_diagC(crossprod_sym(M_coef_raw, QA), psi0)))[,,1L]))
    } else {
      # TODO draw V_raw using invWishart and form Qv by solving --> 1 instead of 2 solves
      draw <- add(draw, quote(Qraw <- rWishart(1L, df=df, Sigma = inverseSPD(psi0 + crossprod_sym(M_coef_raw, QA)))[,,1L]))
    }
    if (usePX && PX$vector)
      draw <- add(draw, quote(inv_xi_qform <- base_tcrossprod(inv_xi)))
    else
      draw <- add(draw, quote(inv_xi_qform <- inv_xi^2))
    draw <- draw |>
      add(quote(Qv <- Qraw * inv_xi_qform)) |>
      # convert precision matrix to standard errors and correlations
      add(quote(p[names_se_cor] <- prec2se_cor(Qv)))
  } else {
    if (var == "diagonal") {
      draw <- add(draw, bquote(SSR <- .colSums(M_coef_raw * (QA %m*m% M_coef_raw), .(l), .(q0))))
    } else {  # scalar variance
      if (q0 == 1L) {
        if (is.null(Q0))
          draw <- add(draw, quote(SSR <- dotprodC(M_coef_raw, QA %m*v% M_coef_raw)))
        else
          draw <- add(draw, quote(SSR <- Q0 * dotprodC(M_coef_raw, QA %m*v% M_coef_raw)))
      } else {
        if (is.null(Q0))
          draw <- add(draw, quote(SSR <- sum(M_coef_raw * (QA %m*m% M_coef_raw))))
        else
          draw <- add(draw, quote(SSR <- sum(M_coef_raw * (QA %m*m% (M_coef_raw %m*m% Q0)))))
      }
    }
    switch(prior$type,
      invchisq = {
        if (is.list(prior$scale))
          draw <- add(draw, bquote(Qraw <- 1 / prior$draw(df.data, SSR, p[[.(name_Qraw)]])))
        else
          draw <- add(draw, quote(Qraw <- 1 / prior$draw(df.data, SSR)))
      },
      exp = {
        apar <- 2 / prior$scale
        draw <- add(draw, quote(Qraw <- 1 / prior$draw(ppar, apar, SSR)))
      }
    )
    # NB if xi is a vector so is Qv, even for scalar Qraw
    draw <- add(draw, quote(Qv <- Qraw * inv_xi^2))
    if (is.null(Q0)) {
      draw <- add(draw, bquote(p[[.(name_sigma)]] <- sqrt(1/Qv)))
    } else {
      draw <- draw |>
        add(quote(sqrtQv <- sqrt(Qv))) |>
        add(quote(Qv <- scale_mat(Q0, sqrtQv))) |>
        add(bquote(p[[.(name_sigma)]] <- 1/sqrtQv))
    }
  }
  if (!is.null(priorA) || is.list(prior$scale) || strucA$update.Q) {
    draw <- add(draw, bquote(p[[.(name_Qraw)]] <- Qraw))
  }
  if (!is.null(e$control$CG) || e$control$expanded.cMVN.sampler ||
      (!is.null(AR1.inferred) && AR1sampler$MH$type != "TN")) {
    draw <- add(draw, bquote(p[[.(name_Qv)]] <- Qv))
  }

  # draw coefficients
  if (gl) {
    i.v <- seq_len(q)  # indices for random effect vector in u=(v, alpha)
    i.alpha <- (q + 1L):(q + glp$p0 * q0)  # indices for group-level effect vector in u=(v, alpha)
  }
  if (in_block) {
    if (gl) {
      draw <- add(draw, bquote(p[[.(name_Q)]] <- kron_prod(glp[["QA.ext"]], Qv, values.only=TRUE)))
    } else {
      if (modus == "regular") {
        draw <- add(draw, bquote(p[[.(name_Q)]] <- kron_prod(QA, Qv, values.only=TRUE)))
      } else {
        draw <- add(draw, bquote(p[[.(name_Q)]] <- rep.int(1 / (a * p[[.(name_sigma)]]^2), .(q))))
      }
    }
    draw <- add(draw, bquote(p[[.(name)]] <- xi * coef_raw))
  } else {
    if (!is.null(AR1.inferred)) {
      # in non-blocked case p[[name_Q]] and kron_prod used only for drawing AR1 parameter
      draw <- add(draw, bquote(p[[.(name_Q)]] <- kron_prod(QA, Qv, values.only=TRUE)))
    }
    if (gl) {  # block sampling of (coef, glp)
      if (e$modeled.Q)
        draw <- add(draw, quote(attr(glp[["XX.ext"]], "x") <- c(XX@x, glp[["Q0"]]@x)))
      if (strucA$update.Q)
        draw <- add(draw, quote(glp[["QA.ext"]] <- crossprod_sym(glp[["IU0"]], QA)))
      draw <- draw |>
        add(quote(Qlist <- update(glp[["XX.ext"]], glp[["QA.ext"]], Qv, 1/tau))) |>
        add(bquote(coef <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=Qlist$Q, Imult=Qlist$Imult, Xy=c(Xy, glp[["Q0b0"]]))[[.(name)]])) |>
        add(bquote(p[[.(name)]] <- coef[i.v])) |>
        add(bquote(p[[.(name_gl)]] <- coef[i.alpha]))
    } else {  # no blocking, no group-level component
      if (modus == "regular") {
        draw <- draw |>
          add(quote(Qlist <- update(XX, QA, Qv, 1/tau))) |>
          add(bquote(p[[.(name)]] <- MVNsampler$draw(p, .(if (e$sigma.fixed) 1 else quote(p[["sigma_"]])), Q=Qlist$Q, Imult=Qlist$Imult, Xy=Xy)[[.(name)]]))
      } else {
        if (modus == "var" || modus == "vargamma")
          draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(e[["n"]]), 0.5, vkappa))))
        if (modus == "gamma")
          draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa))))
        else if (modus == "vargamma")
          draw <- add(draw, bquote(Hz <- Hz + crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa))))
        draw <- add(draw, bquote(Hz <- Hz + rMLiG(.(q), a, a) / (p[[.(name_sigma)]] * sqrt(a))))
        if (modus == "vargamma")
          draw <- add(draw, bquote(cholHH$update(2 * XX, mult=1 / (a * p[[.(name_sigma)]]^2))))
        else
          draw <- add(draw, bquote(cholHH$update(XX, mult=1 / (a * p[[.(name_sigma)]]^2))))
        draw <- add(draw, bquote(p[[.(name)]] <- cholHH$solve(Hz)))
      }
    }
  }

  if (modus == "var" || modus == "vargamma") {
    if (e$single.V.block) {
      if (is_ind_matrix(X) && q < e[["n"]])
        draw <- add(draw, bquote(p[["Q_"]] <- X %m*v% exp(-p[[.(name)]])))
      else
        draw <- add(draw, bquote(p[["Q_"]] <- exp(X %m*v% (-p[[.(name)]]))))
    } else {
      if (is_ind_matrix(X) && q < e[["n"]])
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * (X %m*v% exp(-p[[.(name)]]))))
      else
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(X %m*v% (-p[[.(name)]]))))
    }
  } else {
    if (e$e.is.res)
      draw <- add(draw, bquote(mv_update(p[["e_"]], plus=FALSE, X, p[[.(name)]])))
    else if (e$single.block && length(e$mod) == 1L)
      draw <- add(draw, bquote(p[["e_"]] <- X %m*v% p[[.(name)]]))
    else
      draw <- add(draw, bquote(mv_update(p[["e_"]], plus=TRUE, X, p[[.(name)]])))
  }
  draw <- add(draw, quote(p))
  # END draw function


  start <- function(p) {}

  # TODO check formats of user-provided start values
  if (!in_block) {
    if (gl) {
      # TODO conditional sampling if one of p[[.(name)]] and p[[.(name_gl)]] is provided
      start <- start |>
        add(bquote(coef <- MVNsampler$start(p, e$scale_sigma)[[.(name)]])) |>
        #if (!is.null(glp$R))
        #  start <- add(start, quote(coef <- constrain_cholQ(coef, ops$ch, glp$R)))
        add(bquote(if (is.null(p[[.(name)]])) p[[.(name)]] <- coef[i.v])) |>
        add(bquote(if (is.null(p[[.(name_gl)]])) p[[.(name_gl)]] <- coef[i.alpha]))
    } else {
      if (modus == "regular") {
        start <- add(start, bquote(if (is.null(p[[.(name)]])) p[[.(name)]] <- MVNsampler$start(p, e$scale_sigma)[[.(name)]]))
      } else {
        start <- add(start, bquote(if (is.null(p[[.(name)]])) p[[.(name)]] <- Crnorm(.(q))))
      }
    }
  }

  if (usePX) {
    if (PX$vector)
      start <- add(start, bquote(if (is.null(p[[.(name_xi)]])) p[[.(name_xi)]] <- rep.int(1, .(q0))))
    else
      start <- add(start, bquote(if (is.null(p[[.(name_xi)]])) p[[.(name_xi)]] <- 1))
  }

  if (!is.null(priorA) || is.list(prior$scale) || strucA$update.Q) {
    switch(var,
      unstructured = {
        start <- add(start, bquote(if (is.null(p[[.(name_Qraw)]])) p[[.(name_Qraw)]] <- rWishart(1L, prior$df, if (is.list(prior$scale)) (1/psi0)*diag(q0) else inverseSPD(psi0))[,,1L]))
      },
      diagonal=, scalar = {
        if (prior$type == "exp")
          start <- add(start, bquote(if (is.null(p[[.(name_Qraw)]])) p[[.(name_Qraw)]] <- rexp(.(if (var == "diagonal") q0 else 1L), rate=1/prior$scale)))
        else
          start <- add(start, bquote(if (is.null(p[[.(name_Qraw)]])) p[[.(name_Qraw)]] <- rchisq_scaled(.(if (var == "diagonal") q0 else 1L), prior$df, psi=prior$psi0)))
      }
    )
    if (strucA$update.Q) {
      start <- add(start, quote(p <- strucA$start(p)))
    }
    if (!is.null(priorA)) {
      if (is.list(priorA$df))
        start <- add(start, bquote(if (is.null(p[[.(name_df)]])) p[[.(name_df)]] <- runif(1L, 1, 25)))
      if (is.list(priorA$df) || is.list(priorA$scale) || !is.null(AR1.inferred))
        start <- add(start, bquote(if (is.null(p[[.(name_omega)]])) p[[.(name_omega)]] <- runif(.(lD), 0.75, 1.25)))
    }
  } else if (modus != "regular") {
    start <- add(start, bquote(if (is.null(p[[.(name_sigma)]])) p[[.(name_sigma)]] <- rexp(1L)))
  }

  if (!is.null(AR1.inferred)) {
    start <- add(start, bquote(if (is.null(p[[.(name_AR1)]])) p[[.(name_AR1)]] <- runif(1L, 0.25, 0.75)))
    if (AR1sampler$MH$type == "TN")
      start <- add(start, bquote(if (is.null(p[[.(name_Q)]])) p[[.(name_Q)]] <- kron_prod(QA.template$update(p[[.(name_AR1)]]), generate_Qv(), values.only=TRUE)))
    else {
      start <- add(start, bquote(if (is.null(p[[.(name_Qv)]])) p[[.(name_Qv)]] <- generate_Qv()))
      if (is.null(priorA)) {
        start <- add(start, bquote(if (is.null(p[[.(name_Q)]])) p[[.(name_Q)]] <- kron_prod(QA.template$update(p[[.(name_AR1)]]), generate_Qv(), values.only=TRUE)))
      } else
        start <- add(start, bquote(if (is.null(p[[.(name_Q)]])) p[[.(name_Q)]] <- kron_prod(crossprod_sym(DA.template$update(p[[.(name_AR1)]]), p[[.(name_omega)]]), p[[.(name_Qv)]], values.only=TRUE)))
    }
  }

  start <- add(start, quote(p))

  if (in_block && (!is.null(e$control$CG) || e$control$expanded.cMVN.sampler)) {
    if (q0 == 1L) {
      drawMVNvarQ <- function(p) crossprod_mv(DA, sqrt(p[[name_Qv]]) * Crnorm(lD))
    } else {
      if (var == "unstructured")
        cholQV <- build_chol(rWishart(1L, q0, diag(q0))[,,1L])
      else
        cholQV <- build_chol(runif(q0, 0.5, 1.5))
      drawMVNvarQ <- function(p) {
        y2 <- Crnorm(q0 * lD)
        dim(y2) <- c(q0, lD)
        cholQV$update(p[[name_Qv]])
        cholQV$Ltimes(y2, transpose=FALSE) %m*m% DA
      }
    }

  }

  self
}


setup_priorGMRFsampler <- function(mc, Qv) {
  # TODO
  # - support local scale parameters DA --> diag(omega_i)^-1/2 DA
  # - in some cases perm=TRUE may be more efficient
  mc$cholDD <- build_chol(tcrossprod(mc$DA), control=chol_control(perm=FALSE))
  mc$cholQv <- build_chol(Qv, control=chol_control(perm=FALSE))
  rGMRF <- function(Qv) {}
  rGMRF <- add(rGMRF, quote(cholQv$update(Qv)))
  if (mc[["q0"]] == 1L) {
    rGMRF <- rGMRF |>
      add(bquote(Z <- Crnorm(.(nrow(mc$DA))))) |>
      add(quote(Z <- cholQv$solve(Z, system="Lt"))) |>
      add(quote(crossprod_mv(DA, cholDD$solve(Z))))
  } else {
    rGMRF <- rGMRF |>
      add(bquote(Z <- matrix(Crnorm(.(mc[["q0"]] * nrow(mc$DA))), nrow=.(mc[["q0"]])))) |>
      add(quote(Z <- cholQv$solve(Z, system="Lt"))) |>
      add(quote(coef <- crossprod(DA, cholDD$solve(t.default(Z))))) |>
      add(quote(as.numeric(t.default(coef))))
  }
  mc$rGMRFprior <- rGMRF
  environment(mc$rGMRFprior) <- mc
}

setup_priorMVNsampler <- function(mc, Q) {
  rMVN <- function(p, Q) {}
  if (mc$is.proper) {
    mc$priorMVNsampler <- create_TMVN_sampler(Q, update.Q=TRUE, update.mu=FALSE, name="coef", R=mc$R)
  } else {
    if (is.null(mc$R)) stop("cannot sample from improper GMRF without constraints")
    # for IGMRF add tcrossprod(mc$R) to Q to make it non-singular (--> inefficient)
    # TODO maybe mc$R is removed and stored in (posterior) MVNsampler --> keep a reference(!) to R in mc
    mc$mat_sum_prior <- make_mat_sum(M0=economizeMatrix(tcrossprod(mc$R), symmetric=TRUE), M1=Q)
    mc$priorMVNsampler <- create_TMVN_sampler(mc$mat_sum_prior(Q), update.Q=TRUE, update.mu=FALSE, name="coef", R=mc$R)
    rMVN <- add(rMVN, quote(Q <- mat_sum_prior(Q)))
  }
  rMVN <- add(rMVN, quote(priorMVNsampler$draw(p, Q=Q)[["coef"]]))
  mc$rMVNprior <- rMVN
  environment(mc$rMVNprior) <- mc
}

#' Set computational options for the sampling algorithms used for a 'gen' model component
#'
#' @export
#' @param MHprop MH proposal for the variance component in case of a MLiG prior
#'  on the coefficients. The two options are "GiG" for a generalized inverse gamma
#'  proposal, and "LNRW" for a log-normal random walk proposal. The former should
#'  approximate the conditional posterior quite well provided MLiG parameter \code{a}
#'  is large, such that the coefficients' prior is approximately normal.
#' @returns A list of computational options regarding a 'gen' model component.
gen_control <- function(MHprop = c("GiG", "LNRW")) {
  list(MHprop = match.arg(MHprop))
}
