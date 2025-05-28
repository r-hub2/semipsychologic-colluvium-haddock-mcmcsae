
#' Specify a negative binomial sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a negative binomial sampling distribution.
#'
#' The negative binomial distribution with shape r and probability p
#' has density
#' \deqn{p(y|r, p) = {\Gamma(y + r)\over y!\Gamma(r)}(1-p)^r p^y}
#' with mean \eqn{\mu = E(y|r,p) = {rp\over 1-p}} and variance
#' \eqn{V(y|r,p) = \mu(1 + \mu/r)}. The second term of the variance
#' can be interpreted as overdispersion with respect to a Poisson
#' distribution, which would correspond to the limit \eqn{r \rightarrow \infty}.
#' So the reciprocal shape \eqn{1/r} is an overdispersion parameter,
#' which typically is inferred. It is assigned a default prior, which
#' may be changed through argument \code{inv.shape.prior}.
#'
#' The only supported link function is \code{"log"}. Strictly speaking
#' the relation between mean \eqn{\mu} and linear predictor \eqn{\eta} is
#' \deqn{\log\mu = \log r + \log{p\over 1-p} = \log r + \eta}
#' This way the likelihood function has the same form as that of logistic
#' binomial regression, so that a Polya-Gamma data augmentation sampling
#' algorithm can be employed. Note that the fact that the linear predictor
#' \eqn{\eta} does not include \eqn{\log r} effectively changes the
#' interpretation of its intercept.
# NB negative binomial regression with unknown shape does NOT belong to the
#    class of generalized linear models
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for the negative binomial sampling distribution is \code{"log"}.
#' @param shape.vec optional formula specification of unequal shape values. The
#'  negative binomial (vector) shape parameter is then equal to this vector of
#'  shape values, multiplied by the scalar shape parameter, whose prior is
#'  specified through \code{inv.shape.prior}.
#' @param inv.shape.prior Prior on the (scalar) \emph{reciprocal} shape parameter,
#'  i.e. the overdispersion parameter. Supported prior distributions are
#'  \code{\link{pr_fixed}} with a default value of 1, \code{\link{pr_invchisq}} and
#'  \code{\link{pr_gig}}. The current default is \code{pr_invchisq(df=1, scale=1)}.
#' @param control a list with computational options. These options can
#'  be specified using function \code{\link{negbin_control}}.
#' @returns A family object.
#' @references
#'  N. Polson, J.G. Scott and J. Windle (2013).
#'    Bayesian Inference for Logistic Models Using Polya-Gamma Latent Variables.
#'    Journal of the American Statistical Association 108(504), 1339-1349.
#'
#'  M. Zhou and L. Carin (2015).
#'    Negative Binomial Process Count and Mixture Modeling.
#'    IEEE Transactions on Pattern Analysis and Machine Intelligence 37(2), 307-320.
f_negbinomial <- function(link="log", shape.vec = ~ 1, inv.shape.prior = pr_invchisq(df=1, scale=1),
                          control = negbin_control()) {
  family <- "negbinomial"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
  if (!inherits(shape.vec, "formula"))
    stop("'shape.vec' must be a formula")
  switch(inv.shape.prior[["type"]],
    fixed = {
      if (inv.shape.prior[["value"]] <= 0)
        stop("negative binomial inverse shape (dispersion) parameter must be positive")
    },
    invchisq = {
      if (is.list(inv.shape.prior[["df"]]))
        stop("modelled degrees of freedom parameter not supported for negative binomial (inverse) shape parameter")
    },
    gig = {
    },
    stop("unsupported prior")
  )
  inv.shape.prior$init(n=1L)
  shape.fixed <- inv.shape.prior[["type"]] == "fixed"
  shape.scalar <- intercept_only(shape.vec)
  shape0 <- NULL
  get_shape <- function() stop("please call method 'init' first")
  init <- function(data, y=NULL) {
    get_shape <<- make_get_shape(data)
    if (!is.null(y)) {
      if (!is.numeric(y)) stop("non-numeric target value not allowed in case of negative binomial sampling distribution")
      # NB the algorithm still runs with negative responses, but it probably makes not much sense
      if (any(y < 0)) warn("negative response value(s)")
      if (any(abs(round(y) - y) > sqrt(.Machine$double.eps))) warn("non-integral values modelled by negative binomial sampling distribution")
    }
    y
  }
  make_get_shape <- function(data) {
    if (shape.scalar) {
      shape0 <<- 1
    } else {
      shape0 <<- get_var_from_formula(shape.vec, data)
      if (all(length(shape0) != c(1L, n_row(data))))
        stop("wrong length for shape vector")
      if (anyNA(shape0)) stop("missings in shape vector not allowed")
      if (any(shape0 <= 0)) stop("shape vector must be positive")
    }
    if (shape.fixed) {
      shape0 <<- shape0 / inv.shape.prior[["value"]]
      function(p) shape0
    } else {
      if (shape.scalar)
        function(p) p[["negbin_shape_"]]
      else
        function(p) shape0 * p[["negbin_shape_"]]
    }
  }
  if (!shape.fixed) {
    control <- check_negbin_control(control)
    make_draw_shape <- function(y) {
      # draw latent L_i (i=1:n) from its CRT f.c.
      mCRT <- control[["CRT.approx.m"]]
      draw <- function(p) {L <- CrCRT(y, get_shape(p), mCRT)}
      switch(inv.shape.prior[["type"]],
        # TODO check shape0 appearance in modeled scale invchisq and gig cases
        invchisq = {
          inv.shape.prior$make_draw()
          if (is.list(inv.shape.prior[["scale"]])) {
            if (shape.scalar)
              draw <- add(draw, quote(p$negbin_shape_ <- 1 / inv.shape.prior$draw(2*sum(L), 2*shape0*sum(log1pexpC(p[["e_"]])), p[["negbin_shape_"]])))
            else
              draw <- add(draw, quote(p$negbin_shape_ <- 1 / inv.shape.prior$draw(2*sum(L), 2*sum(shape0*log1pexpC(p[["e_"]])), p[["negbin_shape_"]])))
          } else {
            if (shape.scalar)
              draw <- add(draw, quote(p$negbin_shape_ <- 1 / inv.shape.prior$draw(2*sum(L), 2*shape0*sum(log1pexpC(p[["e_"]])))))
            else
              draw <- add(draw, quote(p$negbin_shape_ <- 1 / inv.shape.prior$draw(2*sum(L), 2*sum(shape0*log1pexpC(p[["e_"]])))))
          }
        },
        gig = {
          prior.p <- inv.shape.prior[["p"]]
          prior.a <- inv.shape.prior[["a"]]
          prior.b <- inv.shape.prior[["b"]]
          rgig <- GIGrvg::rgig
          if (shape.scalar)
            draw <- add(draw, quote(p$negbin_shape_ <- 1 / rgig(1L, prior.p - sum(L), prior.b + 2*shape0*sum(log1pexpC(p[["e_"]])), prior.a)))
          else
            draw <- add(draw, quote(p$negbin_shape_ <- 1 / rgig(1L, prior.p - sum(L), prior.b + 2*sum(shape0*log1pexpC(p[["e_"]])), prior.a)))
        }
      )
      draw
    }
  }
  make_llh <- function(y) {
    if (shape.fixed)
      llh_0 <- sum(negbinomial_coef(shape0, y))
    else
      llh_0 <- -sum(lgamma(y + 1))
    function(p) {
      neg_fitted <- -p[["e_"]]
      if (shape.fixed) {
        llh_0 + sum(shape0 * neg_fitted - (y + shape0) * log1pexpC(neg_fitted))
      } else {
        r <- get_shape(p)
        ny <- y + r
        llh_0 + sum(lgamma(ny) - lgamma(r)) + sum(r * neg_fitted - ny * log1pexpC(neg_fitted))  
      }
    }
  }
  make_llh_i <- function(y) {
    function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      if (shape.fixed) {
        if (shape.scalar)
          rep_each(negbinomial_coef(shape0, y[i]), nr) - shape0 * e_i - rep_each(y[i] + shape0, nr) * log1pexpC(-e_i)
        else
          rep_each(negbinomial_coef(shape0[i], y[i]), nr) - rep_each(shape0[i], nr) * e_i - rep_each(y[i] + shape0[i], nr) * log1pexpC(-e_i)
      } else {
        r <- as.numeric(as.matrix.dc(draws[["negbin_shape_"]], colnames=FALSE))
        if (shape.scalar)
          r <- shape0 * r
        else
          r <- r * rep_each(shape0[i], nr)
        yi <- rep_each(y[i], nr)
        nyi <- yi + r
        negbinomial_coef(r, yi) - r * e_i - nyi * log1pexpC(-e_i)
      }
    }
  }
  make_rpredictive <- function(newdata, weights=NULL) {
    # NB definition of rnbinom has p <-> 1-p
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata
    } else {
      nn <- nrow(newdata)
      get_shape <- make_get_shape(newdata)
    }
    if (is.null(weights)) {
      function(p, lp) {
        size <- get_shape(p)
        rnbinom(nn, size, mu=size*exp(lp))
      }
      #rnbinom(length(lp), size=get_shape(p), prob=1/(1 + exp(lp)))
    } else {
      function(p, lp) {
        size <- weights * get_shape(p)
        rnbinom(nn, size, mu=size*exp(lp))
      }
    }
  }
  environment()
}

#' Set computational options for the sampling algorithms
#'
#' @export
#' @param CRT.approx.m scalar integer specifying the degree of approximation to sampling
#'  from a Chinese Restaurant Table distribution. The approximation is based on Le Cam's theorem.
#'  Larger values yield a slower but more accurate sampler.
#' @returns A list with computational options for the sampling algorithm.
negbin_control <- function(CRT.approx.m=20L) {
  list(CRT.approx.m=CRT.approx.m)
}

check_negbin_control <- function(control) {
  if (is.null(control)) control <- list()
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- negbin_control()
  w <- which(!(names(control) %in% names(defaults)))
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$CRT.approx.m <- as.integer(control$CRT.approx.m)
  if (length(control$CRT.approx.m) != 1L || control$CRT.approx.m < 1L)
    stop("'CRT.approx.m' must be a positive scalar integer")
  control
}
