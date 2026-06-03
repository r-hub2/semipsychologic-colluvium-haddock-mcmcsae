#' Specify a Student-t sampling distribution
#'
#' This function can be used in the \code{family} argument of \code{\link{create_sampler}}
#' or \code{\link{generate_data}} to specify a Student-t sampling distribution.
#'
#' @export
#' @param link the name of a link function. Currently the only allowed link function
#'  for a Student-t distribution is \code{"identity"}.
#' @param df.prior prior specification for the (scalar) number of degrees of
#'  freedom of the Studen-t sampling distribution. Currently allowed priors are
#'  gamma priors specified using \code{\link{pr_gamma}}, or a degenerate prior
#'  fixing the degrees of freedom parameter at a fixed value.
#' @param control sampling options for the degrees of freedom parameter
#'  (if not fixed).
#' @param ... further arguments passed to \code{\link{f_gaussian}}.
#' @returns A family object.
f_student_t <- function(link="identity", df.prior=pr_gamma(2, 0.1),
                        control=student_t_control(), ...) {
  link <- match.arg(link)
  control <- check_student_t_control(control)
  if (df.prior[["type"]] == "fixed") {
    var.model <- ~ vfac(prior=pr_invchisq(df=df.prior[["value"]]))
  } else if (df.prior[["type"]] == "gamma") {
    var.model <- ~ vfac(prior=pr_invchisq(df=list(
      alpha=df.prior[["shape"]], beta0=df.prior[["rate"]],
      proposal=control[["proposal"]], tau=control[["MHscale"]], adapt=control[["adapt"]]
    )))
  } else {
    stop("unsupported prior for Student-t degrees of freedom parameter")
  }
  df.prior$init(n=1L)
  dotargs <- list(...)
  # TODO
  # - allow PX
  # - allow (fixed) vector df
  # - process MH sampling options through student_t_control
  lgau <- do.call(f_gaussian, dotargs)
  if (!is.null(lgau[["prec.mat"]])) {
    stop("argument 'prec.mat' for specifying a non-diagonal sampling ",
         "covariance matrix, is not supported for Student-t family")
  }
  # TODO also keep non-trivial lgau[["var.model"]] separately for llh computation
  if (!is.null(lgau[["var.model"]]))
    var.model <- as.formula(paste(deparse(var.model), deparse(lgau[["var.model"]]), sep=" + "))
  lgau <- lgau[-which(names(lgau) %in% c("family", "var.model"))]
  names(lgau)[names(lgau) == "control"] <- "gaussian.control"
  c(list(family="student_t", df.prior=df.prior, var.model=var.model), lgau)
}

ff_student_t <- function(link="identity", df.prior=pr_gamma(2, 0.1),
                         control=student_t_control(),
                         sm, data, y=NULL, sub=NULL, ...) {
  family <- "student_t"
  linkinv <- make.link(link)$linkinv
  prior.only <- is.null(y)
  dotargs <- list(...)
  y.family <- ff_gaussian(link=link, var.prior=dotargs[["var.prior"]],
    var.vec=dotargs[["var.vec"]], prec.mat=NULL, var.model=dotargs[["var.model"]],
    logJacobian=dotargs[["logJacobian"]], control=dotargs[["gaussian.control"]],
    sm=sm, data=data, y=y, sub=sub
  )
  sc <- sm[["control"]]
  store_default <- function() y.family$store_default()
  if (is.function(y.family[["adapt"]]))
    adapt <- function(ar) y.family$adapt(ar)
  if (is.function(y.family[["drawMVNvarQ"]])) {
    cholQ <- y.family[["cholQ"]]
    drawMVNvarQ <- y.family[["drawMVNvarQ"]]
  }
  llh <- y.family[["llh"]]
  if (FALSE) {
    # TODO true student-t llh (Rao-Blackwellised), new make_rpredictive
    if (df.prior[["type"]] == "fixed") {
      llh_0 <- - 0.5*log(pi*df.prior[["value"]]) +
        lgamma(0.5*(df.prior[["value"]] + 1)) - lgamma(0.5*df.prior[["value"]])
      # TODO SSR should here be computed without the local DA scale parameters
      llh <- function(p) {
        if (sigma.fixed)
          llh_0 - 0.5*(df.prior[["value"]] + 1) * log(1 + p[[SSR.name]]/df.prior[["value"]])
        else
          llh_0 - n * log(p[[sd.name]]) - 0.5*(df.prior[["value"]] + 1) * log(1 + p[[SSR.name]]/(df.prior[["value"]] * p[[sd.name]]^2))
      }
    } else {
      llh <- function(p) {
        out <- - 0.5*log(pi * p[["name_df"]]) +
          lgamma(0.5*(p[[name_df]] + 1)) - lgamma(0.5 * p[[name_df]])
        if (sigma.fixed)
          out - 0.5*(p[[name_df]] + 1) * log(1 + p[[SSR.name]]/p[[name_df]])
        else
          out - n * log(p[[sd.name]]) - 0.5*(p[[name_df]] + 1) * log(1 + p[[SSR.name]]/(p[["value"]] * p[[sd.name]]^2))
      }
    }
  }
  llh_i <- y.family[["llh_i"]]
  self <- environment()
  copy_refs(y.family, self,
    c("n", "y", "sigma.fixed", "modeled.Q", "Q0", "Q0.type",
      "e.is.res", "Q_e", "Vmod", "scale.e", "scale.sigma",
      "rprior", "draw", "start", "MHpars", "make_rpredictive",
      "sd.name", "SSR.name", "compute_SSR")
  )
  rm(data, sm)
  self
}

#' Set computational options for the sampling algorithms
#'
#' @export
#' @param proposal "RW" for random walk Metropolis-Hastings or "mala"
#'  for Metropolis-adjusted Langevin.
#' @param MHscale (starting) scale of Metropolis-Hastings update.
#' @param adapt whether to adapt the scale of the proposal distribution
#'  during burnin to achieve better acceptance rates.
#' @returns A list with computational options.
student_t_control <- function(proposal="RW", MHscale=1, adapt=TRUE) {
  list(proposal=proposal, MHscale=MHscale, adapt=adapt)
}

check_student_t_control <- function(control) {
  if (is.null(control)) return(student_t_control())
  if (!is.list(control)) stop("control options must be specified as a list, preferably using the appropriate control setter function")
  defaults <- student_t_control()
  w <- whichv(names(control) %in% names(defaults), FALSE)
  if (length(w)) stop("unrecognized control parameters ", paste0(names(control)[w], collapse=", "))
  control <- modifyList(defaults, control, keep.null=TRUE)
  control$proposal <- match.arg(control[["proposal"]], c("RW", "mala"))
  if (!is_numeric_scalar(control[["MHscale"]]) || control[["MHscale"]] <= 0)
    stop("'MHscale' must be a positive numeric value")
  if (!is_logical_scalar(control[["adapt"]])) stop("'adapt' must be TRUE or FALSE")
  control
}
