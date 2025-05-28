#' Functions for specifying the method and corresponding options for sampling
#' from a possibly truncated and degenerate multivariate normal distribution
#'
#' These functions are intended for use in the \code{method} argument of \code{\link{create_TMVN_sampler}}.
#'
#' @param use.cholV whether to use the Cholesky factor of the variance instead
#'  of precision matrix for sampling. If \code{NULL} the choice is made based on a
#'  simple heuristic.
#' @param diagnostic whether information about violations of inequalities, bounces
#'  off inequality walls (for 'HMC' and 'HMCZigZag' methods) or gradient events
#'  (for 'HMCZigZag') is printed to the screen.
#' @param slice if \code{TRUE}, a Gibbs within slice sampler is used.
#' @param eps small positive value to control numerical robustness of the algorithm.
#' @param Tsim the duration of a Hamiltonian Monte Carlo simulated particle trajectory.
#'  This can be specified as either a single positive numeric value for a fixed
#'  simulation time, or as a function that is applied in each MCMC iteration to
#'  generates a simulation time.
#' @param max.events maximum number of events (reflections off inequality walls
#'  and for method 'HMCZigZag' also gradient events). Default is unlimited.
#'  Specifying a finite number may speed up the sampling but may also result
#'  in a biased sampling algorithm.
#' @param scale vector of Laplace scale parameters for method 'HMCZigZag'. It must be
#'  a positive numeric vector of length equal to one or the number of variables.
#' @param prec.eq positive numeric vector of length 1 or the number of equality restrictions,
#'  to control the precision by which the equality restrictions are imposed; the larger
#'  \code{prec.eq} the more precisely they will be imposed.
#' @param adapt experimental feature: if \code{TRUE} the rate parameter will be adapted
#'  in an attempt to make the sampling algorithm more efficient.
#' @param sharpness for method 'softTMVN', the sharpness of the soft inequalities; the larger the better
#'  the approximation of exact inequalities. It must be a positive numeric vector of length
#'  one or the number of inequality restrictions.
#' @param useV for method 'softTMVN' whether to base computations on variance instead of precision
#'  matrices.
#' @param CG use a conjugate gradient iterative algorithm instead of Cholesky updates for sampling
#'  the model's coefficients. This must be a list with possible components \code{max.it},
#'  \code{stop.criterion}, \code{verbose}. See the help for function \code{\link{CG_control}},
#'  which can be used to specify these options. Currently the preconditioner and scale
#'  options cannot be set for this use case.
#' @param PG.approx see \code{\link{sampler_control}}.
#' @param PG.approx.m see \code{\link{sampler_control}}.
#' @returns A method object, for internal use only.
#' @name TMVN-methods
NULL


#' @export
#' @rdname TMVN-methods
m_direct <- function(use.cholV=NULL) {
  list(name="direct", use.cholV=use.cholV)
}

#' @export
#' @rdname TMVN-methods
m_Gibbs <- function(slice=FALSE, eps=sqrt(.Machine$double.eps), diagnostic=FALSE) {
  list(name="Gibbs", slice=slice, eps=eps, diagnostic=diagnostic)
}

#' @export
#' @rdname TMVN-methods
m_HMC <- function(Tsim=pi/2, max.events=.Machine$integer.max, diagnostic=FALSE) {
  if (is.function(Tsim)) {
    Tsim_value <- Tsim()
  } else {
    Tsim_value <- as.numeric(Tsim)[1L]
    Tsim <- function() Tsim_value
  }
  if (!is.numeric(Tsim_value) || length(Tsim_value) != 1L || Tsim_value <= 0) stop("'Tsim' should be a single positive numeric value, or a function to generate one")
  max.event <- as.integer(max.events)
  if (length(max.events) != 1L || max.events <= 0) stop("'max.events must be a positive integer'")
  list(name="HMC", Tsim=Tsim, max.events=max.events, diagnostic=diagnostic)
}

#' @export
#' @rdname TMVN-methods
m_HMCZigZag <- function(Tsim=1, scale=1, prec.eq=NULL, diagnostic=FALSE,
                        max.events=.Machine$integer.max, adapt=FALSE) {
  if (is.function(Tsim)) {
    Tsim_value <- Tsim()
  } else {
    Tsim_value <- as.numeric(Tsim)[1L]
    Tsim <- function() Tsim_value
  }
  if (!is.numeric(Tsim_value) || length(Tsim_value) != 1L || Tsim_value <= 0) stop("'Tsim' should be a single positive numeric value, or a function to generate one")
  if (any(scale <= 0)) stop("Laplace distribution scale parameters must be positive")
  max.event <- as.integer(max.events)
  if (length(max.events) != 1L || max.events <= 0) stop("'max.events must be a positive integer'")
  list(name="HMCZigZag", Tsim=Tsim, scale=scale, prec.eq=prec.eq, diagnostic=diagnostic,
       max.events=max.events, adapt=adapt)
}

#' @export
#' @rdname TMVN-methods
m_softTMVN <- function(sharpness=100, useV=FALSE, CG=NULL, PG.approx=TRUE, PG.approx.m=-2L) {
  if (!is.null(CG)) {
    if (useV) stop("conjugate gradients algorithm not supported in combination with 'useV=TRUE'")
    CG <- check_CG_control(CG)
  }
  list(name="softTMVN", sharpness=sharpness, useV=useV, CG=CG, PG.approx=PG.approx, PG.approx.m=PG.approx.m)
}
