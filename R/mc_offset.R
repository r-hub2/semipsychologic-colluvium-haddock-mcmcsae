
#' Create a model component object for an offset, i.e. fixed, non-parametrised term in the linear predictor
#'
#' This function is intended to be used on the right hand side of the \code{formula} argument to
#' \code{\link{create_sampler}} or \code{\link{generate_data}}.
#'
#' @export
#' @param formula model formula.
#' @param value alternative specification of an offset as a single scalar
#'  value that is the same for each data unit.
#' @param name the name of the model component. This name is used in the output of the MCMC simulation
#'  function \code{\link{MCMCsim}}. By default the name will be 'mc_offset' with the number of the model term attached.
#' @returns An model component object with data and methods needed for dealing with an offset term
#'  in model estimation, and prior and posterior prediction. Intended for internal use by
#'  other package functions.
mc_offset <- function(formula, value=NULL, name="") {
  e <- sys.frame(-2L)
  type <- "mc_offset"
  if (name == "") stop("missing model component name")
  n <- e[["n"]]

  if (is.null(value)) {
    offset <- get_var_from_formula(formula, e[["data"]])
  } else {
    offset <- as.numeric(value)
    if (length(offset) != 1L || anyNA(offset)) stop("'value' argument of mc_offset must be a single numeric value")
  }
  if (all(length(offset) != c(1L, n))) stop("offset of unexpected length")

  # offset used for prediction, same as offset except for a Poisson model,
  #   where an internal offset is added to the approximating negbinomial model
  pred.offset <- offset
  internal.offset <- FALSE
  add_internal_offset <- function(value) {
    offset <<- offset + value
    internal.offset <<- TRUE
    # update lp, lp_update methods
    if (length(offset) == 1L)
      lp <<- function(p) copy_vector(rep.int(offset, n))
    else
      lp <<- function(p) copy_vector(offset)
    if (allv(offset, 0))
      lp_update <<- function(x, plus=TRUE, p) NULL
    else
      lp_update <<- function(x, plus=TRUE, p) v_update(x, plus, offset)
  }

  make_predict <- function(newdata=NULL, Xnew=NULL, verbose=TRUE) {
    if (is.null(newdata)) {
      if (is.null(Xnew)) {
        # in-sample prediction
        newoffset <- pred.offset
        nnew <- n
      } else {
        if (!is.vector(Xnew) || ncol(Xnew) != 1L) stop("wrong input for 'Xnew'")
        newoffset <- as.vector(Xnew)
        nnew <- length(newoffset)
      }
    } else {
      if (is.null(value))
        newoffset <- get_var_from_formula(formula, newdata)
      else
        newoffset <- pred.offset
      nnew <- nrow(newdata)
    }
    rm(newdata, Xnew, verbose)
    if (length(newoffset) == 1L && nnew != 1L)
      linpred <- function(p) copy_vector(rep.int(newoffset, nnew))
    else
      linpred <- function(p) copy_vector(newoffset)
    if (allv(newoffset, 0))
      linpred_update <- function(x, plus=TRUE, p) NULL
    else
      linpred_update <- function(x, plus=TRUE, p) v_update(x, plus, newoffset)
    environment()
  }
  # assume that lp, lp_update are used for model fitting only (--> include internal offset)
  if (length(offset) == 1L)
    lp <- function(p) copy_vector(rep.int(offset, n))
  else
    lp <- function(p) copy_vector(offset)
  if (allv(offset, 0))
    lp_update <- function(x, plus=TRUE, p) NULL
  else
    lp_update <- function(x, plus=TRUE, p) v_update(x, plus, offset)
  scalar <- length(pred.offset) == 1L
  # draws_linpred method acts on (subset of) mcdraws object, used in fitted() and pointwise log-likelihood llh_i functions
  # here internal offset (Poisson approx) is excluded
  draws_linpred <- function(obj, units=NULL, chains=NULL, draws=NULL, matrix=FALSE) {
    if (matrix) {
      nr <- (if (is.null(chains)) n_chains(obj) else length(chains)) * (if (is.null(draws)) n_draws(obj) else length(draws))
      if (is.null(units))
        out <- matrix(if (scalar) pred.offset else rep_each(pred.offset, nr), nr, n)
      else
        out <- matrix(if (scalar) pred.offset else rep_each(pred.offset[units], nr), nr, length(units))
    } else {
      nr <- if (is.null(draws)) n_draws(obj) else length(draws)
      out <- list()
      if (is.null(units))
        for (ch in seq_along(chains))
          out[[ch]] <- matrix(if (scalar) pred.offset else rep_each(pred.offset, nr), nr, n)
      else
        for (ch in seq_along(chains))
          out[[ch]] <- matrix(if (scalar) pred.offset else rep_each(pred.offset[units], nr), nr, length(units))
    }
    out
  }

  environment()
}
