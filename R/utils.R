
# faster version of rep(, each)
rep_each <- function(x, each) rep.int(x, rep.int(each, length(x)))

# determine whether a formula consists of an intercept only
intercept_only <- function(formula) {
  identical(as.character(update.formula(formula, ~ .)), as.character(~ 1))
  #isTRUE(all.equal(formula, ~ 1))
}

# if data is a (scalar) integer it is interpreted as sample size
n_row <- function(data) if (is.integer(data)) data else nrow(data)

# set up binary file for writing and write header
write_header <- function(con, n.iter, n.chain, n.par, parnames=NULL, single.prec=FALSE) {
  writeBin(as.integer(c(n.iter, n.chain, n.par)), con)
  if (is.null(parnames)) {
    writeBin(FALSE, con)
  } else {
    parnames <- as.character(parnames)
    if (length(parnames) != n.par) stop("wrong length for vector of parameter names")
    writeBin(TRUE, con)
    writeBin(parnames, con)
  }
  if (single.prec) {
    writeBin(4L, con)
  } else {
    writeBin(8L, con)
  }
}

# set up a binary file for reading and read header
# returns a vector with numbers of draws, chains and parameters 
read_header <- function(con) {
  out <- list(
    n.iter = readBin(con, "integer", n=1L),
    n.chain = readBin(con, "integer", n=1L),
    n.par = readBin(con, "integer", n=1L)
  )
  if (readBin(con, "logical", n=1L))
    out$parnames <- readBin(con, "character", n=out$n.par)
  precision <- readBin(con, "integer", n=1L)
  if (all(precision != c(4L, 8L))) stop("unexpected file content")
  out$single <- precision == 4L
  out
}

#' Generate artificial data according to an additive spatio-temporal model
#'
#' This function is used to generate data for several examples.
#'
#' @examples
#' \donttest{
#' ex <- mcmcsae_example()
#' str(ex)
#' }
#'
#' @export
#' @param n the size of the generated dataset.
#' @param family sampling distribution family, see \code{\link{create_sampler}}.
#' @returns A \code{list} containing the generated dataset, the values of the model
#'   parameters, and the model specification as a formula.
mcmcsae_example <- function(n=100L, family="gaussian") {
  if (n < 3L) stop("choose a size n >= 3")
  nA <- as.integer(sqrt(n))
  nT <- max(3L, as.integer(sqrt(n)))
  dat <- data.frame(
    x=rnorm(n),
    fA=qF(sample.int(nA, n, replace=TRUE)),
    fT=qF(sample.int(nT, n, replace=TRUE))
  )
  # fake data simulation
  model <- ~ reg(~x, prior=pr_normal(precision=1), name="beta") +
             gen(formula=~x, factor=~fA, name="v") +
             gen(factor=~RW2(fT), name="u")
  gd <- generate_data(model, data=dat, family=family)
  dat$y <- gd[["y"]]
  list(dat=dat, pars=gd[["pars"]], model=update.formula(model, y ~ .), family=family)
}

# extend a function by appending a line to its body
# expr must be an unevaluated (possibly multi-line) expression (use (b)quote or substitute)
add <- function(f, expr) {
  body(f)[[length(body(f)) + 1L]] <- expr
  f
}

#' Alphabetically order labels that may be composed of multiple factor labels separated by ':'
#'
#' @noRd
#' @param names character vector of labels of the categories of an interaction term.
#' @returns Names with each component alphabetically ordered by factor label separated by ':'.
order_interactions <- function(names)
  s_apply(strsplit(names, ":", fixed=TRUE), \(s) paste0(sort.int(s), collapse=":"))

warn <- function(..., immediate. = TRUE) warning(..., call. = FALSE, immediate. = immediate.)

# copy objects from one env to another by reference
# see https://stackoverflow.com/questions/9965577/copy-move-one-environment-to-another
copy_objects <- function(from, to, names=ls(from, all.names=TRUE)) {
  mapply(assign, names, mget(names, from), list(to), SIMPLIFY=FALSE, USE.NAMES=FALSE)
  invisible(NULL)
}

is_integer_scalar <- function(x) is.integer(x) && length(x) == 1L
is_numeric_scalar <- function(x) any(typeof(x) == c("integer", "double")) && length(x) == 1L
is_character_scalar <- function(x) is.character(x) && length(x) == 1L

unlst <- function(x, recursive=TRUE) unlist(x, recursive, FALSE)

s_apply <- function(x, f, ..., USE.NAMES=FALSE) vapply(x, f, "", ..., USE.NAMES=USE.NAMES)
b_apply <- function(x, f, ..., USE.NAMES=FALSE) vapply(x, f, TRUE, ..., USE.NAMES=USE.NAMES)
i_apply <- function(x, f, ..., USE.NAMES=FALSE) vapply(x, f, 0L, ..., USE.NAMES=USE.NAMES)
