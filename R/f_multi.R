
f_multi <- function(fam.list, fam.ind) {
  if (!is.list(fam.list) || length(fam.list) == 0L)
    stop("'fam.list' must be a (non-empty) list")
  if (is.null(names(fam.list)) || any_duplicated(names(fam.list)) || any(names(fam.list) == ""))
    stop("'fam.list' must have unique (non-empty) names")
  if (!inherits(fam.ind, "formula")) stop("'fam.ind' must be specified as a formula")
  list(family="multi", link="identity", fam.list=fam.list, fam.ind=fam.ind,
       `_raw_`=TRUE)
}

ff_multi <- function(link, fam.list, fam.ind, sm, data, y=NULL) {
  family <- "multi"
  linkinv <- identity
  f_mean <- identity  # mean function acting on linear predictor
  e.is.res <- FALSE
  sigma.fixed <- TRUE
  prior.only <- is.null(y)
  n <- n_row(data)
  prior.only <- is.null(y)

  fam.fac <- fdroplevels(as.factor(get_var_from_formula(fam.ind, data)))
  if (!all(levels(fam.fac) %in% names(fam.list)))
    stop("some levels of fam.ind variable have no corresponding entry in fam.list")

  for (f in seq_along(fam.list)) {
    fam <- fam.list[[f]]
    fam$sm <- sm
    sub <- which(fam.fac == names(fam.list)[f])
    fam$data <- data[sub, , drop=FALSE]
    if (!prior.only) fam$y <- y[sub]
    fam$famid <- names(fam.list)[f]
    fam$sub <- sub
    fam$`_raw_` <- NULL
    fam.list[[f]] <- do.call(
      getFromNamespace(paste0("ff_", fam[["family"]]), "mcmcsae"),
      fam[-1L], envir=environment(fam.ind)
    )
    if (!prior.only) y[sub] <- fam.list[[f]][["y"]]
  }
  sc <- sm[["control"]]
  store_default <- function() {
    out <- fam.list[[1L]]$store_default()
    for (fam in fam.list[-1L]) out <- c(out, fam$store_default())
    out
  }

  #modeled.Q <- any(sapply(fam.list, `[[`, "modeled.Q"))
  modeled.Q <- TRUE
  if (any(sapply(fam.list, `[[`, "Q0.type") == "symm")) {
    Q0.type <- "symm"
  } else if (any(sapply(fam.list, `[[`, "Q0.type") == "diag")) {
    Q0.type <- "diag"
  } else {
    Q0.type <- "unit"
  }
  # TODO bdiag method for ddi matrices
  Q0 <- economizeMatrix(do.call(bdiag, lapply(fam.list, `[[`, "Q0")), sparse=TRUE, symmetric=TRUE)

  scale.e <- 1  # TODO some weighted average of fam.list scale.e's
  scale.sigma <- scale.e

  Q_e <- function(p) {
    out <- numeric(n)
    for (fam in fam.list) out[fam[["sub"]]] <- fam$Q_e(p)
    out
  }

  rprior <- function(p) {
    for (fam in fam.list)
      if (is.function(fam$rprior)) p <- fam$rprior(p)
    p
  }
  MHpars <- NULL
  for (fam in fam.list) MHpars <- c(MHpars, fam[["MHpars"]])
  if (any(sapply(fam.list, function(fam) is.function(fam[["adapt"]])))) {
    adapt <- function(ar)
      for (fam in fam.list)
        if (is.function(fam$adapt)) fam$adapt(ar)
  }
  draw <- function(p) {
    p$llh_ <- 0
    for (fam in fam.list)
      if (is.function(fam$draw)) p <- fam$draw(p)
    p
  }
  start <- function(p) {
    if (modeled.Q)
      p$Q_ <- check_and_get(p, "Q_", n, \() runif(n, 0.9, 1.1), pos=TRUE)
    for (fam in fam.list)
      if (is.function(fam$start)) p <- fam$start(p)
    p
  }

  llh <- function(p) {
    out <- 0
    for (fam in fam.list) out <- out + fam$llh(p)
    out
  }
  llh_i <- function(draws, i, e_i) {
    # TODO pass subsetted i and e_i to families' llh_i
    out <- fam.list[[1L]]$llh_i(draws, i, e_i)
    for (fam in fam.list[-1L]) out <- out + fam$llh_i(draws, i, e_i)
    out
  }
  # weights: passed from predict.mcdraws
  #   can be either a numeric scalar, or a vector of length n or nrow(newdata) if the latter is provided
  make_rpredictive <- function(newdata, weights=NULL) {
    if (is.integer(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      nn <- newdata
      nfam.fac <- fam.fac
      # TODO error if nn != n ?
    } else {
      nn <- nrow(newdata)
      nfam.fac <- as.factor(get_var_from_formula(fam.ind, newdata))
      if (!all(levels(nfam.fac) %in% names(fam.list)))
        stop("unexpected levels of response family variable")
    }
    rpredictive.list <- NULL
    fam.sub <- NULL
    for (f in levels(nfam.fac)) {
      sub <- which(nfam.fac == f)
      rpredictive.list[[f]] <- fam.list[[f]]$make_rpredictive(
        if (is.integer(newdata)) length(sub) else newdata[sub, , drop=FALSE],
        if (is.null(weights)) NULL else if (length(weights) == 1L) weights else weights[sub]
      )
      fam.sub[[f]] <- sub
    }
    function(p, lp) {
      out <- numeric(nn)
      for (f in levels(nfam.fac))
        out[fam.sub[[f]]] <- rpredictive.list[[f]](p, lp[fam.sub[[f]]])
      out
    }
  }

  rm(data, sm, sub, fam, f)
  environment()
}
