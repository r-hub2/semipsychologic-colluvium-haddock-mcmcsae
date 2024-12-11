
has_response <- function(formula) length(formula) == 3L

get_response <- function(formula, data=NULL) {
  tf <- terms(formula, data=data)
  ind <- attr(tf, "response")
  if (ind > 0L) {
    vars <- as.list(attr(tf, "variables"))[-1L]
    eval(vars[[ind]], data, environment(formula))
  } else {
    NULL
  }
}

get_types <- function(mod) {
  if (length(mod))
    s_apply(mod, function(x) match.arg(as.character(x[[1L]]), .mod.specials), USE.NAMES=TRUE)
  else
    NULL
}

has_explicit_intercept <- function(formula) {
  attr(terms(update.formula(formula, . ~ 0 + .)), "intercept") > 0L
}

standardize_formula <- function(formula, default="reg", data=NULL, internal.offset=FALSE) {
  # interpret everything not in special terms as a default component
  tf <- terms(formula, keep.order=TRUE, specials=.mod.specials, data=data)
  # NB ~ . - var does not warn if var is not in data
  idx <- unlst(attr(tf, "specials"))  # variable indices of special terms
  if (length(idx)) {
    fac <- attr(tf, "factors")
    for (i in seq_along(idx)) {
      term.idx <- which(fac[idx[i], ] > 0L)  # translate to term indices
      if (length(term.idx) != 1L) stop("cannot parse formula")
      idx[i] <- term.idx
    }
    remainder <- attr(tf, "term.labels")[-idx]
  } else {
    remainder <- attr(tf, "term.labels")
  }
  e <- environment(formula)
  # build up standardized formula
  if (length(remainder)) {
    if (attr(tf, "intercept") == 0L)
      out <- paste0(default, "(~ 0 +", paste0(remainder, collapse=" + "), ")")
    else
      out <- paste0(default, "(~ ", paste0(remainder, collapse=" + "), ")")
  } else {
    if (has_explicit_intercept(formula) ||
        (!length(idx) && !length(attr(tf, "offset")) && attr(tf, "intercept") == 1L))
      out <- paste0(default, "(~ 1)")
    else
      out <- NULL
  }
  if (length(idx)) {
    funpart <- paste0(attr(tf, "term.labels")[idx], collapse=" + ")
    out <- paste0(c(out, funpart), collapse=" + ")
  }
  if (length(attr(tf, "offset"))) {
    for (o in attr(tf, "offset")) {
      out <- paste0(c(out,
        paste0("mc_",
          paste0(gsub("offset(", "offset(~ I(", deparse(attr(tf, "variables")[[o + 1L]]), fixed = TRUE), ")")
        )), collapse = " + "
      )
    }
  } else if (internal.offset) {
    # internal offset used for Poisson family approximated by negbinomial with large negative offset
    out <- paste0(c(out, "mc_offset(value=0)"), collapse=" + ")
    # the internal offset is set later in samplers.R
  }
  if (attr(tf, "response") > 0L) {
    as.formula(paste0(deparse(attr(tf, "variables")[[attr(tf, "response") + 1L]]), " ~ ", out), env=e)
  } else {
    as.formula(paste0("~ ", out), env=e)
  }
}

get_vars <- function(formula, rhs.only=TRUE) {
  tf <- terms(formula)
  vars <- as.list(attr(tf, "variables"))[-1L]
  if (rhs.only) {
    if (attr(tf, "response") > 0L || !is.null(attr(tf, "offset"))) {
      if (attr(tf, "response") > 0L)
        vars <- vars[-c(attr(tf, "response"), attr(tf, "offset"))]
      else
        vars <- vars[-attr(tf, "offset")]
    }
  }
  vars
}

# use prefix to prevent duplicate names in automatic naming of
# model components in different model parts, e.g. mean and variance model
to_mclist <- function(formula, prefix="") {
  vars <- get_vars(formula)
  if (length(vars)) {
    parnames <- s_apply(vars, function(x) if (is.null(x$name)) NA_character_ else x$name)
    types <- get_types(vars)
    if (prefix == "v") {  # backward compatible naming for vfac, vreg components
      prefix <- ifelse(any(types[is.na(parnames)] == c("reg", "gen")), "v", "")
    }
    parnames[is.na(parnames)] <- paste0(prefix, types[is.na(parnames)], which(is.na(parnames)))
    check_mod_names(parnames)
  }
  mod <- list()
  for (m in seq_along(vars)) {
    mc <- vars[[m]]
    mod[[parnames[m]]] <- mc
  }
  mod
}
