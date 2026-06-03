
has_response <- function(formula) length(formula) == 3L

get_response <- function(formula, data=NULL) {
  tf <- terms(formula, data=data)
  ind <- attr(tf, "response")
  if (ind > 0L) {
    vars <- as.list(attr(tf, "variables"))[-1L]
    y <- eval(vars[[ind]], if (is_integer_scalar(data)) NULL else data, environment(formula))
    if (is.numeric(y))
      if (!all(is.finite(y))) stop(sum(!is.finite(y)), " missing or infinite value(s) in response variable ", deparse(vars[[ind]]))
    else
      if (anyNA(y)) stop(sum(is.na(y)), " missing value(s) in response variable ", deparse(vars[[ind]]))
    y
  } else {
    NULL
  }
}

get_types <- function(mod) {
  if (length(mod))
    s_apply(mod, \(x) match.arg(as.character(x[[1L]]), .mod.specials), USE.NAMES=TRUE)
  else
    NULL
}

has_explicit_intercept <- function(formula) {
  attr(terms(update.formula(formula, . ~ 0 + .)), "intercept") > 0L
}

standardise_formula <- function(formula, default="reg", data=NULL, internal.offset=FALSE) {
  # interpret everything not in special terms as a default component
  tf <- terms(formula, keep.order=TRUE, specials=.mod.specials, data=data)
  term.labels <- attr(tf, "term.labels")
  idx <- unlst(attr(tf, "specials"))  # variable indices of special terms
  if (length(idx)) {
    fac <- attr(tf, "factors")
    for (i in seq_along(idx)) {
      term.idx <- whichv(fac[idx[i], ] > 0L, TRUE)  # translate to term indices
      if (length(term.idx) != 1L) stop("cannot parse formula")
      idx[i] <- term.idx
    }
    remainder <- term.labels[-idx]
  } else {
    remainder <- term.labels
  }
  # look for lme4-style random effect specifications
  idx.lme4.re <- grep("|", term.labels, fixed = TRUE)
  if (length(idx.lme4.re))
    remainder <- setdiff(remainder, term.labels[idx.lme4.re])
  e <- environment(formula)
  # build up standardised formula
  if (length(remainder)) {
    if (attr(tf, "intercept") == 0L)
      out <- paste0(default, "(~ 0 +", paste0(remainder, collapse=" + "), ")")
    else
      out <- paste0(default, "(~ ", paste0(remainder, collapse=" + "), ")")
  } else {
    if (has_explicit_intercept(formula) ||
        (!length(idx) && !length(idx.lme4.re) && !length(attr(tf, "offset")) && attr(tf, "intercept") == 1L))
      out <- paste0(default, "(~ 1)")
    else
      out <- NULL
  }
  if (length(idx)) {
    funpart <- paste0(term.labels[idx], collapse=" + ")
    out <- paste0(c(out, funpart), collapse=" + ")
  }
  if (length(idx.lme4.re)) {
    for (i in idx.lme4.re) {
      diag.var <- grepl("||", term.labels[i], fixed = TRUE)
      eff <- strsplit(term.labels[i], if (diag.var) "||" else "|", fixed=TRUE)[[1L]]
      # expand multiple grouping factors and drop empty levels as lme4 would do
      tryCatch(
        grps <- attr(terms(as.formula(paste("~", eff[2L]))), "term.labels"),
        error = function(e) stop("unsupported model term ", term.labels[i], call.=FALSE)
      )
      re.term <- paste0(
        "gen(formula = ~ ", trimws(eff[[1L]]), if (diag.var) ", var = 'diagonal'" else NULL,
        ", drop.empty.levels=TRUE, factor = ~ ", grps, ")", collapse=" + "
      )
      out <- paste0(c(out, re.term), collapse=" + ")
    }
  }
  if (length(attr(tf, "offset"))) {
    for (o in attr(tf, "offset")) {
      expr <- attr(tf, "variables")[[o + 1L]]
      if (is.call(expr[[2L]]) && identical(expr[[2L]][[1L]], as.name("~")))
        stop("formulas are not allowed in offset", call. = FALSE)
      out <- paste0(c(out,
        paste0("mc_",
          paste0(gsub("offset(", "offset(~ I(", deparse(expr), fixed = TRUE), ")")
        )), collapse = " + "
      )
    }
  } else if (internal.offset) {
    # internal offset used for Poisson family approximated by negbinomial with large negative offset
    out <- paste0(c(out, "mc_offset(value=0)"), collapse=" + ")
    # the internal offset is set later in samplers.R
  }
  if (attr(tf, "response") > 0L) {
    as.formula(paste0(deparse(attr(tf, "variables")[[attr(tf, "response") + 1L]]),
      " ~ ", if (is.null(out)) "0" else out), env=e)
  } else {
    as.formula(paste0("~ ", if (is.null(out)) "0" else out), env=e)
  }
}

get_var_from_formula <- function(f, data) {
  tf <- terms(f, data=data)
  vars <- as.list(attr(tf, "variables"))[-1L]
  if (length(vars) != 1L) stop("formula with single variable expected, but found ", length(vars))
  out <- eval(vars[[1L]], if (is_integer_scalar(data)) NULL else data, environment(f))
  if (inherits(out, "AsIs")) out <- unclass(out)
  if (is.factor(out) || is.character(out)) {
    if (anyNA(out)) stop("total of ", sum(is.na(out)), " NAs in variable ", vars)
  } else {
    if (!all(is.finite(out))) stop("total of ", sum(!is.finite(out)), " NA/NaN/Inf in variable ", vars)
  }
  out
}

# use prefix to prevent duplicate names in automatic naming of
# model components in different model parts, e.g. mean and variance model
to_mclist <- function(formula, prefix="") {
  tf <- terms(formula)
  vars <- as.list(attr(tf, "variables"))[-1L]
  # drop lhs and offset(s)
  drp <- attr(tf, "offset")  # NULL if no offset
  if (attr(tf, "response") > 0L) drp <- c(attr(tf, "response"), drp)
  if (length(drp)) vars <- vars[-drp]
  if (length(vars)) {
    parnames <- s_apply(vars, \(x) if (is.null(x[["name"]])) NA_character_ else x[["name"]])
    types <- get_types(vars)
    if (prefix == "v") {  # backward compatible naming for vfac, vreg components
      prefix <- ifelse(any(types[is.na(parnames)] == c("reg", "gen")), "v", "")
    }
    parnames[is.na(parnames)] <- paste0(prefix, types[is.na(parnames)], whichNA(parnames))
    check_mod_names(parnames)
  }
  mod <- list()
  for (m in seq_along(vars)) {
    mc <- vars[[m]]
    mod[[parnames[m]]] <- mc
  }
  mod
}
