
# TODO
# - allow multi-response family object to be passed using helper function, or as list
# - add init method for multi-response case calling individual families' init methods,
#   used e.g. to process response variable
process_family <- function(family) {
  if (is.function(family)) {
    fam <- family()
    if (!is.list(fam) || !is_character_scalar(fam[["family"]])) stop("unrecognised input for argument 'family'")
    family <- fam[["family"]]
  }
  if (class(family)[1L] == "family") {
    family <- switch(family[["family"]],
      gaussian = f_gaussian(link = family[["link"]]),
      binomial = f_binomial(link = family[["link"]]),
      poisson = f_poisson(link = family[["link"]]),
      Gamma = f_gamma(link = family[["link"]])
    )
  } else if (is_character_scalar(family)) {
    if (startsWith(family, "f_")) family <- substring(family, 3L)
    family <- tolower(family)  # allows e.g. both "gamma" and "Gamma"
    if (family == "multi") stop("multi response family must be specified as a list")
    family <- match.arg(family, c("gaussian", "student_t", "binomial", "negbinomial", "poisson", "multinomial", "gamma", "gaussian_gamma"))
    family <- eval(call(paste0("f_", family)))
  } else {
    if (!is.list(family) || is.null(family[["_raw_"]]))
      stop("'family' must be a correctly specified list, most conveniently ",
        "created by one of the family-object creating functions whose names ",
        "start with 'f_', see the help for functions 'f_gaussian', 'f_binomial', etc.; ",
        "use 'f_multi' to specify multiple response families")
    if (family[["family"]] == "multi") {
      for (f in seq_along(family[["fam.list"]])) {
        family$fam.list[[f]] <- process_family(family$fam.list[[f]])
      }
    }
  }
  family
}
