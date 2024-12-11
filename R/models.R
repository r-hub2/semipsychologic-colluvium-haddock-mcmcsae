
.mod.specials <- c("reg", "gen", "mec", "brt", "vreg", "vfac", "mc_offset")

#' Check names of model components
#'
#' @noRd
#' @param x vector of names of model components.
#' @returns \code{TRUE} if all names are OK; throws an error otherwise.
check_mod_names <- function(x) {
  if (any(startsWith(x, "_")) || any(endsWith(x, "_"))) stop("\'_\' as first or last character of a model component name is reserved for internal use")
  # check for names ending with extensions like _sigma, _rho, _xi etc since they might clash
  reserved.exts <- c("_df", "_gl", "_GMRFext", "_omega", "_rho", "_sigma", "_xi")
  if (any(grepl(paste0("(", paste0(reserved.exts, collapse="|"), ")$"), x)))
    stop("model component names ending with any of ('", paste0(reserved.exts, collapse="', '"), "') are not allowed")
  if (any(duplicated(x))) stop("duplicate model component name(s): '", paste0(x[duplicated(x)], collapse="', '"), "'")
  TRUE
}

#' Compute a list of design matrices for all terms in a model formula,
#' or based on a sampler environment
#'
#' If \code{sampler} is provided instead of \code{formula}, the design matrices
#' are based on the model used to create the sampler environment. In that case, if
#' \code{data} is \code{NULL}, the design matrices stored in \code{sampler} are returned,
#' otherwise the design matrices are computed for the provided data based on the sampler's model.
#' The output is a list of dense or sparse design matrices for the model components
#' with respect to \code{data}.
#'
#' @examples
#' n <- 1000
#' dat <- data.frame(
#'   x = rnorm(n),
#'   f = factor(sample(1:50, n, replace=TRUE))
#' )
#' str(computeDesignMatrix(~ x, dat)[[1]])
#' model <- ~ reg(~x, name="beta") + gen(~x, factor=~f, name="v")
#' X <- computeDesignMatrix(model, dat)
#' str(X)
#'
#' @export
#' @param formula model formula.
#' @param data data frame to be used in deriving the design matrices.
#' @param labels if \code{TRUE}, column names are assigned.
#' @returns A list of design matrices.
computeDesignMatrix <- function(formula=NULL, data=NULL, labels=TRUE) {
  if (!inherits(formula, "formula")) stop("'formula' must be a formula")
  formula <- standardize_formula(formula, data=data)
  out <- to_mclist(formula)
  for (m in seq_along(out)) {
    out[[m]] <- match.call(eval(out[[m]][[1L]]), out[[m]])
    mc <- as.list(out[[m]])[-1L]
    if (is.null(mc$formula))
      mc$formula <- ~ 1
    else
      mc$formula <- as.formula(eval(mc$formula))
    environment(mc$formula) <- environment(formula)
    if (any("|" == all.names(mc$formula))) {
      # assume this is a mec component; we use as design matrix the covariate matrix subject to error
      vs <- as.list(attr(terms(mc$formula), "variables")[-1L])
      formula.X <- as.formula(paste0("~ 0 + ", paste0(s_apply(vs, function(x) deparse(x[[2L]])), collapse=" + ")), env=environment(formula))
      out[[m]] <- model_matrix(formula.X, data, sparse=mc$sparse)
    } else {
      if (!is.null(mc$factor))
        mc$factor <- as.formula(mc$factor, env=environment(formula))
      if (is.null(mc$remove.redundant)) mc$remove.redundant <- FALSE
      if (is.null(mc$drop.empty.levels)) mc$drop.empty.levels <- FALSE
      out[[m]] <- compute_X(mc$formula, mc$factor, mc$remove.redundant,
          mc$drop.empty.levels, mc$sparse, data)
    }        
    if (!labels) colnames(out[[m]]) <- NULL
  }
  out
}

#' Compute design matrix for a model component
#'
#' The output is a dense or sparse design matrix for the specified model component with respect to
#' the specified data frame. The design matrix contains appropriate column names.
#'
#' @noRd
#' @param formula the formula part of a model component.
#' @param factor the factor part of a model component (also a formula object).
#' @param remove.redundant whether to remove redundant columns.
#' @param drop.empty.levels whether to remove factor levels without observations.
#' @param sparse whether the design matrix based on \code{formula} alone should be in a sparse matrix format.
#' @param data data frame to be used in deriving the design matrix.
#' @returns A design matrix including column names.
compute_X <- function(formula=~1, factor=NULL,
                      remove.redundant=TRUE, drop.empty.levels=FALSE,
                      sparse=NULL, data=NULL) {
  X0 <- model_matrix(formula, data, sparse=sparse)
  if (remove.redundant) X0 <- remove_redundancy(X0)
  if (is.null(factor) || (inherits(factor, "formula") && intercept_only(factor)))
    return(X0)
  XA <- compute_XA(factor, data)
  if (drop.empty.levels) {
    cols2remove <- which(zero_col(XA))
    if (length(cols2remove))
      XA <- XA[, -cols2remove]
  } else {
    cols2remove <- NULL
  }
  X <- combine_X0_XA(X0, XA)
  if (length(cols2remove))
    attr(X, "factor.cols.removed") <- cols2remove
  if (ncol(X0) > 1L)
    attr(X, "formula.colnames") <- colnames(X0)
  X
}

combine_X0_XA <- function(X0, XA) {
  n <- nrow(X0)
  q0 <- ncol(X0)
  qA <- ncol(XA)
  if ((is_unit_diag(XA) || class(XA)[1L] == "tabMatrix") && !large_and_sparse(X0)) {
    fac <- if (is_unit_diag(XA)) seq_len(n) else XA@perm + 1L
    x <- X0
    if (class(XA)[1L] == "tabMatrix") {
      if (XA@num) x <- x * XA@x
      if (XA@reduced) x <- x * as.numeric(XA@perm != -1L)
    }
    out <- drop0(sparseMatrix(
        i = rep.int(seq_len(n), q0),
        j = rep.int(q0 * (fac - 1L) + 1L, q0) + rep_each(0L:(q0 - 1L), n),
        x = as.numeric(x),  # inefficient for large sparse X0
        dims = c(n, q0 * qA), check=FALSE
      ), is.Csparse=TRUE
    )
  } else {
    # use Matrix::KhatriRao as it is more memory-efficient
    out <- t(KhatriRao(t(XA), t(X0)))  # a dgCMatrix
  }
  if (q0 > 1L) {
    if (qA * q0 > 1e6L) {
      labs <- NULL
    } else {  # this may be too much; forget about column names
      labs <- paste(rep.int(colnames(X0), qA), rep_each(colnames(XA), q0), sep=":")
    }
  } else {
    labs <- colnames(XA)
  }
  attr(out, "Dimnames") <- list(NULL, labs)
  economizeMatrix(out, strip.names=FALSE)
}

compute_XA <- function(factor=NULL, data=NULL) {
  enclos <- environment(factor)
  factor.info <- get_factor_info(factor, data)
  if (is.null(factor.info)) return(NULL)
  n <- n_row(data)
  if (any("spline" == factor.info$types)) {  # B-spline design matrix components
    fs <- as.list(attr(terms(factor), "variables"))[-1L]
    out <- matrix(1, nrow=1L, ncol=n)
    labs <- ""
    for (f in seq_along(factor.info$types)) {
      variable <- eval_in(factor.info$variables[f], data, enclos)
      switch(factor.info$types[f],
        spline = {
          # P-splines, i.e. penalized B-splines
          # RW1/2(t) are special cases with degree 1, knots (tmin+1):(tmax-1), and the same precision matrix
          if (is.null(fs[[f]][["knots"]])) fs[[f]]$knots <- min(variable) + ((-4:35)/30) * (max(variable) - min(variable))
          if (is.null(fs[[f]][["degree"]])) fs[[f]]$degree <- 3L
          if (length(fs[[f]][["knots"]]) == 1L) {  # assumed to be the number of knots
            if (fs[[f]][["knots"]] <= 2L * (fs[[f]][["degree"]] + 2L)) stop("spline: too few knots")
            fs[[f]]$knots <- min(variable) + ((seq_len(fs[[f]][["knots"]]) - (fs[[f]][["degree"]] + 1L))/(fs[[f]][["knots"]] - 2L*(fs[[f]][["degree"]] + 2L))) * (max(variable) - min(variable))
          }
          Xf <- drop0(splines::splineDesign(fs[[f]][["knots"]], variable, fs[[f]][["degree"]] + 1L, sparse=TRUE), tol=sqrt(.Machine$double.eps), is.Csparse=TRUE)
          labs <- as.vector(outer(labs, paste0("bs", seq_len(ncol(Xf))), FUN=paste, sep=if (identical(labs, "")) "" else ":"))
        },
        {
          Xf <- aggrMatrix(variable, facnames=TRUE)
          labs <- as.vector(outer(labs, colnames(Xf), FUN=paste, sep=if (identical(labs, "")) "" else ":"))
        }
      )
      out <- KhatriRao(t(Xf), out)
    }
    out <- economizeMatrix(t(out), strip.names=FALSE)  # names of @i and @x removed by t()
    dimnames(out) <- list(NULL, labs)
    out
  } else {
    fac <- combine_factors(factor.info$variables, data, enclos=enclos)
    if (anyNA(fac)) stop("NA's in 'factor' not allowed")
    aggrMatrix(fac, facnames=TRUE)
  }
}

# extract information from the factor formula component of a model component list
get_factor_info <- function(formula, data) {
  enclos <- environment(formula)
  if (is.null(formula) || intercept_only(formula)) return(NULL)
  fs <- as.list(attr(terms(formula), "variables"))[-1L]
  fs <- lapply(fs, function(x) if (is.symbol(x)) call("iid", x) else x)
  variables <- vapply(fs, function(x) deparse(x[[2L]]), "")
  types <- vapply(fs, function(x) as.character(x[[1L]]), "")
  ncols <- integer(length(fs))
  # extra info, currently used to hold AR1 prior and MH options
  extra <- vector("list", length(fs))
  for (f in seq_along(ncols)) {
    if (types[f] == "spline") {
      if (is.null(fs[[f]][["knots"]])) {
        n.knots <- 40L
      } else {
        if (length(fs[[f]][["knots"]]) == 1L)
          n.knots <- fs[[f]][["knots"]]
        else
          n.knots <- length(eval(fs[[f]][["knots"]]))
      }
      if (is.null(fs[[f]][["degree"]])) fs[[f]]$degree <- 3L
      ncols[f] <- n.knots - fs[[f]][["degree"]] - 1L
    } else {
      if (variables[[f]] == "local_")
        ncols[f] <- n_row(data)
      else
        ncols[f] <- nlevels(as.factor(eval_in(variables[f], data, enclos)))
      if (types[f] == "AR1") {
        fs[[f]] <- match.call(function(variable, phi=NULL, w=NULL, control=NULL) {}, fs[[f]])
        prior.AR1 <- eval(fs[[f]][["phi"]])
        if (is.null(prior.AR1)) prior.AR1 <- pr_unif(-1, 1)
        if (is.environment(prior.AR1)) {
          if (prior.AR1$type == "fixed") {
            fs[[f]]$phi <- prior.AR1$value
          } else {
            prior.AR1$init(1L)
            control.AR1 <- eval(fs[[f]][["control"]])
            extra[[f]] <- list(prior = prior.AR1, control = control.AR1)
          }
        } else {
          if (!is_numeric_scalar(prior.AR1)) stop("AR1 argument 'phi' must either be a single numeric value or a prior specification")
        }
      }
    }
  }
  out <- list(factors=fs, types=types, variables=variables, n=ncols, extra=extra)
  attr(out, ".Environment") <- enclos
  out
}

eval_in <- function(text, data, enclos=emptyenv()) {
  if (text == "local_" && all("local_" != colnames(data))) return(seq_len(n_row(data)))
  if (text == "global_" && all("global_" != colnames(data))) return(rep.int(1, n_row(data)))
  if (is.integer(data) && length(data) == 1L)
    eval(substitute(evalq(expr, NULL, enclos), list(expr=str2lang(text))))
  else
    eval(substitute(evalq(expr, data, enclos), list(expr=str2lang(text))))
}

#' Compute (I)GMRF incidence, precision and restriction matrices corresponding to a generic model component
#'
#' This function computes incidence, precision and restriction matrices, or
#' a subset thereof, for a Gaussian Markov Random Field (GMRF).
#' A GMRF is specified by a formula passed to the \code{factor} argument,
#' in the same way as for the \code{factor} argument of \code{\link{gen}}.
#'
#' @examples
#' n <- 1000
#' dat <- data.frame(
#'   x = rnorm(n),
#'   f1 = factor(sample(1:50, n, replace=TRUE)),
#'   f2 = factor(sample(1:10, n, replace=TRUE))
#' )
#' mats <- compute_GMRF_matrices(~ f1 * RW1(f2), dat)
#' str(mats)
#'
#' @noRd
#' @param factor factor formula of a generic model component,
#'  see \code{\link{gen}}. Alternatively a factor info object as returned by get_factor_info.
#' @param data data frame to be used in deriving the matrices. Not required if \code{factor} is
#'  a factor info object.
#' @param D if \code{TRUE} compute the incidence matrix.
#' @param Q if \code{TRUE} compute the precision matrix.
#' @param R if \code{TRUE} compute the restriction matrix.
#' @param cols2remove if an integer vector is passed, the dimensions (columns of D,
#'  rows and columns of Q and rows of R) that are removed. This can be useful in the
#'  case of empty domains.
#' @param remove.redundant.R.cols whether to test for and remove redundant
#'  restrictions from restriction matrix R.
#' @param scale.precision whether to scale incidence and precision matrices such
#'  that the corresponding marginal variances are close to 1.
#' @param ... further arguments passed to \code{economizeMatrix}.
#' @returns A list containing some or all of the components \code{D} (incidence matrix),
#'  \code{Q} (precision matrix) and \code{R} (restriction matrix).
compute_GMRF_matrices <- function(factor, data, D=TRUE, Q=TRUE, R=TRUE, cols2remove=NULL,
                                  remove.redundant.R.cols=TRUE, scale.precision=FALSE, ...) {
  out <- list()
  if (!D && !Q && !R) return(out)
  if (D || Q) DQ <- CdiagU(1L)
  if (is.null(factor)) {
    if (D) out$D <- DQ
    if (Q) out$Q <- DQ
    return(out)
  }
  if (!D && !Q) scale.precision <- FALSE
  needD <- D || scale.precision
  needR <- R || scale.precision
  if (needR) {
    out$R <- NULL
    nIGMRF <- 0L  # number of (singular) IGMRF factors
  }
  if (inherits(factor, "formula"))
    info <- get_factor_info(factor, data)
  else
    info <- factor
  env <- environment(factor)
  for (f in seq_along(info[["types"]])) {
    fcall <- info$factors[[f]]
    if (inherits(fcall, "name")) {
      fcall <- call("iid", fcall)  # if no GMRF type name is used assume iid
    }
    fcall[[2L]] <- NULL  # remove first argument, which is assumed to be the factor variable's name
    if (info$types[f] != "spatial") fcall$n <- info$n[f]
    if (needD || Q) {
      DQcall <- fcall
      DQcall[[1L]] <- as.name(paste0(if (needD) "D_" else "Q_", DQcall[[1L]]))
    }
    if (needR) {
      Rcall <- fcall
      Rcall[[1L]] <- as.name(paste0("R_", Rcall[[1L]]))
    }
    switch(info$types[f],
      spline = {
        penalty <- fcall[["penalty"]]
        if (is.null(penalty)) penalty <- "RW2"
        if (all(penalty != c("iid", "RW1", "RW2"))) stop("unsupported spline penalty argument")
        if (needD || Q) {
          DQcall[[1L]] <- as.name(paste0(if (needD) "D_" else "Q_", penalty))
          DQcall$knots <- DQcall$degree <- DQcall$penalty <- NULL
          DQf <- eval(DQcall)
        }
        if (needR) {
          if (penalty == "iid") {
            Rf <- NULL
          } else {
            Rcall[[1L]] <- as.name(paste0("R_", penalty))
            Rcall$knots <- Rcall$degree <- Rcall$penalty <- NULL
            Rf <- eval(Rcall)
          }
        }
      },
      spatial = {
        deprecated <- FALSE
        if (needD || Q) {
          DQcall <- match.call(spatial, DQcall)
          if (!is.null(DQcall[["graph"]]))
            graph <- eval(DQcall[["graph"]], env)
          else if (!is.null(DQcall[["poly.df"]])) {
            graph <- eval(DQcall[["poly.df"]], env)
            deprecated <- TRUE
          } else
            graph <- data
          snap <- eval(DQcall[["snap"]], env)
          queen <- eval(DQcall[["queen"]], env)
        } else {
          Rcall <- match.call(spatial, Rcall)
          if (!is.null(Rcall[["graph"]]))
            graph <- eval(Rcall[["graph"]], env)
          else if (!is.null(Rcall[["poly.df"]])) {
            graph <- eval(Rcall[["poly.df"]], env)
            deprecated <- TRUE
          } else
            graph <- data
          snap <- eval(Rcall[["snap"]], env)
          queen <- eval(Rcall[["queen"]], env)
        }
        if (deprecated) warn("argument 'poly.df' is deprecated; please use argument 'graph' instead")
        graph <- get_neighbours_list(graph, snap, queen)
        if (length(graph) != info$n[f]) stop("number of nodes of 'graph', ", length(graph), ", does not equal number of levels of variable '", info$variables[f], "', ", info$n[f])
        # TODO
        #   aggregate spatial object to factor, if the spatial object has more detail
        #   check that level order is the same: first try variable of the same name as factor, then try rownames of @data; if unsuccessful issue a warning
        if (needD)
          DQf <- D_spatial(graph)
        else if (Q)
          DQf <- Q_spatial(graph)
        if (needR) {
          if (isTRUE(Rcall$derive.constraints)) warn("argument 'derive.constraints' is deprecated as it is no longer needed")
          Rf <- R_spatial(graph)
        }
      },
      custom = {
        if (needD || Q) {
          DQcall <- match.call(custom, DQcall)
          if (needD) {
            if (is.null(DQcall[["D"]])) stop("missing argument 'D' in custom factor")
            DQf <- economizeMatrix(eval(DQcall[["D"]], env), sparse=TRUE)
            if (!is.null(DQcall[["Q"]])) warn("argument 'Q' in custom() ignored")
          } else if (is.null(DQcall[["Q"]])) {
            if (is.null(DQcall[["D"]])) stop("custom factor must specify either incidence matrix 'D' or precision matrix 'Q'")
            DQf <- economizeMatrix(crossprod(eval(DQcall[["D"]], env)), sparse=TRUE)
          } else {
            DQf <- economizeMatrix(eval(DQcall[["Q"]], env), sparse=TRUE, symmetric=TRUE)
          }
          if (ncol(DQf) != info$n[f]) stop("custom factor matrix has unexpected number of columns")
        }
        if (needR) {
          Rcall <- match.call(custom, Rcall)
          if (is.null(Rcall[["R"]])) {
            if (isTRUE(Rcall[["derive.constraints"]])) {
              if (!needD && !Q) warn("cannot derive constraints as custom incidence or precision matrix unavailable")
              Rf <- derive_constraints(if (needD) crossprod(DQf) else DQf)
            } else {
              Rf <- NULL
            }
          } else {
            Rf <- economizeMatrix(eval(Rcall[["R"]], env), allow.tabMatrix=FALSE)
            if (dim(Rf)[1L] != info$n[f]) stop("wrong dimension of custom restriction matrix")
            if (isTRUE(Rcall[["derive.constraints"]])) warn("argument 'derive.constraints' ignored as restriction matrix 'R' is specified")
          }
        }
      },
      AR1 = {
        if (needD || Q) {
          if (!is_numeric_scalar(info$factors[[f]]$phi) && !is.null(info$extra[[f]])) {
            # inferred AR1 parameter, use value 0.5 for template construction
            DQcall$phi <- 0.5
          }
          DQcall$control <- NULL
          DQf <- eval(DQcall)
        }
        if (needR) Rf <- NULL
      },
      {
        if (needD || Q) DQf <- eval(DQcall)
        if (needR) Rf <- if (info$types[f] == "iid") NULL else eval(Rcall)
      }
    )
    if (needD || Q) DQ <- cross(DQ, DQf)
    if (needR) {
      if (!is.null(out$R)) {
        out$R <- kronecker(CdiagU(info$n[f]), out$R)
      }
      if (!is.null(Rf)) {
        nIGMRF <- nIGMRF + 1L
        if (is.null(out$R))
          out$R <- zeroMatrix(prod(info$n[seq_len(f)]), 0L)
        out$R <- cbind(out$R, kronecker(Rf, CdiagU(prod(info$n[seq_len(f - 1L)]))))
      }
    }
  }  # END for (f in seq_along(factor.info$types))
  if (needR && !is.null(out$R)) {
    if (!is.null(cols2remove)) out$R <- out$R[-cols2remove, ]
    if ((remove.redundant.R.cols && nIGMRF >= 2L) || !is.null(cols2remove)) {
      # for multiple IGMRF factors R as constructed has redundant columns,
      #   because of duplicate inclusion of cross-poducts of null-vectors
      out$R <- remove_redundancy(out$R)
    }
    out$R <- economizeMatrix(out$R, allow.tabMatrix=FALSE, ...)
  }
  if (needD) {
    if (!is.null(cols2remove)) {
      DQ <- DQ[, -cols2remove]
      # then see which rows become all-zero and remove them too
      DQ <- DQ[-which(rowSums(DQ * DQ) == 0), ]
    }
    out$D <- economizeMatrix(DQ, ...)
    if (scale.precision && !is_unit_diag(out$D)) {
      # scale with geometric mean of marginal variances, as proposed in Riebler ea
      mvar <- sim_marg_var(D=out$D, R=out$R)
      if (!is.null(out$R) & ncol(out$R) > 1L) {
        # assume that each column of R corresponds to a connected component of the graph
        comp <- apply(out$R, 2L, function(x) which(x != 0))
        mvar <- lapply(comp, function(ind) exp(mean(log(mvar[ind]))))
        f <- numeric(ncol(out$D))
        for (i in seq_along(comp)) {
          if (length(comp[[i]]) == 1L) {
            # deal with singletons: set precision to 1
            e_i <- numeric(ncol(out$D))
            e_i[comp[[i]]] <- 1
            out$D <- rbind(out$D, e_i)
            f[comp[[i]]] <- 1
          } else {
            f[comp[[i]]] <- mvar[[i]]
          }
        }
        if (any(i_apply(comp, length) == 1L)) {
          # associate independent random effects with singletons, do not constrain to 0
          out$R <- economizeMatrix(
            out$R[, -which(i_apply(comp, length) == 1L), drop=FALSE],
            allow.tabMatrix=FALSE, ...
          )
        }
        out$D <- economizeMatrix(out$D %*% Cdiag(f), ...)
      } else {
        mvar <- exp(mean(log(mvar)))
        if (mvar > 0)
          out$D <- sqrt(mvar) * out$D
        else warn("zero marginal variance detected; will not scale precision matrix")
      }
    }
  }
  if (Q) {
    if (needD) {
      out$Q <- crossprod(out$D)
    } else {
      out$Q <- DQ
      if (!is.null(cols2remove)) out$Q <- out$Q[-cols2remove, -cols2remove]
    }
    out$Q <- economizeMatrix(out$Q, symmetric=TRUE, ...)
  }
  if (needD && !D) out$D <- NULL
  if (needR && !R) out$R <- NULL
  out
}

#######
# functions to compute sparse precision matrices Q_, incidence matrices D_
#   and for IGMRFs associated constraint matrices R_


#' Correlation factor structures in generic model components
#'
#' Element 'factor' of a model component created using function
#' \code{\link{gen}} is a formula composed of several possible terms described
#' below. It is used to derive a (typically sparse) precision matrix for a set of
#' coefficients, and possibly a matrix representing a set of linear constraints
#' to be imposed on the coefficient vector.
#' \describe{
#'   \item{iid(f)}{Independent effects corresponding to the levels of factor \code{f}.}
#'   \item{RW1(f, circular=FALSE, w=NULL)}{First-order random walk over the levels of factor \code{f}.
#'     The random walk can be made circular and different (fixed) weights can be attached to the innovations.
#'     If specified, \code{w} must be a positive numeric vector of length one less than the number of
#'     factor levels. For example, if the levels correspond to different times, it would often be
#'     reasonable to choose \code{w} proportional to the reciprocal time differences. For equidistant
#'     times there is generally no need to specify \code{w}.}
#'   \item{RW2(f)}{Second-order random walk.}
#'   \item{AR1(f, phi, w=NULL, control=NULL)}{First-order autoregressive correlation structure among
#'     the levels of \code{f}. Argument \code{phi} can be a single numerical value of the
#'     autoregressive parameter, or an appropriate prior specification if phi should be inferred.
#'     If not supplied, a uniform prior on (-1, 1] is assumed.
#'     For irregularly spaced AR(1) processes weights can be specified, in the same way as for
#'     \code{RW1}.}
#'   \item{season(f, period)}{Dummy seasonal with period \code{period}.}
#'   \item{spatial(f, graph, snap, queen)}{CAR spatial correlation.
#'     Argument \code{graph} can either be an object of (S4) class \code{SpatialPolygonsDataFrame}
#'     or an object of (S3) class \code{sf}. The latter can be obtained, e.g., by reading in a
#'     shape file using function \code{\link[sf]{st_read}}. Arguments \code{snap} and \code{queen}
#'     are passed to \code{\link[spdep]{poly2nb}}, which computes a neighbours list.
#'     Alternatively, a neighbours list object of class \code{nb} can be passed directly
#'     as argument \code{graph}.}
#'   \item{spline(f, knots, degree)}{P-splines, i.e. penalized B-splines structure over
#'     the domain of a quantitative variable f. Arguments knots and degree are passed to
#'     \code{\link[splines]{splineDesign}}. If \code{knots} is a single value it is interpreted as
#'     the number of knots, otherwise as a vector of knot positions. By default 40 equally spaced
#'     knots are used, and a degree of 3.}
#'   \item{custom(f, D=NULL, Q=NULL, R=NULL, derive.constraints=NULL)}{Either a custom precision or incidence
#'     matrix associated with factor f can be passed to argument \code{Q} or \code{D}. Optionally a
#'     constraint matrix can be supplied as \code{R}, or constraints can be derived from the null space
#'     of the precision matrix by setting \code{derive.constraints=TRUE}.}
#' }
#'
#' @examples
#' \donttest{
#' # example of CAR spatial random effects
#' if (requireNamespace("sf")) {
#'   # 1. load a shape file of counties in North Carolina
#'   nc <- sf::st_read(system.file("shape/nc.shp", package="sf"))
#'   # 2. generate some data according to a model with a few regression
#'   # effects, as well as spatial random effects
#'   gd <- generate_data(
#'     ~ reg(~ AREA + BIR74, prior=pr_normal(precision=1), name="beta") +
#'       gen(factor = ~ spatial(NAME, graph=nc), name="vs"),
#'     sigma.mod = pr_invchisq(df=10, scale=1),
#'     data = nc
#'   )
#'   # add the generated target variable and the spatial random effects to the
#'   # spatial dataframe object
#'   nc$y <- gd$y
#'   nc$vs_true <- gd$pars$vs
#'   # 3. fit a model to the generated data, and see to what extent the
#'   #    parameters used to generate the data, gd$pars, are reproduced
#'   sampler <- create_sampler(
#'     y ~ reg(~ AREA + BIR74, prior=pr_normal(precision=1), name="beta") +
#'     gen(factor = ~ spatial(NAME, graph=nc), name="vs"),
#'     data=nc
#'   )
#'   sim <- MCMCsim(sampler, store.all=TRUE, n.iter=600, n.chain=2, verbose=FALSE)
#'   (summ <- summary(sim))
#'   nc$vs <- summ$vs[, "Mean"]
#'   plot(nc[c("vs_true", "vs")])
#'   plot(gd$pars$vs, summ$vs[, "Mean"]); abline(0, 1, col="red")
#' }
#' }
#'
#' @param name name of a variable, unquoted.
#' @param circular whether the random walk is circular.
#' @param w a vector of weights.
#' @param phi prior distribution, or fixed value, for an
#'  autoregressive parameter. The default is a uniform prior over
#'  the interval [-1, 1]. A single numeric value is interpreted as
#'  a fixed value, corresponding to a degenerate prior, which can also
#'  be specified as \code{pr_fixed(value)}. Alternatively,
#'  \code{link{pr_truncnormal}} can be used to specify a truncated
#'  normal prior.
#' @param control options for Metropolis-Hastings sampling from the
#'  conditional posterior for an autoregressive parameter. These options
#'  can be set using function \code{\link{set_MH}}. Supported proposal
#'  types are "TN" and "RWTN". By default an independence truncated
#'  normal proposal (type="TN"), or a random walk truncated normal
#'  proposal (type="RWTN") with adaptive scale initialised at 0.025
#'  is used, depending on whether the specified random effects'
#'  distribution is Gaussian or non-Gaussian.
#' @param period a positive integer specifying the seasonal period.
#' @param graph either a spatial object of class \code{SpatialPolygons},
#'  \code{sf}, \code{sfc}, or a neighbours list of class \code{nb}.
#' @param snap passed to \code{\link[spdep]{poly2nb}}.
#'  Ignored if \code{graph} is a neighbours list.
#' @param queen passed to \code{\link[spdep]{poly2nb}}.
#'  Ignored if \code{graph} is a neighbours list.
#' @param poly.df a spatial data frame. DEPRECATED, use argument
#'  \code{graph} instead.
#' @param derive.constraints whether to derive the constraint matrix for an
#'  IGMRF model component numerically from the precision matrix.
#'  The use of \code{derive.constraints} in function \code{spatial}
#'  is DEPRECATED, as it is no longer needed.
#' @param knots passed to \code{\link[splines]{splineDesign}}.
#' @param degree passed to \code{\link[splines]{splineDesign}}.
#' @param D custom incidence matrix.
#' @param Q custom precision matrix.
#' @param R custom restriction matrix.
#'
#' @name correlation
#' @references
#'  B. Allevius (2018).
#'    On the precision matrix of an irregularly sampled AR(1) process.
#'    arXiv:1801.03791.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
NULL

dont_call <- function() stop("this function should only be used inside a formula")

#' @export
#' @rdname correlation
iid <- function(name) dont_call()

#' @export
#' @rdname correlation
RW1 <- function(name, circular=FALSE, w=NULL) dont_call()

#' @export
#' @rdname correlation
RW2 <- function(name) dont_call()

#' @export
#' @rdname correlation
AR1 <- function(name, phi, w=NULL, control=NULL) dont_call()

#' @export
#' @rdname correlation
season <- function(name, period) dont_call()

#' @export
#' @rdname correlation
spatial <- function(name, graph=NULL, snap=sqrt(.Machine$double.eps), queen=TRUE,
                    poly.df=NULL, derive.constraints=FALSE) dont_call()

#' @export
#' @rdname correlation
spline <- function(name, knots, degree) dont_call()

#' @export
#' @rdname correlation
custom <- function(name, D=NULL, Q=NULL, R=NULL, derive.constraints=NULL) dont_call()


D_iid <- Q_iid <- function(n) CdiagU(n)

# use weights w for irregular spacing, see Allevius 2018 research report
Q_AR1 <- function(n, phi=NULL, w=NULL) {
  if (n < 2L) stop("AR1 model needs a sequence of at least 2 periods")
  if (!is_numeric_scalar(phi) || phi < -1 || phi > 1) stop("autoregressive parameter must be a single numeric value between -1 and +1")
  if (is.null(w)) {
    bandSparse(n, n, 0:1, list(c(1, rep.int(1 + phi^2, n-2L), 1), rep.int(-phi, n-1L)), symmetric=TRUE)
  } else {
    # w_i to be interpreted as 1/(t_{i+1} - t_{i})
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    phi2 <- phi^2
    iw <- 1/w
    d0 <- numeric(n)
    d0[1L] <- (1 - phi2) / (1 - phi2^iw[1L])
    d0[2:(n-1L)] <- (1 - phi2)*(1 - phi2^(iw[1:(n-2L)] + iw[2:(n-1L)])) / ((1 - phi2^iw[1:(n-2L)])*(1 - phi2^iw[2:(n-1L)]))
    d0[n] <- (1 - phi2) / (1 - phi2^iw[n-1L])
    d1 <- - (1 - phi2) * phi^iw / (1 - phi2^iw)
    bandSparse(n, n, 0:1, list(d0, d1), symmetric=TRUE)
  }
}

D_AR1 <- function(n, phi=NULL, w=NULL) {
  if (n < 2L) stop("AR1 model needs a sequence of at least 2 periods")
  if (!is_numeric_scalar(phi) || phi < -1 || phi > 1) stop("autoregressive parameter must be a single numeric value between -1 and +1")
  if (is.null(w)) {
    as(bandSparse(n, n, 0:1, list(c(rep.int(1, n-1L), sqrt(1 - phi^2)), rep.int(-phi, n-1L))), "generalMatrix")
  } else {
    # w_i to be interpreted as 1/(t_{i+1} - t_{i})
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    phi2 <- phi^2
    iw <- 1/w
    d0 <- c(sqrt((1 - phi2)/(1 - phi2^iw)), sqrt(1 - phi2))
    d1 <- - sqrt((1 - phi2)/(1 - phi2^iw)) * phi^iw
    as(bandSparse(n, n, 0:1, list(d0, d1)), "generalMatrix")
  }
}

Q_RW1 <- function(n, circular=FALSE, w=NULL) {
  if (n < 2L) stop("RW1 model needs a sequence of at least 2 periods")
  if (is.null(w)) {
    if (circular)
      bandSparse(n, n, c(0L, 1L, n-1L), list(rep.int(2, n), rep.int(-1, n-1L), -1L), symmetric=TRUE)
    else
      bandSparse(n, n, 0:1, list(c(1, rep.int(2, n-2L), 1), rep.int(-1, n-1L)), symmetric=TRUE)
  } else {
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    if (circular) {
      bandSparse(n, n, c(0L, 1L, n-1L), list(c(w, w[n-1L]) + c(w[n-1L], w), rep.int(-w, n-1L), -w[n-1L]), symmetric=TRUE)
    } else {
      d0 <- c(w, w[n-1L])
      d0[2:(n-1L)] <- d0[2:(n-1L)] + w[1:(n-2L)]
      bandSparse(n, n, 0:1, list(d0, rep.int(-w, n-1L)), symmetric=TRUE)
    }
  }
}

# incidence matrix for first-order random walk
D_RW1 <- function(n, w=NULL, circular=FALSE) {
  if (n < 2L) stop("RW1 model needs a sequence of at least 2 periods")
  if (is.null(w)) {
    if (circular)
      bandSparse(n, n, c(0L, 1L, 1L-n), list(rep.int(-1, n), rep.int(1, n), 1))
    else
      bandSparse(n-1L, n, 0:1, list(rep.int(-1, n-1L), rep.int(1, n-1L)))
  } else {
    w <- as.numeric(w)
    if (any(w <= 0) || length(w) != n-1L) stop("w must be a positive vector of length n-1")
    sqrt.w <- sqrt(w)
    if (circular)
      bandSparse(n, n, c(0L, 1L, 1L-n), list(-c(sqrt.w, sqrt.w[n-1L]), c(sqrt.w, sqrt.w[n-1L]), sqrt.w[n-1L]))
    else
      bandSparse(n - 1L, n, 0:1, list(-sqrt.w, sqrt.w))
  }
}

R_RW1 <- function(n, w=NULL, circular=FALSE) matrix(1, n, 1L)

Q_RW2 <- function(n) {
  if (n < 4L) stop("RW2 model needs a sequence of at least 4 periods")
  bandSparse(n, n, 0:2, list(c(1, 5, rep.int(6, n-4L), 5, 1), c(-2, rep.int(-4, n-3L), -2), rep.int(1, n-2L)), symmetric=TRUE)
}

D_RW2 <- function(n) {
  if (n < 4L) stop("RW2 model needs a sequence of at least 4 periods")
  # equal to D_RW1(n-1) %*% D_RW1(n)
  bandSparse(n - 2L, n, 0:2, list(rep.int(1, n - 2L), rep.int(-2, n - 2L), rep.int(1, n - 1L)))
}

R_RW2 <- function(n) {
  cbind(rep.int(1, n), seq_len(n))
}

Q_season <- function(n, period) {
  if (n < 2L * period) stop(paste("Seasonal component with period", period, "needs a sequence of at least", 2L * period, "periods"))
  band_list <- list()
  for (b in seq_len(period))
    band_list <- c(band_list, list(c(seq_len(period - b + 1L), rep.int(period - b + 1L, n - 2L * period + b - 1L), (period - b + 1L):1L)))
  bandSparse(n, n, 0:(period - 1L), band_list, symmetric=TRUE)
}

D_season <- function(n, period) {
  if (n < 2L * period) stop(paste("Seasonal component with period", period, "needs a sequence of at least", 2L * period, "periods"))
  band_list <- list()
  for (b in seq_len(period))
    band_list <- c(band_list, list(rep.int(1, n - period + 1L)))
  bandSparse(n - period + 1L, n, 0:(period - 1L), band_list)
}

R_season <- function(n, period) {
  CM <- matrix(0, n, period - 1L)
  for (b in seq_len(period - 1L)) {
    CM[seq.int(b, n, period), b] <- 1
    CM[seq.int(period, n, period), b] <- -1
  }
  CM
}

get_neighbours_list <- function(graph, snap=NULL, queen=NULL, poly.df) {
  if (is.null(graph) && !is.null(poly.df)) {
    warn("argument 'poly.df' is deprecated; please use argument 'graph' instead")
    graph <- poly.df
  }
  if (!inherits(graph, "nb")) {
    if (!inherits(graph, c("sf", "sfc", "SpatialPolygons"))) stop("no spatial structure has been provided")
    if (!requireNamespace("spdep", quietly=TRUE)) stop("package spdep required to construct a spatial precision matrix")
    if (is.null(snap)) snap <- sqrt(.Machine$double.eps)
    if (is.null(queen)) queen <- TRUE
    graph <- spdep::poly2nb(graph, snap=snap, queen=queen)
  }
  if (length(graph) < 2L) stop("spatial component needs at least 2 areas")
  graph
}

# CAR spatial precision matrix
Q_spatial <- function(graph=NULL, snap=sqrt(.Machine$double.eps), queen=TRUE, poly.df=NULL) {
  # transform into spatial precision matrix
  # cannot use spdep::nb2mat, as it returns a dense matrix
  nb2Q(get_neighbours_list(graph, snap, queen, poly.df))
}

# CAR spatial incidence matrix
D_spatial <- function(graph=NULL, snap=sqrt(.Machine$double.eps), queen=TRUE, poly.df=NULL) {
  nb2D(get_neighbours_list(graph, snap, queen, poly.df))
}

R_spatial <- function(graph=NULL, snap=sqrt(.Machine$double.eps), queen=TRUE, poly.df=NULL, derive.constraints=FALSE) {
  nbs <- get_neighbours_list(graph, snap, queen, poly.df)
  if (derive.constraints) warn("argument 'derive.constraints' is deprecated as it is no longer needed")
  comp <- attr(nbs, "ncomp")
  if (is.null(comp)) {
    if (!requireNamespace("spdep", quietly=TRUE)) stop("Package spdep required to construct a spatial precision matrix. Please install it.")
    # older versions of spdep may not add ncomp attribute to neighbours list(?)
    comp <- spdep::n.comp.nb(nbs)  # compute disconnected parts
  }
  if (comp[["nc"]] == 1L)
    matrix(1, length(nbs), 1L)
  else
    Ctab2dgC(fac2tabM("comp.id", comp))
}

# convert neighbourhood structure to a sparse precision matrix
nb2Q <- function(nb) {
  d <- spdep::card(nb)  # numbers of neighbours, diagonal of the precision matrix
  if (sum(d) %% 2L != 0L) warn("not a symmetric neighbourhood structure")
  i <- j <- c(seq_along(nb), rep.int(NA_integer_, sum(d) %/% 2L))  # row and column indices
  r <- length(nb)
  for (g in seq_along(nb)) {
    if (d[g] > 0L) {
      nbg <- nb[[g]]
      nbg <- nbg[nbg > g]  # upper triangular part
      l <- length(nbg)
      if (l > 0L) {
        int <- (r+1L):(r+l)
        i[int] <- g
        j[int] <- nbg
        r <- r + l
      }
    }
  }
  sparseMatrix(i=i, j=j, x=c(d, rep.int(-1, length(i) - length(d))), dims=c(length(nb), length(nb)), symmetric=TRUE)
}

# convert neighbourhood structure to sparse oriented incidence matrix
nb2D <- function(nb) {
  d <- spdep::card(nb)  # numbers of neighbours, diagonal of the precision matrix
  if (sum(d) %% 2L != 0L) warn("not a symmetric neighbourhood structure")
  nr <- sum(d) %/% 2; nc <- length(d)  # dimension of incidence matrix
  i <- rep_each(seq_len(nr), 2L)
  j <- rep.int(NA_integer_, length(i))
  r <- 1L
  for (g in seq_along(nb)) {
    if (d[g] > 0L) {
      nbg <- nb[[g]]
      nbg <- nbg[nbg > g]  # otherwise each edge is counted twice
      for (eg in nbg) {
        j[r] <- g
        j[r + 1L] <- eg
        r <- r + 2L
      }
    }
  }
  sparseMatrix(i=i, j=j, x=rep.int(c(-1, 1), nr), dims=c(nr, nc))
}

#' Derive constraints for, i.e. null space of, a precision matrix
#' 
#' @noRd
#' @param Q l x l precision matrix.
#' @param tol tolerance in regarding a small eigenvalue to be zero.
#' @returns An l x r Matrix with r the number of singular values.
derive_constraints <- function(Q, tol=sqrt(.Machine$double.eps)) {
  test <- eigen(Q)
  zeros <- which(test$values < tol)
  if (length(zeros))
    drop0(Matrix(test$vectors[, zeros, drop=FALSE]), tol=tol)
  else
    NULL
}

# create a template for incidence matrix of a GMRF with an AR1 factor with inferred phi
# info is factor info object created with get_factor_info
# assume that DA0.5 has been created with phi=0.5
# nr: number of the AR1 factor
DA_AR1_template <- function(info, DA0.5, nr) {
  info$factors[[nr]]$phi <- 0.25
  DA0.25 <- compute_GMRF_matrices(info, D=TRUE, Q=FALSE, R=FALSE)$D
  # indices for -phi:
  ind1 <- which(abs(DA0.25@x - 0.5 * DA0.5@x) < sqrt(.Machine$double.eps))
  # indices for 1 + phi^2:
  ind2 <- which(abs(DA0.25@x - (sqrt(1 - 0.25^2) / sqrt(1 - 0.5^2)) * DA0.5@x) < sqrt(.Machine$double.eps))
  rm(DA0.25, info)
  update <- function(phi) {
    DA <- DA0.5
    attr(DA, "x")[ind1] <- 2 * phi * DA@x[ind1]
    attr(DA, "x")[ind2] <- 2 * sqrt((1 - phi^2)/3) * DA@x[ind2]
    DA
  }
  environment()
}

# create a template for precision matrix of a GMRF with an AR1 factor with inferred phi
# info is factor info object created with get_factor_info
# assume that QA0.5 has been created with phi=0.5
# nr: number of the AR1 factor
QA_AR1_template <- function(info, QA0.5, nr) {
  info$factors[[nr]]$phi <- 0.25
  QA0.25 <- compute_GMRF_matrices(info, D=FALSE, Q=TRUE, R=FALSE, sparse=TRUE)$Q
  # indices for -phi:
  ind1 <- which(abs(QA0.25@x - 0.5 * QA0.5@x) < sqrt(.Machine$double.eps))
  # indices for 1 + phi^2:
  ind2 <- which(abs(QA0.25@x - ((1 + 0.25^2) / (1 + 0.5^2)) * QA0.5@x) < sqrt(.Machine$double.eps))
  rm(QA0.25, info)
  update <- function(phi) {
    QA <- QA0.5
    attr(QA, "x")[ind1] <- 2 * phi * QA@x[ind1]
    attr(QA, "x")[ind2] <- 0.8 * (1 + phi^2) * QA@x[ind2]
    QA
  }
  environment()
}


#' Maximize the log-likelihood or log-posterior as defined by a sampler closure
#'
#' @examples
#' \donttest{
#' n <- 1000
#' dat <- data.frame(
#'   x = rnorm(n),
#'   f = factor(sample(1:50, n, replace=TRUE))
#' )
#' df <- generate_data(
#'   ~ reg(~x, name="beta", prior=pr_normal(precision=1)) + gen(~x, factor=~f, name="v"),
#'   sigma.fixed=TRUE, data=dat
#' )
#' dat$y <- df$y
#' sampler <- create_sampler(y ~ x + gen(~x, factor=~f, name="v"), data=dat)
#' opt <- maximize_log_lh_p(sampler)
#' str(opt)
#' plot(df$par$v, opt$par$v); abline(0, 1, col="red")
#' }
#'
#' @export
#' @param sampler sampler function closure, i.e. the return value of a call to \code{\link{create_sampler}}.
#' @param type either "llh" (default) or "lpost", for optimization of the log-likelihood,
#'  or the log-posterior, respectively.
#' @param method optimization method, passed to \code{\link[stats]{optim}}.
#' @param control control parameters, passed to \code{\link[stats]{optim}}.
#' @param ... other parameters passed to \code{\link[stats]{optim}}.
#' @returns A list of parameter values that, provided the optimization was successful, maximize the (log-)likelihood
#'  or (log-)posterior.
maximize_log_lh_p <- function(sampler, type=c("llh", "lpost"), method="BFGS", control=list(fnscale=-1), ...) {
  type <- match.arg(type)
  e <- sampler
  if (e$has.bart) stop("optimization not supported for models with a 'brt' component")
  if (!is.null(e$Vmod)) stop("optimization not yet implemented for models with sampling variance model specified in 'formula.V'")
  ind_sigma <- e$vec_list[["sigma_"]]
  if (type == "llh") {
    f <- function(x)
      if (!is.null(ind_sigma) && x[ind_sigma] < 0) -Inf else e$llh_opt(x)
  } else {
    logposterior_opt <- make_logposterior_opt(e)
    f <- function(x)
      if (!is.null(ind_sigma) && x[ind_sigma] < 0) -Inf else logposterior_opt(x)
  }
  res <- optim(par=e$list2vec(sampler$start()), f, method=method, control=control, ...)
  res$par <- e$vec2list(res$par)
  res
}

make_logposterior_opt <- function(sampler) {
  # logprior for optimization of logposterior = logprior + llh; llh is defined in sampler
  log_prior <- function(p) {out <- 0}
  for (k in seq_along(sampler$mod)) {
    mc <- sampler$mod[[k]]
    switch(mc[["type"]],
      reg = {
        mc$prior$setup_logprior(mc$name)
        log_prior <- add(log_prior, bquote(out <- out + mod[[.(k)]]$prior$logprior(p)))
      },
      stop("TBI: log-prior for components other than 'reg'")
    )
  }
  # TODO add log-priors of Vmod components
  log_prior <- add(log_prior, quote("out"))
  assign("log_prior", log_prior, environment(sampler))
  # log-posterior function for optimization: function of a vector instead of list
  f <- function(x) {}
  f <- f |>
    add(quote(p <- vec2list(x))) |>
    add(quote(p$e_ <- compute_e(p))) |>
    add(quote(log_prior(p) + llh(p)))
  environment(f) <- environment(sampler)
  f
}
