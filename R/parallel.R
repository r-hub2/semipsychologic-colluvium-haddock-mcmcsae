#' Set up a cluster for parallel computing
#'
#' The cluster is set up for a number of workers by loading the \pkg{mcmcsae}
#' package and setting up independent RNG streams.
#'
#' @export
#' @param n.cores the number of cpu cores to use.
#' @param seed optional random seed for reproducibility.
#' @param export a character vector with names of objects to export to the workers.
#' @returns An object representing the cluster.
setup_cluster <- function(n.cores=NULL, seed=NULL, export=NULL) {
  if (is.null(n.cores)) {
    n.cores <- max(1L, parallel::detectCores() - 1L)
  } else {
    n.cores <- min(n.cores, parallel::detectCores())
  }
  message("setting up cluster on ", n.cores, " cores")
  cl <- parallel::makeCluster(n.cores)
  # make sure the same libPaths are set, so that packages are found
  parallel::clusterCall(cl, function(x) .libPaths(x), .libPaths())
  parallel::clusterEvalQ(cl, library(mcmcsae))
  # set up independent RNG streams, selecting L'Ecuyer-CMRG RNG
  if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
  if (!is.null(export)) parallel::clusterExport(cl, export)
  cl
}

#' Stop a cluster
#'
#' Stop a cluster set up by \code{\link{setup_cluster}}.
#'
#' @export
#' @param cl the cluster object.
#' @returns \code{NULL}.
stop_cluster <- function(cl) {
  parallel::stopCluster(cl)
}

#' Combine multiple mcdraws objects into a single one by combining their chains
#'
#' This function can be used to combine the results of parallel simulations.
#'
#' @export
#' @param ... objects of class \code{mcdraws}.
#' @returns A combined object of class \code{mcdraws} where the number of stored
#'  chains equals the sum of the numbers of chains in the input objects.
combine_chains <- function(...) {
  dotargs <- list(...)
  out <- list()
  out[["_info"]] <- dotargs[[1L]][["_info"]]
  out[["_info"]]$n.chain <- sum(i_apply(dotargs, n_chains))
  out[["_model"]] <- dotargs[[1L]][["_model"]]
  out[["_state"]] <- do.call("c", lapply(dotargs, `[[`, "_state"))
  if (!is.null(dotargs[[1L]][["_accept"]])) {
    out[["_accept"]] <- list()
    for (v in names(dotargs[[1L]][["_accept"]]))
      out[["_accept"]][[v]] <- do.call("c", lapply(dotargs, function(x) x[["_accept"]][[v]]))
  }
  if (!is.null(dotargs[[1L]][["_means"]])) {
    out[["_means"]] <- list()
    for (v in names(dotargs[[1L]][["_means"]]))
      out[["_means"]][[v]] <- do.call("c", lapply(dotargs, function(x) x[["_means"]][[v]]))
  }
  if (!is.null(dotargs[[1L]][["_sds"]])) {
    out[["_sds"]] <- list()
    for (v in names(dotargs[[1L]][["_sds"]]))
      out[["_sds"]][[v]] <- do.call("c", lapply(dotargs, function(x) x[["_sds"]][[v]]))
  }
  for (v in par_names(dotargs[[1L]])) {
    out[[v]] <- do.call("c", lapply(dotargs, `[[`, v))
    if (!is.null(attr(dotargs[[1L]][[v]], "labels"))) attr(out[[v]], "labels") <- attr(dotargs[[1L]][[v]], "labels")
    class(out[[v]]) <- "dc"
  }
  if (!is.null(dotargs[[1L]][["_info"]][["list.pars"]])) {
    for (v in dotargs[[1L]][["_info"]][["list.pars"]])
      out[[v]] <- do.call("c", lapply(dotargs, `[[`, v))
  }
  class(out) <- "mcdraws"
  out
}

# obj is a list of dc objects
combine_chains_dc <- function(obj) {
  out <- unlst(obj, recursive=FALSE)
  if (!is.null(attr(obj[[1L]], "ppp"))) {
    l <- length(attr(obj[[1L]], "ppp"))
    if (l == 1L)
      attr(out, "ppp") <- mean(vapply(obj, attr, 0, "ppp"))
    else
      attr(out, "ppp") <- rowMeans(vapply(obj, attr, numeric(l), "ppp"))
  }
  if (!is.null(attr(obj[[1L]], "labels"))) attr(out, "labels") <- attr(obj[[1L]], "labels")
  class(out) <- "dc"
  out
}

# obj can be of class mcdraws or dc (the latter may also be a list)
split_chains <- function(obj, parts=NULL) {
  n.chain <- n_chains(obj)
  if (is.null(parts)) {
    parts <- n.chain  # by default, split into single chains
  } else {
    parts <- as.integer(parts)
  }
  if (parts > n.chain || parts < 1L) stop("wrong input for 'parts'")
  chs <- rep(n.chain %/% parts, parts) + rep(1:0, c(n.chain %% parts, parts - n.chain %% parts))
  chs <- split(seq_len(n.chain), rep(seq_len(parts), chs))
  out <- vector(mode="list", parts)
  obj_class <- class(obj)[1L]
  if (all(obj_class != c("mcdraws", "dc", "list"))) stop("unsupported class '", obj_class, "'")
  for (k in seq_along(out)) {
    if (obj_class == "mcdraws") {
      out[[k]][["_info"]] <- obj[["_info"]]
      out[[k]][["_info"]]$n.chain <- length(chs[[k]])
      out[[k]][["_model"]] <- obj[["_model"]]
      out[[k]][["_state"]] <- obj[["_state"]][chs[[k]]]
      if (!is.null(obj[["_accept"]])) {
        out[[k]][["_accept"]] <- list()
        for (v in names(obj[["_accept"]]))
          out[[k]][["_accept"]][[v]] <- obj[["_accept"]][[v]][chs[[k]]]
      }
      if (!is.null(obj[["_means"]])) {
        out[[k]][["_means"]] <- list()
        for (v in names(obj[["_means"]]))
          out[[k]][["_means"]][[v]] <- obj[["_means"]][[v]][chs[[k]]]
      }
      if (!is.null(obj[["_sds"]])) {
        out[[k]][["_sds"]] <- list()
        for (v in names(obj[["_sds"]]))
          out[[k]][["_sds"]][[v]] <- obj[["_sds"]][[v]][chs[[k]]]
      }
      for (v in par_names(obj))
        out[[k]][[v]] <- obj[[v]][chs[[k]]]
    } else {
      out[[k]] <- obj[chs[[k]]]
    }
    class(out[[k]]) <- obj_class
  }
  out
}

#' Combine multiple mcdraws objects into a single one by combining their draws
#'
#' This function is used to combine the results of parallel posterior predictive simulations.
#'
#' @export
#' @param ... objects of class \code{mcdraws}
#' @returns A combined object of class \code{mcdraws} where the number of stored
#'  draws equals the sum of the numbers of draws in the input objects.
combine_iters <- function(...) {
  dotargs <- list(...)
  out <- list()
  if (!is.null(dotargs[[1L]][["_info"]])) {
    out[["_info"]] <- dotargs[[1L]][["_info"]]
    out[["_info"]]$n.draw <- sum(i_apply(dotargs, n_draws))
  }
  out[["_model"]] <- dotargs[[1L]][["_model"]]
  out[["_state"]] <- dotargs[[length(dotargs)]][["_state"]]
  if (!is.null(dotargs[[1L]][["_means"]])) {
    out[["_means"]] <- list()
    for (v in names(dotargs[[1L]][["_means"]]))
      out[["_means"]][[v]] <- do.call("c", lapply(dotargs, function(x) x[["_means"]][[v]]))
  }
  if (!is.null(dotargs[[1L]][["_sds"]])) {
    out[["_sds"]] <- list()
    for (v in names(dotargs[[1L]][["_sds"]]))
      out[["_sds"]][[v]] <- do.call("c", lapply(dotargs, function(x) x[["_sds"]][[v]]))
  }
  for (v in par_names(dotargs[[1L]])) {
    out[[v]] <- dotargs[[1L]][[v]]
    for (obj.i in dotargs[-1L]) {
      for (ch in seq_along(out[[v]]))
        out[[v]][[ch]] <- rbind(out[[v]][[ch]], obj.i[[v]][[ch]])
    }
    class(out[[v]]) <- "dc"
  }
  class(out) <- "mcdraws"
  out
}

# obj is a list of dc objects
combine_iters_dc <- function(obj) {
  chains <- seq_along(obj[[1L]])
  n.draw <- sum(i_apply(obj, function(x) nrow(x[[1L]])))
  n.var <- ncol(obj[[1L]][[1L]])
  out <- list()
  if (is.integer(obj[[1L]][[1L]][1L, 1L]))
    for (ch in chains) out[[ch]] <- matrix(NA_integer_, nrow=n.draw, ncol=n.var)
  else
    for (ch in chains) out[[ch]] <- matrix(NA_real_, nrow=n.draw, ncol=n.var)
  r <- 1L
  for (obj.i in obj) {
    ind <- r:(r - 1L + nrow(obj.i[[1L]]))
    for (ch in chains) out[[ch]][ind, ] <- obj.i[[ch]]
    r <- ind[length(ind)] + 1L
  }
  if (!is.null(attr(obj[[1L]], "ppp"))) {
    if (length(attr(obj[[1L]], "ppp")) == 1L)
      attr(out, "ppp") <- sum(vapply(obj, attr, 0, "ppp")) / (length(chains) * n.draw)
    else
      attr(out, "ppp") <- rowSums(vapply(obj, attr, numeric(n.var), "ppp")) / (length(chains) * n.draw)
  }
  if (!is.null(attr(obj[[1L]], "labels"))) attr(out, "labels") <- attr(obj[[1L]], "labels")
  class(out) <- "dc"
  out
}

# obj can be of class dc or mcdraws
# parts: in how many parts
# iters: subset of iterations; default is all iterations in obj
split_iters <- function(obj, iters=NULL, parts=NULL) {
  parts <- as.integer(parts)
  if (length(parts) != 1L || parts < 1L) stop("wrong input for 'parts'")
  if (is.null(iters))
    iters <- seq_len(n_draws(obj))
  else
    if (!all(iters %in% seq_len(n_draws(obj)))) stop("non-existing iterations selected")
  n.iter <- length(iters)
  if (n.iter < parts) stop("cannot split ", n.iter, " draws in ", parts, " parts")
  its <- rep(n.iter %/% parts, parts) + rep(1:0, c(n.iter %% parts, parts - n.iter %% parts))
  its <- split(iters, rep(seq_len(parts), its))
  out <- vector(mode="list", parts)
  obj_class <- class(obj)[1L]
  if (all(obj_class != c("dc", "mcdraws"))) stop("unsupported class '", obj_class, "'")
  for (k in seq_along(out)) {
    if (obj_class == "mcdraws") {
      out[[k]][["_info"]] <- obj[["_info"]]
      out[[k]][["_info"]]$n.draw <- length(its[[k]])
      out[[k]][["_model"]] <- obj[["_model"]]
      for (v in par_names(obj))
        out[[k]][[v]] <- get_from(obj[[v]], draws=its[[k]])
      if (!is.null(obj[["_info"]][["list.pars"]]))
        for (v in obj[["_info"]][["list.pars"]]) out[[k]][[v]] <- lapply(obj[[v]], function(x) x[its[[k]]])
    } else {
      out[[k]] <- get_from(obj, draws=its[[k]])
    }
    class(out[[k]]) <- obj_class
  }
  out
}
