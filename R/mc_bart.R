#' Create a model component object for a BART (Bayesian Additive Regression Trees)
#' component in the linear predictor
#'
#' This function is intended to be used on the right hand side of the
#' \code{formula} argument to \code{\link{create_sampler}} or
#' \code{\link{generate_data}}. It creates a BART term in the
#' model's linear predictor. To use this model component one needs
#' to have R package \pkg{dbarts} installed.
#'
#' @examples
#' # generate data, based on an example in Friedman (1991)
#' gendat <- function(n=200L, p=10L, sigma=1) {
#'   x <- matrix(runif(n * p), n, p)
#'   mu <- 10*sin(pi*x[, 1] * x[, 2]) + 20*(x[, 3] - 0.5)^2 + 10*x[, 4] + 5*x[, 5]
#'   y <- mu + sigma * rnorm(n)
#'   data.frame(x=x, mu=mu, y=y)
#' }
#'
#' train <- gendat()
#' test <- gendat(n=25)
#'
#' # keep trees for later prediction based on new data
#' sampler <- create_sampler(
#'   y ~ brt(~ . - y, name="bart", keepTrees=TRUE),
#'   family = f_gaussian(var.prior=pr_invchisq(df=3, scale=var(train$y))),
#'   data = train
#' )
#' # increase burnin and n.iter below to improve MCMC convergence
#' sim <- MCMCsim(sampler, n.chain=2, burnin=200, n.iter=500, thin=2,
#'   store.all=TRUE, verbose=FALSE)
#' (summ <- summary(sim))
#' plot(train$mu, summ$bart[, "Mean"]); abline(0, 1)
#' # NB prediction is currently slow
#' \donttest{
#' pred <- predict(sim, newdata=test,
#'   iters=sample(seq_len(n_draws(sim)), 100),
#'   show.progress=FALSE
#' )
#' (summpred <- summary(pred))
#' plot(test$mu, summpred[, "Mean"]); abline(0, 1)
#' }
#'
#' @export
#' @param formula a formula specifying the predictors to be used in the BART
#'  model component. Variable names are looked up in the data frame
#'  passed as \code{data} argument to \code{\link{create_sampler}} or
#'  \code{\link{generate_data}}, or in \code{environment(formula)}.
#' @param X a design matrix can be specified directly, as an alternative
#'  to the creation of one based on \code{formula}. If \code{X} is
#'  specified \code{formula} is ignored.
#' @param n.trees number of trees used in the BART ensemble.
#' @param name the name of the model component. This name is used in the output of the
#'  MCMC simulation function \code{\link{MCMCsim}}. By default the name will be 'bart'
#'  with the number of the model term attached.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @param keepTrees whether to store the trees ensemble for each Monte Carlo draw. This
#'  is required for prediction based on new data. The default is \code{FALSE} to save
#'  memory.
#' @param ... parameters passed to \code{\link[dbarts]{dbarts}}.
#' @returns An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Intended
#'  for internal use by other package functions.
#' @references
#'  H.A. Chipman, E.I. Georgea and R.E. McCulloch (2010).
#'    BART: Bayesian additive regression trees.
#'    The Annals of Applied Statistics 4(1), 266-298.
#'    
#'  J.H. Friedman (1991).
#'    Multivariate adaptive regression splines.
#'    The Annals of Statistics 19, 1-67.
brt <- function(formula, X=NULL, n.trees=75L,
                name="", debug=FALSE, keepTrees=FALSE, ...) {

  e <- sys.frame(-2L)
  type <- "brt"
  if (name == "") stop("missing model component name")

  if (!requireNamespace("dbarts", quietly=TRUE)) stop("package dbarts required for a model including Bayesian Additive Regression Trees")

  if (e[["modeled.Q"]] && e[["Q0.type"]] == "symm") stop("BART component not compatible with non-diagonal residual variance matrix")

  if (is.null(X)) {
    # no contrasts applied, no (urgent) need to remove redundancy
    X <- model_matrix(formula, e[["data"]], contrasts.arg="contr.none", sparse=FALSE)
  } else {
    if (is.null(dimnames(X)[[2L]])) colnames(X) <- seq_len(ncol(X))
  }
  X <- economizeMatrix(X, sparse=FALSE, strip.names=FALSE, check=TRUE)
  q <- ncol(X)

  control <- dbarts::dbartsControl(n.chains=1L, updateState=FALSE, n.trees=n.trees)

  in_block <- FALSE

  name_sampler <- paste0(name, "_sampler_")  # trailing '_' --> not stored by MCMCsim even if store.all=TRUE
  lp <- function(p) copy_vector(p[[name]])
  lp_update <- function(x, plus=TRUE, p) v_update(x, plus, p[[name]])
  draws_linpred <- function(obj, units=NULL, chains=NULL, draws=NULL, matrix=FALSE) {
    if (is.null(obj[[name]])) stop("no fitted values for 'brt' component; please re-run MCMCsim with 'store.all=TRUE'")
    if (matrix)
      out <- as.matrix.dc(get_from(obj[[name]], vars=units, chains=chains, draws=draws))
    else {
      out <- list()
      for (ch in seq_along(chains))
        out[[ch]] <- get_from(obj[[name]], vars=units, chains=chains[ch], draws=draws)[[1L]]
    }
    out
  }

  # TODO if sigma.fixed set prior to (almost) fix sigma to 1
  create_BART_sampler <- function() {
    if (e[["Q0.type"]] == "diag") {
      dbarts::dbarts(formula=X, data=e$y_eff(),
        weights = e$Q0@x, sigma = if (e[["sigma.fixed"]]) 1 else NA_real_,
        control=control, ...
      )
    } else {
      dbarts::dbarts(formula=X, data=e$y_eff(),
        sigma = if (e[["sigma.fixed"]]) 1 else NA_real_,
        control=control, ...
      )
    }
  }

  start <- function(p) {
    # TODO multinomial family
    if (is.null(p[[name_sampler]])) {
      p[[name_sampler]] <- create_BART_sampler()
    }
    p[[name]] <- p[[name_sampler]]$run(0L, 1L, FALSE)$train[, 1L]
    p
  }

  if (e$family$family == "multinomial") {
    edat <- new.env(parent = environment(formula))
    environment(formula) <- edat
  }

  if (keepTrees) name_trees <- paste0(name, "_", "trees_")
  make_predict <- function(newdata=NULL, Xnew=NULL, verbose=TRUE) {
    if (is.null(newdata) && is.null(Xnew)) {
      # in-sample prediction
      linpred <- function(p) copy_vector(p[[name]])
      linpred_update <- function(x, plus=TRUE, p) v_update(x, plus, p[[name]])
    } else {
      if (is.null(newdata)) {
        if (missing(Xnew)) stop("one of 'newdata' and 'Xnew' should be supplied")
      } else {
        if (!keepTrees) stop("out-of-sample prediction requires setting keepTrees=TRUE in brt() model specification")
        nnew <- nrow(newdata)
        if (e$family[["family"]] == "multinomial") {
          Xnew <- NULL
          for (k in seq_len(e[["Km1"]])) {
            edat$cat_ <- factor(rep.int(e$cats[k], nnew), levels=e$cats[-length(e$cats)])
            Xnew <- rbind(Xnew, model_matrix(formula, data=newdata, contrasts.arg="contr.none", sparse=FALSE))
          }
          edat$cat_ <- NULL
        } else {
          Xnew <- model_matrix(formula, newdata, contrasts.arg="contr.none", sparse=FALSE)
        }
      }
      if (is.null(dimnames(Xnew)[[2L]])) {
        if (ncol(Xnew) != q) stop("wrong number ", ncol(Xnew), " predictor column(s) for model term '", name, "' versus ", q, " originally")
      } else {
        Xnew <- Xnew[, dimnames(X)[[2L]], drop=FALSE]
      }
      Xnew <- economizeMatrix(Xnew, sparse=FALSE, strip.names=TRUE, check=TRUE)
      rm(newdata, verbose)
      nnew <- nrow(Xnew)

      # predict based on Xnew using stored tree; see dbarts' 'Working with saved trees' vignette
      # TODO rewrite in C++
      tree_predict <- function(p) {
        getPredictionsForTreeRecursive <- function(tree, indices) {
          if (tree$var[1L] == -1L) {
            predictions[indices] <<- predictions[indices] + tree$value[1L]
            return(1L)
          }
          goesLeft <- Xnew[indices, tree$var[1L]] <= tree$value[1L]
          headOfLeftBranch <- tree[-1L, ]
          n_nodes.left <- getPredictionsForTreeRecursive(headOfLeftBranch, indices[goesLeft])
          headOfRightBranch <- tree[seq.int(2L + n_nodes.left, dim(tree)[1L]), ]
          n_nodes.right <- getPredictionsForTreeRecursive(headOfRightBranch, indices[!goesLeft])
          return(1L + n_nodes.left + n_nodes.right)
        }
        trees <- p[[name_trees]]
        predictions <- numeric(nnew)
        for (t in seq_len(n.trees))
          getPredictionsForTreeRecursive(trees[trees$tree == t, ], seq_len(nnew))
        predictions
      }
      linpred <- function(p) tree_predict(p)
      linpred_update <- function(x, plus=TRUE, p) v_update(x, plus, tree_predict(p))
    }
    environment()
  }

  rprior <- function(p) {
    if (is.null(p[[name_sampler]])) {
      p[[name_sampler]] <- create_BART_sampler()
    }
    p[[name_sampler]]$sampleTreesFromPrior()
    p[[name_sampler]]$sampleNodeParametersFromPrior()
    p[[name]] <- p[[name_sampler]]$predict(X)
    if (keepTrees) {
      # getTrees returns a data.frame; 2nd column 'n' is not needed for prediction
      p[[name_trees]] <- p[[name_sampler]]$getTrees()[, -2L]
    }
    p
  }

  if (!e[["prior.only"]]) {
    draw <- if (debug) function(p) {browser()} else function(p) {}
    if (!e[["single.block"]]) {
      if (e[["e.is.res"]])
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] + p[[.(name)]]))
      else
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] - p[[.(name)]]))
      draw <- add(draw, bquote(p[[.(name_sampler)]]$setResponse(p[["e_"]])))
    }
    if (e[["modeled.Q"]])
      draw <- add(draw, bquote(p[[.(name_sampler)]]$setWeights(p[["Q_"]])))
    draw <- draw |>
      add(bquote(p[[.(name_sampler)]]$setSigma(.(if (e$sigma.fixed) 1 else quote(p[["sigma_"]]))))) |>
      add(bquote(p[[.(name)]] <- p[[.(name_sampler)]]$run(0L, 1L)$train[, 1L]))
    if (keepTrees) {
      # getTrees returns a data.frame; 2nd column 'n' is not needed for prediction
      draw <- add(draw, bquote(p[[.(name_trees)]] <- p[[.(name_sampler)]]$getTrees()[, -2L]))
      # undo internal scaling of response vector
      # dbarts scales the 'data' to lie within -0.5 and 0.5
      # define an inverse transformation for prediction purposes
      # min.y and range.y are set below, if setResponse is used should be updated in each iteration
      # undo_scaling <- function(pred) min.y + range.y * (pred + 0.5)
      if (e[["single.block"]]) {
        min.y <- min(e$y_eff())
        range.y <- max(e$y_eff()) - min.y
      } else {
        draw <- add(draw, quote(min.y <- min(p[["e_"]])))
        draw <- add(draw, quote(range.y <- max(p[["e_"]]) - min.y))
      }
      draw <- add(draw, bquote(p[[.(name_trees)]]$value[p[[.(name_trees)]]$var == -1L] <- min.y/n.trees + range.y * (p[[.(name_trees)]]$value[p[[.(name_trees)]]$var == -1L] + 0.5/n.trees)))
    }
    if (e[["single.block"]]) {
      if (e[["e.is.res"]])
        draw <- add(draw, bquote(p$e_ <- e$y_eff() - p[[.(name)]]))
      else
        draw <- add(draw, bquote(p$e_ <- p[[.(name)]]))
    } else {
      if (e[["e.is.res"]])
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] - p[[.(name)]]))
      else
        draw <- add(draw, bquote(p$e_ <- p[["e_"]] + p[[.(name)]]))
    }
    draw <- add(draw, quote(p))
  }

  environment()
}
