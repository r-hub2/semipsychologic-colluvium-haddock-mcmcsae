
#' Simulation based calibration
#'
#' @examples
#' \dontrun{
#' # this example may take a long time
#' n <- 10L
#' dat <- data.frame(x=runif(n))
#' ranks <- SBC_test(~ reg(~ 1 + x, prior=pr_normal(mean=c(0.25, 1), precision=1), name="beta"),
#'   sigma.mod=pr_invchisq(df=1, scale=list(df=1, scale=1)), data=dat,
#'   pars=list(mu="beta[1]", beta_x="beta[2]", sigma="sigma_"),
#'   n.draw=9L, n.sim=10L*20L, thin=2L, burnin=20L
#' )
#' ranks
#' }
#'
#' @export
#' @param ... passed to \code{\link{create_sampler}} (can be all parameters except \code{prior.only})
#' @param pars named list with univariate functions of the parameters to use in test. This list
#'   is passed to argument \code{pred} of \code{\link{MCMCsim}}.
#' @param n.draw number of posterior draws to retain in posterior simulations.
#' @param n.sim number of simulation iterations.
#' @param burnin burnin to use in posterior simulations, passed to \code{\link{MCMCsim}}.
#' @param thin thinning to use in posterior simulations, passed to \code{\link{MCMCsim}}.
#' @param show.progress whether a progress bar should be shown.
#' @param verbose set to \code{FALSE} to suppress messages.
#' @param n.cores the number of cpu cores to use. Default is one, i.e. no parallel computation.
#'  If an existing cluster \code{cl} is provided, \code{n.cores} will be set to the number
#'  of workers in that cluster.
#' @param cl an existing cluster can be passed for parallel computation. If \code{NULL} and
#'  \code{n.cores > 1}, a new cluster is created.
#' @param seed a random seed (integer). For parallel computation it is used to independently
#'  seed RNG streams for all workers.
#' @param export a character vector with names of objects to export to the workers. This may
#'  be needed for parallel execution if expressions in the model formulae depend on global variables.
#' @returns A matrix with ranks.
#' @references
#'   M. Modrak, A.H. Moon, S. Kim, P. Buerkner, N. Huurre, K. Faltejskova,
#'     A. Gelman and A. Vehtari (2023).
#'     Simulation-based calibration checking for Bayesian computation:
#'     The choice of test quantities shapes sensitivity.
#'     Bayesian Analysis, 1(1), 1-28.
SBC_test <- function(..., pars, n.draw=25L, n.sim=20L*n.draw, burnin=25L, thin=2L,
                      show.progress=TRUE, verbose=TRUE,
                      n.cores=1L, cl=NULL, seed=NULL, export=NULL) {
  sampler_args <- list(...)
  sampler_args$prior.only <- TRUE
  sampler <- do.call("create_sampler", sampler_args)
  # TODO: in some cases like TMVN prior draws may not be iid; in that case use burnin and thin here as well
  if (names(sampler_args)[1L] == "") names(sampler_args)[1L] <- "formula"
  sampler_args$formula <- update.formula(sampler_args$formula, y.tilde ~ .)

  n.pars <- length(pars)

  run_sbc <- function(n.sim) {
    ranks <- matrix(0, n.draw + 1L, n.pars, dimnames=list(NULL, names(pars)))
    for (i in seq_len(n.sim)) {
      # 1. draw from prior
      sim0 <- MCMCsim(sampler, burnin=0L, n.iter=1L, n.chain=1L, pred=pars, from.prior=TRUE, verbose=FALSE)
      # 2. generate data using the prior draw
      pred <- predict(sim0, show.progress=FALSE)
      y.tilde <- pred[[1L]][1L, ]
      # 3. run the posterior sampler based on the generated data
      sampler_args$prior.only <- FALSE
      environment(sampler_args$formula) <- environment()
      sampler <- do.call("create_sampler", sampler_args)
      sim <- MCMCsim(sampler, pred=pars, burnin=burnin, n.chain=1L, n.iter=thin*n.draw, thin=thin, verbose=FALSE)
      # 4. compute ranks of prior draws in the sets of posterior draws
      for (v in seq_along(pars)) {
        vname <- names(pars)[v]
        r <- sum(sim[[vname]][[1L]] < sim0[[vname]][[1L]][1L, ]) + 1L
        ranks[r, v] <- ranks[r, v] + 1L
      }
      if (show.progress) setTxtProgressBar(pb, i)
    }
    ranks
  }

  if (is.null(cl))
    n.cores <- max(as.integer(n.cores)[1L], 1L)
  else
    n.cores <- length(cl)
  if (n.cores > 1L) {
    show.progress <- FALSE
    if (n.cores > parallel::detectCores()) stop("too many cores")
    if (is.null(cl)) {
      cl <- setup_cluster(n.cores, seed)
      on.exit(stop_cluster(cl))
    } else {
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    }
    if (!is.null(export)) {
      if (!is.character(export)) stop("'export' must be a character vector")
      parallel::clusterExport(cl, export, parent.frame())
    }
    if (verbose) message(n.sim, " model fits distributed over ", n.cores, " cores")
    sim_list <- rep.int(n.sim %/% n.cores, n.cores) + c(rep.int(1L, n.sim %% n.cores), integer(n.cores - n.sim %% n.cores))
    out <- Reduce(`+`, parallel::parLapply(cl, X=sim_list, fun=run_sbc))
  } else {
    if (!is.null(seed)) set.seed(seed)
    environment(sampler_args$formula) <- environment()
    if (show.progress) pb <- txtProgressBar(min=1L, max=n.sim, style=3L)
    out <- run_sbc(n.sim)
    if (show.progress) close(pb)
  }
  out
}
