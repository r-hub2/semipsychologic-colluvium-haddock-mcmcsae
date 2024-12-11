#' Create a model component object for a variance factor component in the variance function of a
#' gaussian sampling distribution
#'
#' This function is intended to be used on the right hand side of the \code{formula.V} argument to
#' \code{\link{create_sampler}} or \code{\link{generate_data}}.
#'
#' @export
#' @param factor The name of a factor variable. The name \code{"local_"} has a special meaning,
#'  and assigns a different variance scale parameter to each data unit.
#'  In case of inverse chi-squared priors this implies that the marginal sampling distribution
#'  is a t distribution. In case of exponential priors the marginal sampling distribution
#'  is a Laplace or double exponential distribution.
#' @param prior the prior assigned to the variance factors. Currently the prior can be inverse chi-squared
#'  or exponential, specified by a call to \code{\link{pr_invchisq}} or \code{\link{pr_exp}}, respectively.
#'  The default priors are inverse chi-squared with 1 degree of freedom. See the help pages of the
#'  prior specification functions for details on how to set non-default priors.
#' @param name The name of the variance model component. This name is used in the output of the MCMC simulation
#'  function \code{\link{MCMCsim}}. By default the name will be 'vfac' with the number of the variance model
#'  term attached.
#' @param debug If \code{TRUE} a breakpoint is set at the beginning of the posterior
#'  draw function associated with this model component. Mainly intended for developers.
#' @returns An object with precomputed quantities and functions for sampling from
#'  prior or conditional posterior distributions for this model component. Intended
#'  for internal use by other package functions.
# TODO generalize inverse chi-squared (including beta-prime) and exponential priors to generalized hyperbolic
vfac <- function(factor="local_",
                 prior=pr_invchisq(df=1, scale=1),
                 name="", debug=FALSE) {

  e <- sys.frame(-2L)
  type <- "vfac"
  if (name == "") stop("missing model component name")

  switch(factor,
    local_ = {
      X <- CdiagU(e[["n"]])
      nh <- 1
    },
    global_ = {
      X <- aggrMatrix(rep.int(1L, e[["n"]]))
      nh <- e[["n"]]
    },
    {
      X <- aggrMatrix(e$data[[factor]])
      nh <- as.numeric(colSums(X))
      e$coef.names[[name]] <- levels(e$data[[factor]])
    }
  )
  q <- ncol(X)

  if (e$Q0.type == "symm") {
    # check that scale factor is compatible with block-diagonal structure of Q0
    if (sum(abs(commutator(e[["Q0"]], Cdiag(X %m*v% sin(1.23*seq_len(q)))))) > sqrt(.Machine$double.eps)) stop("'formula.V' incompatible with 'Q0'")
  }

  switch(prior$type,
    invchisq = prior$init(q, !e$prior.only),
    exp = prior$init(q, !e$prior.only)
  )

  compute_Qfactor <- function(p) X %m*v% (1 / p[[name]])

  # function that creates an oos prediction function (closure) based on new data
  make_predict_Vfactor <- function(newdata) {
    if (factor == "local_") {
      # "local_" corresponds to unit-level, e.g. Student-t sampling distribution
      qnew <- nrow(newdata)
    } else {
      # TODO allow random generation for oos-levels more generally (currently only for Student-t sampling)
      Xnew <- aggrMatrix(newdata[[factor]])
    }
    rm(newdata)
    pred <- function(p) {}
    if (factor == "local_") {
      switch(prior$type,
        invchisq = {
          # NB this assumes that all hyperparameters are scalar! TODO check this
          # Student t marginal sampling distribution, provided scale is fixed
          if (is.list(prior$df))
            pred <- add(pred, bquote(df <- p[[.(name_df)]]))
          else
            pred <- add(pred, quote(df <- prior$df))
          if (is.list(prior$scale)) {
            if (prior$scale$common) {
              # for this case we need to store the common scale parameter!
              stop("TBI: oos prediction for common scale")
            } else {
              pred <- add(pred, bquote(draw_betaprime(.(qnew), 0.5*prior$scale$df, 0.5*df, df/prior$psi0)))
            }
          } else {
            pred <- add(pred, bquote(1 / rchisq_scaled(.(qnew), df, prior$scale)))
          }
        },
        exp = {  # Laplace marginal sampling distribution
          pred <- add(pred, bquote(rexp(.(qnew))))
        }
      )
    } else {
      pred <- add(pred, bquote(Xnew %m*v% p[[.(name)]]))
    }
    pred
  }

  rprior <- function(p) {}
  switch(prior$type,
    invchisq = {
      if (is.list(prior$df)) {
        name_df <- paste0(name, "_df")
        rprior <- rprior |>
          add(bquote(p[[.(name_df)]] <- prior$rprior_df())) |>
          add(bquote(p[[.(name)]] <- prior$rprior(p[[.(name_df)]])))
      } else {
        rprior <- add(rprior, bquote(p[[.(name)]] <- prior$rprior()))
      }
    },
    exp = {
      rprior <- add(rprior, bquote(p[[.(name)]] <- prior$rprior()))
    }
  )
  rprior <- add(rprior, quote(p))

  if (e$prior.only) return(environment())

  # BEGIN draw function
  draw <- if (debug) function(p) {browser()} else function(p) {}

  # define function to compute partial variance factors (cf. partial residuals)
  if (e$single.V.block) {
    switch(e$Q0.type,
      unit = {
        if (e$sigma.fixed)
          get_partial_factor <- function(p) crossprod_mv(X, p[["e_"]]^2)
        else
          get_partial_factor <- function(p) crossprod_mv(X, p[["e_"]]^2) * (1 / p[["sigma_"]]^2)
      },
      diag = {
        if (e$sigma.fixed)
          get_partial_factor <- function(p) crossprod_mv(X, e[["Q0"]]@x * p[["e_"]]^2)
        else
          get_partial_factor <- function(p) crossprod_mv(X, e[["Q0"]]@x * p[["e_"]]^2) * (1 / p[["sigma_"]]^2)
      },
      symm = {
        # create list with blocks of Q0 corresponding to the subdivision by factor
        Q0.list <- list()
        X.ind <- list()
        fac <- numeric(q)
        for (i in seq_len(q)) {
          X.ind[[i]] <- which(X@perm == i - 1L)  # NB tabMatrix 0-based
          Q0.list[[i]] <- e[["Q0"]][X.ind[[i]], X.ind[[i]]]
        }
        get_partial_factor <- function(p) {
          for (i in seq_along(Q0.list)) {
            res <- p[["e_"]][X.ind[[i]]]
            fac[i] <- dotprodC(res, Q0.list[[i]] %m*v% res)
          }
        }
        get_partial_factor <- add(get_partial_factor, bquote(.(if (e$sigma.fixed) quote(fac) else quote(fac * (1 / p[["sigma_"]]^2)))))
      }
    )
  } else {
    switch(e$Q0.type,
      unit=, diag = {
        if (e$sigma.fixed)
          get_partial_factor <- function(p) crossprod_mv(X, p[["Q_"]] * p[["e_"]]^2) * p[[name]]
        else
          get_partial_factor <- function(p) crossprod_mv(X, p[["Q_"]] * p[["e_"]]^2) * p[[name]] * (1 / p[["sigma_"]]^2)
      },
      stop("TBI")
    )
  }

  switch(prior$type,
    invchisq = {
      if (is.list(prior$df)) {
        draw <- add(draw, bquote(p[[.(name_df)]] <- prior$draw_df(p[[.(name_df)]], 1 / p[[.(name)]])))
        if (is.list(prior$scale))
          draw <- add(draw, bquote(lambda <- prior$draw(p[[.(name_df)]], nh, get_partial_factor(p), 1 / p[[.(name)]])))
        else
          draw <- add(draw, bquote(lambda <- prior$draw(p[[.(name_df)]], nh, get_partial_factor(p))))
      } else {
        if (is.list(prior$scale))
          draw <- add(draw, bquote(lambda <- prior$draw(nh, get_partial_factor(p), 1 / p[[.(name)]])))
        else
          draw <- add(draw, bquote(lambda <- prior$draw(nh, get_partial_factor(p))))
      }
    },
    exp = {
      ph <- 1 - nh/2
      ah <- 2 / prior$scale
      draw <- add(draw, quote(lambda <- prior$draw(ph, ah, get_partial_factor(p))))
    }
  )

  if (e$single.V.block) {
    switch(e$Q0.type,
      unit = {
        draw <- add(draw, quote(p$Q_ <- X %m*v% (1 / lambda)))
      },
      diag = {
        draw <- add(draw, quote(p$Q_ <- e[["Q0"]]@x * (X %m*v% (1 / lambda))))
      },
      symm = {
        draw <- draw |>
          add(quote(p$Q_ <- X %m*v% (1 / lambda))) |>
          add(quote(p$QM_ <- block_scale_dsCMatrix(e[["Q0"]], p[["Q_"]])))
      }
    )
  } else {
    switch(e$Q0.type,
      unit=, diag = {
        draw <- add(draw, bquote(p$Q_ <- p[["Q_"]] * (X %m*v% (p[[.(name)]] / lambda))))
      },
      stop("TBI")
    )
  }

  draw <- draw |>
    add(bquote(p[[.(name)]] <- lambda)) |>
    add(quote(p))
  # END draw function

  if (prior$type == "invchisq" && is.list(prior$df) && prior$df$adapt) {
    # adaptation function that tunes MH proposal based on acceptance rates
    adapt <- function(ar) {
      if (ar[[name_df]] < .2)
        prior$df$tau <<- prior$df$tau * runif(1L, 0.6, 0.9)
      else if (ar[[name_df]] > .7)
        prior$df$tau <<- prior$df$tau * runif(1L, 1.1, 1.5)
    }
  }

  start <- function(p) {}
  if (prior$type == "invchisq" && is.list(prior$df)) {
    start <- start |>
      add(bquote(if (is.null(p[[.(name_df)]])) p[[.(name_df)]] <- rgamma(1L, prior$df$alpha0, prior$df$beta0))) |>
      add(bquote(if (length(p[[.(name_df)]]) != 1L) stop("wrong length for start value '", name_df, "'")))
  }
  # starting values 1 seem better than drawing from prior
  start <- start |>
    add(bquote(if (is.null(p[[.(name)]])) p[[.(name)]] <- rep.int(1, .(q)))) |>
    add(bquote(if (length(p[[.(name)]]) != .(q)) stop("wrong length for start value '", name, "'"))) |>
    add(quote(p))

  if (!e$single.V.block || e$Q0.type == "unit") rm(e)
  environment()
}
