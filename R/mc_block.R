
# mcs: sublist of mod components to be sampled in a block
# e: environment of create_sampler
create_mc_block <- function(mcs, e=parent.frame()) {
  type <- "block"
  name <- "coef_"
  debug <- any(b_apply(mcs, `[[`, "debug"))

  if (e$family[["family"]] == "gamma") {
    modus <- "gamma"       # model for log(mean) of gamma
  } else if (all(s_apply(mcs, `[[`, "name") %in% names(e[["Vmod"]]))) {
    if (e$family[["family"]] == "gaussian_gamma")
      modus <- "vargamma"  # model for log(var) of gaussian and log(mean) of gamma
    else
      modus <- "var"       # model for log(var) of gaussian
  } else {
    modus <- "regular"     # model for mean of gaussian/binomial/...
  }

  if (e$control[["auto.order.block"]]) {
    mcs <- local({
      # order the components such that sparse matrices come first; may help find a better Cholesky permutation
      o <- unname(which(vapply(mcs, \(mc) isDiagonal(mc[["X"]]), TRUE)))  # start with diagonal matrices
      if (length(o))
        o <- c(o, seq_along(mcs)[-o][order(vapply(mcs[-o], \(mc) sparsity(mc[["X"]]), 1), decreasing=TRUE)])
      else
        o <- order(vapply(mcs, \(mc) sparsity(mc[["X"]]), 1), decreasing=TRUE)
      mcs[o]
    })
  }

  X <- matrix(0, e[["n"]], 0L)  # block design matrix
  ind <- 0L
  for (mc in mcs) {
    if (mc[["type"]] == "gen" && mc[["gl"]])
      X <- cbind(X, mc[["X"]], zeroMatrix(e[["n"]], mc$glp[["q"]]))
    else
      X <- cbind(X, mc[["X"]])
    mc$block.i <- (ind + 1L):ncol(X)
    ind <- ncol(X)
  }
  rm(ind)
  X <- economizeMatrix(X, check=FALSE)
  q <- ncol(X)

  # template for (updating) blocked precision matrix
  # each mc$Q for gen is either ddi or dsC; the same holds true for mc$Q0 for reg and mec
  get_Q <- function(mc) if (mc[["type"]] == "gen") mc[["Q"]] else mc[["Q0"]]
  if (all(b_apply(mcs, \(mc) class(get_Q(mc))[[1L]] == "ddiMatrix"))) {
    QT <- Cdiag(unlst(lapply(mcs, \(mc) ddi_diag(get_Q(mc)))))
  } else {
    QT <- local({
      x <- NULL
      i <- NULL
      size <- 0L
      p <- 0L
      for (mc in mcs) {
        Q <- get_Q(mc)
        if (class(Q)[1L] == "ddiMatrix") {
          x <- c(x, ddi_diag(Q))
          i <- c(i, size + seq_len(nrow(Q)) - 1L)
          p <- c(p, p[length(p)] + seq_len(nrow(Q)))
        } else {
          x <- c(x, Q@x)
          i <- c(i, size + Q@i)
          p <- c(p, p[length(p)] + Q@p[-1L])
        }
        size <- size + nrow(Q)
      }
      new("dsCMatrix", i=i, p=p, x=x, uplo="U", Dim=c(size, size))
    })
  }
  rm(get_Q)
  # individual Q matrices no longer needed (we still have kron_prod closures)
  for (mc in mcs) if (mc[["type"]] == "gen") rm("Q", envir=mc)
  get_Qvector <- function(p, tau) {
    Qvector <- NULL
    for (mc in mcs)
      if (mc[["type"]] == "gen")
        Qvector <- c(Qvector, p[[mc$name_Q]])
      else
        Qvector <- c(Qvector, tau * mc[["Q0"]]@x)
    Qvector
  }

  if (modus == "regular") {
    # non-zero prior means of reg or mec components
    nonzero.mean <- any(b_apply(mcs, \(mc) any(mc[["type"]] == c("reg", "mec")) && !mc[["zero.mean"]]))
    if (nonzero.mean) {
      Q0b0 <- numeric(q)
      for (mc in mcs) {
        if (any(mc[["type"]] == c("reg", "mec")) && !mc[["zero.mean"]])
          Q0b0[mc$block.i] <- mc[["Q0b0"]]
      }
    }
  }

  if (is.null(e$control[["CG"]])) {
    if (e[["modeled.Q"]]) {
      XX <- crossprod_sym(X, crossprod_sym(Cdiag(runif(e[["n"]], 0.9, 1.1)), e[["Q0"]]))
    } else {
      XX <- economizeMatrix(crossprod_sym(X, e[["Q0"]]),
        symmetric=TRUE, drop.zeros=TRUE
      )
    }

    # derive constraint matrix, if any
    if (any(b_apply(mcs, \(mc) !is.null(mc[["R"]])))) {
      R <- zeroMatrix(0L, 0L)
      r <- NULL
      for (mc in mcs) {
        if (is.null(mc[["R"]])) {
          if (mc[["type"]] == "gen" && mc[["gl"]])
            R <- rbind(R, zeroMatrix(mc[["q"]] + mc$glp[["q"]], ncol(R)))
          else
            R <- rbind(R, zeroMatrix(mc[["q"]], ncol(R)))
        } else {
          if (mc[["type"]] == "gen" && mc[["gl"]]) {
            R <- bdiag(R, mc$glp[["R"]])
            r <- c(r, rep(if (is.null(mc$glp[["r"]])) 0 else mc$glp[["r"]], ncol(mc$glp[["R"]])))
          } else {
            R <- bdiag(R, mc[["R"]])
            r <- c(r, if (is.null(mc[["r"]])) rep(0, ncol(mc[["R"]])) else mc[["r"]])
          }
        }
      }
      if (nrow(R) != q) stop("incompatible dimensions of design and constraint matrices")
      # TODO remove individual R matrices as they are not needed in the single block sampler
      #      add support for additional constraints defined over the whole coefficient vector
      # In the case of constraints the XX + Q matrix often becomes singular
      # sampling from such a IGMRF can be done by first adding a multiple of tcrossprod(R) to it (Rue and Held, 2005)
      # most convenient is to add it to XX;
      # But this is not always required. Better add as little as is necessary to get pd (takes up lots of memory for lengthy random walks ...)
      # Another option is to add a multiple of I to XX + Q and correct with MH
    } else {
      R <- NULL
    }

    if (any(b_apply(mcs, \(mc) !is.null(mc[["S"]])))) {
      S <- zeroMatrix(0L, 0L)
      s <- NULL
      for (mc in mcs) {
        if (is.null(mc[["S"]]))
          S <- rbind(S, zeroMatrix(mc[["q"]], ncol(S)))
        else {
          S <- bdiag(S, mc[["S"]])
          s <- c(s, if (is.null(mc[["s"]])) rep(0, ncol(mc[["S"]])) else mc[["s"]])
        }
      }
      if (nrow(S) != q) stop("incompatible dimensions of design and constraint matrices")
      # TODO remove individual S matrices as they are not needed in the single block sampler
      #      add support for additional constraints defined over the whole coefficient vector
    } else {
      S <- NULL
    }

    if (modus == "regular") {
      sparse_template(environment(), update.XX=e[["modeled.Q"]] || any(s_apply(mcs, `[[`, "type") == "mec"),
                      control=e[["control"]])
    } else {
      if (is.null(R)) {
        # TODO include in X0 the fixed part of QT (from reg components)
        mat_sum <- make_mat_sum(M0 = if (modus == "vargamma") 2 * XX else XX, M1=QT)
        cholQ <- build_chol(mat_sum(QT))
      } else {
        stop("not supported: blocked sampler for variance model components or gamma/gaussian_gamma family with constraints")
      }
    }
  } else {
    # TODO check that the only constraints are IGMRF equality constraints
    S <- NULL
    CG_sampler <- setup_CG_sampler(mbs=mcs, X=X, sampler=e, control=e$control[["CG"]])
  }

  if (e[["compute.weights"]]) {
    # form the q_all x m matrix corresponding to the linear predictor as represented componentwise in linpred
    # TODO do we really need to store both X and t(X) in this case?
    linpred <- if (is.null(e[["linpred"]]))
      economizeMatrix(t(X), strip.names=FALSE)
    else
      economizeMatrix(t(do.call(cbind, lapply(e[["linpred"]][names(mcs)], \(x) x[["Xnew"]]))), allow.tabMatrix=FALSE)
  }

  if (e[["prior.only"]]) return(environment())

  # BEGIN draw function
  draw <- if (debug) function(p) {browser()} else function(p) {}
  if (modus == "var" || modus == "vargamma") {
    if (!e[["single.V.block"]])
      for (m in seq_along(mcs)) {
        # TODO linpred function (NB name clash) --> apply exp() only once
        #      or exploit index design matrices to reduce cost of exp()
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(mcs[[.(m)]]$lp(p))))
      }
  } else if (e[["single.block"]]) {
    if (e[["e.is.res"]])
      draw <- add(draw, quote(p$e_ <- e$y_eff()))
    # otherwise p$e_ = e$y_eff() = 0
  } else {
    for (m in seq_along(mcs)) {
      # residuals could also be computed using the block X,
      #   but this way it is usually faster due to more efficient matrix types
      if (e[["e.is.res"]])
        draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], TRUE, p)))
      else
        draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], FALSE, p)))
    }
  }
  draw <- draw |>
    add(bquote(tau <- .(if (e[["sigma.fixed"]]) 1 else quote(1 / p[["sigma_"]]^2)))) |>
    # update the block-diagonal joint precision matrix
    add(quote(attr(QT, "x") <- get_Qvector(p, tau)))

  if (modus == "regular") {
    if (!is.null(S)) {  # need to reconstruct coef_ as input to TMVN sampler
      # TODO store coef_ component and only replace the subcomponents with PX
      #      and check whether this works in case of gen component with gl=TRUE
      draw <- add(draw, quote(
        for (mc in mcs) p[["coef_"]][mc$block.i] <- p[[mc$name]]
      ))
    }
    # update mec component columns of X
    # TODO more efficient update of only those elements that can change (for dgC or matrix X)
    for (mc in mcs)
      if (mc[["type"]] == "mec")
        draw <- add(draw, bquote(X[, mcs[[.(mc[["name"]])]]$block.i] <- p[[.(mc[["name_X"]])]]))
    if (e[["single.block"]] && !e[["modeled.Q"]] && !any(s_apply(mcs, `[[`, "type") == "mec") && e$family[["link"]] != "probit") {
      Xy <- crossprod_mv(X, e[["Q0"]] %m*v% e$y_eff())
      if (nonzero.mean) {
        Xy <- Xy + Q0b0
        rm(Q0b0)
      }
    } else {
      if (nonzero.mean)
        draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p)) + Q0b0))
      else
        draw <- add(draw, quote(Xy <- crossprod_mv(X, e$Q_e(p))))
    }
    if (is.null(e$control[["CG"]])) {
      if (e[["modeled.Q"]]) {
        if (e[["Q0.type"]] == "symm")
          draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
        else {
          cps_template <- NULL
          if (inherits(X, "dgCMatrix")) {
            tryCatch(
              cps_template <- sparse_crossprod_sym_template(X, e$control[["max.size.cps.template"]]),
              error = function(e) {
                # template too large
                NULL
              }
            )
          }
          if (is.null(cps_template)) {
            draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
          } else {
            draw <- add(draw, quote(XX <- cps_template(p[["Q_"]])))
          }
        }
      } else if (any(s_apply(mcs, `[[`, "type") == "mec")) {
        draw <- add(draw, quote(XX <- crossprod_sym(X, e[["Q0"]])))
      }
      draw <- add(draw, quote(Qlist <- update(XX, QT, 1, 1/tau)))
      if (e$control[["cMVN.sampler"]]) {
        draw <- add(draw, bquote(coef <- MVNsampler$draw(p, Xy=Xy, X=X, QT=Qlist[["Q"]])[[.(name)]]))
      } else {
        draw <- add(draw, bquote(coef <- MVNsampler$draw(p, .(if (e[["sigma.fixed"]]) 1 else quote(p[["sigma_"]])), Q=Qlist[["Q"]], Imult=Qlist[["Imult"]], Xy=Xy)[[.(name)]]))
      }
    } else {
      draw <- add(draw, bquote(CGstart <- numeric(.(q))))
      for (mc in mcs)
        draw <- add(draw, bquote(CGstart[mcs[[.(mc[["name"]])]]$block.i] <- p[[.(mc[["name"]])]]))
      draw <- add(draw, quote(coef <- CG_sampler$draw(p, Xy, X, QT, e, start=CGstart)))
    }
  } else {
    if (modus == "var" || modus == "vargamma") {  # variance modelling
      if (e[["single.V.block"]])
        draw <- add(draw, bquote(vkappa <- .(if (e[["sigma.fixed"]]) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2))
      else
        draw <- add(draw, bquote(vkappa <- .(if (e[["sigma.fixed"]]) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2 * p[["Q_"]]))
    }
    if (modus == "gamma") {
      if (e$family[["alpha.fixed"]]) {
        alpha <- e$family$get_shape()
        if (e[["single.block"]]) {
          kappa <- alpha * e[["y"]]
        } else {
          kappa0 <- alpha * e[["y"]]
          draw <- add(draw, quote(kappa <- kappa0 * exp(-p[["e_"]])))
        }
      } else {
        draw <- add(draw, quote(alpha <- e[["family"]]$get_shape(p)))
        if (e[["single.block"]])
          draw <- add(draw, quote(kappa <- alpha * e[["y"]]))
        else
          draw <- add(draw, quote(kappa <- alpha * e[["y"]] * exp(-p[["e_"]])))
      }
    } else if (modus == "vargamma") {
      if (e$family[["alpha.fixed"]]) {
        alpha <- e[["family"]]$get_shape()
        if (e[["single.V.block"]]) {
          kappa <- alpha * e$family[["sigmasq"]]
        } else {
          kappa0 <- alpha * e$family[["sigmasq"]]
          draw <- add(draw, quote(kappa <- kappa0 * p[["Q_"]]))
        }
      } else {
        draw <- add(draw, quote(alpha <- e[["family"]]$get_shape(p)))
        if (e[["single.V.block"]]) {
          draw <- add(draw, quote(kappa <- alpha * e$family[["sigmasq"]]))
        } else {
          draw <- add(draw, quote(kappa <- alpha * e$family[["sigmasq"]] * p[["Q_"]]))
        }
      }
    }
    draw <- add(draw, quote(cholQ$update(mat_sum(QT))))  # TODO if only reg components cholQ is fixed
    if (modus == "var" || modus == "vargamma")
      draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(e[["n"]]), 0.5, vkappa))))
    if (modus == "gamma")
      draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa))))
    else if (modus == "vargamma")
      draw <- add(draw, bquote(Hz <- Hz + crossprod_mv(X, rMLiG(.(e[["n"]]), alpha, kappa))))
    # prior contributions from mcs
    draw <- add(draw, quote(
      for (m in seq_along(mcs)) {
        mc <- mcs[[m]]
        switch(mc[["type"]],
          reg =
            if (mc[["informative.prior"]]) {
              if (mc[["zero.mean"]])
                z <- rMLiG(mc[["q"]], mc$prior[["a"]], mc$prior[["a"]])
              else
                z <- rMLiG(mc[["q"]], mc$prior[["a"]], mc$prior[["a"]] * exp(sqrt(mc$prior$precision/mc$prior$a) * mc$prior$mean))
              Hz[mc$block.i] <- Hz[mc$block.i] + sqrt(mc$prior$precision/mc$prior[["a"]]) * z
            },
          gen = {
            z <- rMLiG(mc[["q"]], mc[["a"]], mc[["a"]])
            Hz[mc$block.i] <- Hz[mc$block.i] + z / (p[[mc$name_sigma]] * sqrt(mc[["a"]]))
          },
          stop("only reg and gen components supported")
        )
      }
    ))
    draw <- add(draw, quote(coef <- cholQ$solve(Hz)))
  }

  # split coef and assign to the separate coefficient batches
  for (m in seq_along(mcs)) {
    if (mcs[[m]]$type == "gen" && mcs[[m]]$gl) {
      draw <- draw |>
        add(bquote(u <- coef[mcs[[.(m)]]$block.i])) |>
        add(bquote(p[[.(mcs[[m]]$name)]] <- u[mcs[[.(m)]]$i.v])) |>
        add(bquote(p[[.(mcs[[m]]$name_gl)]] <- u[mcs[[.(m)]]$i.alpha]))
    } else {
      draw <- add(draw, bquote(p[[.(mcs[[m]]$name)]] <- coef[mcs[[.(m)]]$block.i]))
    }
    if (modus == "var" || modus == "vargamma") {
      if (e[["single.V.block"]])
        draw <- add(draw, bquote(p[["Q_"]] <- exp(-mcs[[.(m)]]$lp(p))))
      else
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(-mcs[[.(m)]]$lp(p))))
    } else {
      if (e[["e.is.res"]]) {
        draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], FALSE, p)))
      } else {
        if (m == 1L && e[["single.block"]]) {
          # in this case p$e_ = e$y_eff() = 0
          draw <- add(draw, quote(p$e_ <- mcs[[1L]]$lp(p)))
        } else {
          draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], TRUE, p)))
        }
      }
    }
  }
  if (e[["compute.weights"]]) {
    # TODO solve-sparse method that returns dense
    draw <- add(draw, quote(p$weights_ <- X %m*m% as.matrix(MVNsampler$cholQ$solve(linpred))))
    if (e[["modeled.Q"]]) {
      if (e[["Q0.type"]] == "symm")
        draw <- add(draw, quote(p$weights_ <- p[["QM_"]] %m*m% p[["weights_"]]))
      else
        draw <- add(draw, quote(p$weights_ <- p[["Q_"]] * p[["weights_"]]))
    } else {
      if (e[["Q0.type"]] != "unit") {
        draw <- add(draw, quote(p$weights_ <- e[["Q0"]] %m*m% p[["weights_"]]))
      }
    }
  }
  draw <- add(draw, quote(p))
  # END draw function

  start <- function(p) {}
  if (is.null(e$control[["CG"]]) && !e$control[["cMVN.sampler"]]) {
    if (modus == "regular") {
      start <- add(start, bquote(coef <- MVNsampler$start(p, e[["scale.sigma"]])[[.(name)]]))
    } else {
      # account for scaling of covariates
      start <- add(start, bquote(coef <- Crnorm(.(q)) / colwise_maxabs(X)))
    }
    start <- add(start, quote(
      for (mc in mcs) {
        if (mc[["type"]] == "gen" && mc[["gl"]]) {
          u <- coef[mc$block.i]
          if (is.null(p[[mc$name]])) p[[mc$name]] <- u[mc$i.v]
          if (is.null(p[[mc$name_gl]])) p[[mc$name_gl]] <- u[mc$i.alpha]
        } else {
          if (is.null(p[[mc$name]])) p[[mc$name]] <- coef[mc$block.i]
        }
      }
    ))
  } else {
    start <- add(start, quote(
      for (mc in mcs) {
        if (mc[["type"]] == "gen" && mc[["fastGMRFprior"]]) {
          Qv <- rexp(1L)
          if (is.null(mc$rGMRFprior))
            setup_priorGMRFsampler(mc, Qv)
          p[[mc[["name"]]]] <- mc$rGMRFprior(Qv)
        } else {
          p[[mc[["name"]]]] <- Crnorm(mc[["q"]], sd=e[["scale.sigma"]])
        }
      }
    ))
  }
  # TODO check formats of user-provided start values
  start <- add(start, quote(p))

  rm(mc)
  environment()
}
