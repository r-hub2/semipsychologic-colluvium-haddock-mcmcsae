
# mcs: sublist of mod components to be sampled in a block
create_mc_block <- function(mcs, fam, sc, prior.only, 
                            compute.weights=FALSE, linpred=NULL) {
  type <- "block"
  name <- "coef_"
  debug <- any(b_apply(mcs, `[[`, "debug"))

  if (fam[["family"]] == "gamma") {
    modus <- "gamma"       # model for log(mean) of gamma
  } else if (all(s_apply(mcs, `[[`, "name") %in% names(fam[["Vmod"]]))) {
    if (fam[["family"]] == "gaussian_gamma")
      modus <- "vargamma"  # model for log(var) of gaussian and log(mean) of gamma
    else
      modus <- "var"       # model for log(var) of gaussian
  } else {
    modus <- "regular"     # model for mean of gaussian/binomial/...
  }

  if (sc[["auto.order.block"]]) {
    mcs <- local({
      # order the components such that sparse matrices come first, to help find a better Cholesky permutation
      o <- whichv(vapply(mcs, \(mc) isDiagonal(mc[["X"]]), TRUE), TRUE)  # start with diagonal matrices
      if (length(o))
        o <- c(o, seq_along(mcs)[-o][order(vapply(mcs[-o], \(mc) sparsity(mc[["X"]]), 1), decreasing=TRUE)])
      else
        o <- order(vapply(mcs, \(mc) sparsity(mc[["X"]]), 1), decreasing=TRUE)
      mcs[o]
    })
  }

  X <- matrix(0, fam[["n"]], 0L)  # block design matrix
  ind <- 0L
  for (mc in mcs) {
    if (mc[["type"]] == "gen" && mc[["gl"]])
      X <- cbind(X, mc[["X"]], zeroMatrix(fam[["n"]], mc$glp[["q"]]))
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
  QT <- bdiag_ddidsC(
    lapply(mcs, \(mc) if (mc[["type"]] == "gen") mc[["Q"]] else mc[["Q0"]])
  )
  # individual Q matrices no longer needed (we still have kron_prod closures)
  for (mc in mcs) if (mc[["type"]] == "gen") rm("Q", envir=mc)
  if (any(b_apply(mcs, \(mc) is.function(mc[["get_Q"]])))) {
    QT.ind <- 0L  # 0-based indices of components' contributions to QT@x
    if (length(mcs) > 1L) QT.ind <- c(QT.ind,
      if (class(QT)[1L] == "ddiMatrix")
        cumsum(i_apply(mcs[seq_len(length(mcs) - 1L)], `[[`, "q"))
      else
        QT@p[cumsum(i_apply(mcs[seq_len(length(mcs) - 1L)], `[[`, "q")) + 1L]
    )
    i.getQ <- b_apply(mcs, function(x) is.function(x[["get_Q"]]))
    Q.funs <- lapply(mcs[i.getQ], function(x) x$get_Q)  # direct refs to existing get_Q
    QT.ind <- QT.ind[i.getQ]
    rm(i.getQ)
    get_Qvector <- function(p) {
      Qvector <- copy_obj(QT@x)
      for (i in seq_along(QT.ind)) set_in_place(Qvector, QT.ind[i], Q.funs[[i]](p))
      Qvector
    }
  } else {
    get_Qvector <- NULL
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

  if (is.null(sc[["CG"]])) {
    if (fam[["modeled.Q"]]) {
      XX <- crossprod_sym(X, crossprod_sym(Cdiag(runif(fam[["n"]], 0.9, 1.1)), fam[["Q0"]]))
    } else {
      XX <- economizeMatrix(crossprod_sym(X, fam[["Q0"]]),
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
      sparse_template(environment(), update.XX=fam[["modeled.Q"]] || any(s_apply(mcs, `[[`, "type") == "mec"),
                      control=sc)
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
    CGsampler <- setup_CG_sampler(mbs=mcs, X=X, fam=fam, control=sc[["CG"]])
  }

  if (compute.weights) {
    # form the q_all x m matrix corresponding to the linear predictor as represented componentwise in linpred
    # TODO do we really need to store both X and t(X) in this case?
    linpred <- if (is.null(linpred))
      economizeMatrix(t(X), strip.names=FALSE)
    else
      economizeMatrix(t(do.call(cbind, lapply(linpred[names(mcs)], \(x) x[["Xnew"]]))), allow.tabMatrix=FALSE)
  }

  if (prior.only) return(environment())

  # BEGIN draw function
  draw <- if (debug) function(p) {browser()} else function(p) {}
  if (modus == "var" || modus == "vargamma") {
    if (!fam[["single.V.block"]])
      for (m in seq_along(mcs)) {
        # TODO linpred function (NB name clash) --> apply exp() only once
        #      or exploit index design matrices to reduce cost of exp()
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(mcs[[.(m)]]$lp(p))))
      }
  } else if (sc[["single.block"]]) {
    if (fam[["e.is.res"]])
      draw <- add(draw, quote(p$e_ <- copy_obj(fam[["y"]])))
    # otherwise p$e_ = 0
  } else {
    for (m in seq_along(mcs)) {
      # residuals could also be computed using the block X,
      #   but this way it is usually faster due to more efficient matrix types
      if (fam[["e.is.res"]])
        draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], TRUE, p)))
      else
        draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], FALSE, p)))
    }
  }
  # update the block-diagonal joint precision matrix
  # this creates a new local matrix QT, with everything except x-slot shared with original QT
  if (is.function(get_Qvector)) draw <- add(draw, quote(attr(QT, "x") <- get_Qvector(p)))

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
    if (is.null(sc[["CG"]])) {
      if (fam[["modeled.Q"]]) {
        if (fam[["Q0.type"]] == "symm")
          draw <- add(draw, quote(XX <- crossprod_sym(X, p[["QM_"]])))
        else {
          cps_template <- NULL
          if (inherits(X, "dgCMatrix")) {
            tryCatch(
              cps_template <- sparse_crossprod_sym_template(X, sc[["max.size.cps.template"]]),
              error = function(e) {
                # template too large
                NULL
              }
            )
          }
          if (is.null(cps_template))
            draw <- add(draw, quote(XX <- crossprod_sym(X, p[["Q_"]])))
          else
            draw <- add(draw, quote(XX <- cps_template(p[["Q_"]])))
        }
      } else if (any(s_apply(mcs, `[[`, "type") == "mec")) {
        draw <- add(draw, quote(XX <- crossprod_sym(X, fam[["Q0"]])))
      }
      # NB multi-response family always has sigma.fixed=TRUE
      draw <- add(draw, bquote(update(XX, QT, 1, .(if (fam[["sigma.fixed"]]) 1 else quote(p[["sigma_"]]^2)))))
    }
    if (fam[["link"]] == "probit") {
      if (!is.null(sc[["CG"]])) fam$control[["probit.HaarPXDA"]] <- FALSE
      if (fam$control[["probit.HaarPXDA"]]) {
        if (!all(sapply(mcs, function(mc) if (any(c("reg", "mec") == mc[["type"]])) mc[["zero.mean"]] else TRUE)))
          warn("in case of coefficients with non-zero prior mean, the current implementation of ",
               "the default sampler for probit models may be biased; you may want to revert to ",
               "the standard Albert-Chib sampler in this case by setting probit.HaarPXDA=FALSE ",
               "via create_sampler's control argument")
        # Haar PX-DA sandwich step
        if (!is.null(sc[["cMVN.sampler"]])) {
          draw <- add(draw, quote(
            rate <- 0.5 * dotprodC(p[["z_"]], p[["z_"]] - X %m*v% MVNsampler$smplr$cholQ$solve(crossprod_mv(X, p[["z_"]]))[MVNsampler$smplr$Iq])
          ))
        } else {
          draw <- add(draw, quote(
            rate <- 0.5 * dotprodC(p[["z_"]], p[["z_"]] - X %m*v% MVNsampler$cholQ$solve(crossprod_mv(X, p[["z_"]])))
          ))
        }
        draw <- add(draw, quote(
          p[["z_"]] <- p[["z_"]] * sqrt(
            rgamma(1L, 
              shape = 0.5*fam[["n"]],
              rate = rate
            )
          )
        ))
      }
    }
    if (sc[["single.block"]] && !fam[["modeled.Q"]] && !any(s_apply(mcs, `[[`, "type") == "mec") && fam[["link"]] != "probit") {
      # necessarily e.is.res=TRUE, i.e. gaussian(-derived) family
      Xy <- crossprod_mv(X, fam[["Q0"]] %m*v% fam[["y"]])
      if (nonzero.mean) {
        Xy <- Xy + Q0b0
        rm(Q0b0)
      }
    } else {
      if (nonzero.mean)
        draw <- add(draw, quote(Xy <- crossprod_mv(X, fam$Q_e(p)) + Q0b0))
      else
        draw <- add(draw, quote(Xy <- crossprod_mv(X, fam$Q_e(p))))
    }
    if (is.null(sc[["CG"]])) {
      if (!is.null(sc[["cMVN.sampler"]]))
        draw <- add(draw, bquote(coef <- MVNsampler$draw(p, Xy=Xy, X=X)[[.(name)]]))
      else
        draw <- add(draw, bquote(coef <- MVNsampler$draw(p, .(if (fam[["sigma.fixed"]]) 1 else quote(p[["sigma_"]])), Xy=Xy)[[.(name)]]))
    } else {
      draw <- add(draw, bquote(CGstart <- numeric(.(q))))
      for (mc in mcs)
        draw <- add(draw, bquote(CGstart[mcs[[.(mc[["name"]])]]$block.i] <- p[[.(mc[["name"]])]]))
      draw <- add(draw, quote(coef <- CGsampler$draw(p, Xy, X, QT, fam, start=CGstart)))
    }
  } else {
    if (modus == "var" || modus == "vargamma") {  # variance modelling
      if (fam[["single.V.block"]])
        draw <- add(draw, bquote(vkappa <- .(if (fam[["sigma.fixed"]]) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2))
      else
        draw <- add(draw, bquote(vkappa <- .(if (fam[["sigma.fixed"]]) 0.5 else quote(0.5/p[["sigma_"]]^2)) * p[["e_"]]^2 * p[["Q_"]]))
    }
    if (modus == "gamma") {
      if (fam[["alpha.fixed"]]) {
        alpha <- fam$get_shape()
        if (sc[["single.block"]]) {
          kappa <- alpha * fam[["y"]]
        } else {
          kappa0 <- alpha * fam[["y"]]
          draw <- add(draw, quote(kappa <- kappa0 * exp(-p[["e_"]])))
        }
      } else {
        draw <- add(draw, quote(alpha <- fam$get_shape(p)))
        if (sc[["single.block"]])
          draw <- add(draw, quote(kappa <- alpha * fam[["y"]]))
        else
          draw <- add(draw, quote(kappa <- alpha * fam[["y"]] * exp(-p[["e_"]])))
      }
    } else if (modus == "vargamma") {
      if (fam[["alpha.fixed"]]) {
        alpha <- fam$get_shape()
        if (fam[["single.V.block"]]) {
          kappa <- alpha * fam[["sigmasq"]]
        } else {
          kappa0 <- alpha * fam[["sigmasq"]]
          draw <- add(draw, quote(kappa <- kappa0 * p[["Q_"]]))
        }
      } else {
        draw <- add(draw, quote(alpha <- fam$get_shape(p)))
        if (fam[["single.V.block"]]) {
          draw <- add(draw, quote(kappa <- alpha * fam[["sigmasq"]]))
        } else {
          draw <- add(draw, quote(kappa <- alpha * fam[["sigmasq"]] * p[["Q_"]]))
        }
      }
    }
    draw <- add(draw, quote(cholQ$update(mat_sum(QT))))  # TODO if only reg components cholQ is fixed
    if (modus == "var" || modus == "vargamma")
      draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(fam[["n"]]), 0.5, vkappa))))
    if (modus == "gamma")
      draw <- add(draw, bquote(Hz <- crossprod_mv(X, rMLiG(.(fam[["n"]]), alpha, kappa))))
    else if (modus == "vargamma")
      draw <- add(draw, bquote(Hz <- Hz + crossprod_mv(X, rMLiG(.(fam[["n"]]), alpha, kappa))))
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
    if (mcs[[m]][["type"]] == "gen" && mcs[[m]][["gl"]]) {
      draw <- draw |>
        add(bquote(u <- coef[mcs[[.(m)]][["block.i"]]])) |>
        add(bquote(p[[.(mcs[[m]]$name)]] <- u[mcs[[.(m)]][["i.v"]]])) |>
        add(bquote(p[[.(mcs[[m]]$name_gl)]] <- u[mcs[[.(m)]][["i.alpha"]]]))
    } else if (mcs[[m]][["type"]] == "s") {
      if (mcs[[m]][["qf"]] > 0L) {
        draw <- draw |>
          add(bquote(u <- coef[mcs[[.(m)]][["block.i"]]])) |>
          add(bquote(p[[.(mcs[[m]]$name.f)]] <- u[mcs[[.(m)]][["i.f"]]])) |>
          add(bquote(p[[.(mcs[[m]]$name.r)]] <- u[mcs[[.(m)]][["i.r"]]]))
      } else {
        draw <- add(draw, bquote(p[[.(mcs[[m]]$name.r)]] <- coef[mcs[[.(m)]][["block.i"]]]))
      }
    } else {
      draw <- add(draw, bquote(p[[.(mcs[[m]]$name)]] <- coef[mcs[[.(m)]]$block.i]))
    }
    if (modus == "var" || modus == "vargamma") {
      if (fam[["single.V.block"]])
        draw <- add(draw, bquote(p[["Q_"]] <- exp(-mcs[[.(m)]]$lp(p))))
      else
        draw <- add(draw, bquote(p[["Q_"]] <- p[["Q_"]] * exp(-mcs[[.(m)]]$lp(p))))
    } else {
      if (fam[["e.is.res"]]) {
        draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], FALSE, p)))
      } else {
        if (m == 1L && sc[["single.block"]]) {
          # case !e.is.res where p$e_ = 'y_eff' = 0
          draw <- add(draw, quote(p$e_ <- mcs[[1L]]$lp(p)))
        } else {
          draw <- add(draw, bquote(mcs[[.(m)]]$lp_update(p[["e_"]], TRUE, p)))
        }
      }
    }
  }
  if (compute.weights) {
    # TODO solve-sparse method that returns dense
    draw <- add(draw, quote(p$weights_ <- X %m*m% as.matrix(MVNsampler[["cholQ"]]$solve(linpred))))
    if (fam[["modeled.Q"]]) {
      if (fam[["Q0.type"]] == "symm")
        draw <- add(draw, quote(p$weights_ <- p[["QM_"]] %m*m% p[["weights_"]]))
      else
        draw <- add(draw, quote(p$weights_ <- p[["Q_"]] * p[["weights_"]]))
    } else {
      if (fam[["Q0.type"]] != "unit") {
        draw <- add(draw, quote(p$weights_ <- fam[["Q0"]] %m*m% p[["weights_"]]))
      }
    }
  }
  draw <- add(draw, quote(p))
  # END draw function

  start <- function(p) {}
  if (sc[["cMVN.or.CG"]]) {
    start <- add(start, quote(
      for (mc in mcs) {
        if (mc[["type"]] == "gen" && mc[["fastGMRFprior"]]) {
          Qv <- rexp(1L)
          if (is.null(mc$rGMRFprior))
            setup_priorGMRFsampler(mc, Qv)
          p[[mc[["name"]]]] <- check_and_get(p, mc[["name"]], mc[["q"]],
            \() if (is.null(mc[["priorA"]])) mc$rGMRFprior(Qv) else mc$rGMRFprior(Qv, rep.int(1, mc[["lD"]]))
          )
        } else {
          p[[mc[["name"]]]] <- check_and_get(p, mc[["name"]], mc[["q"]], \() Crnorm(mc[["q"]], sd=fam[["scale.sigma"]]))
        }
      }
    ))
  } else {
    if (modus == "regular") {
      start <- add(start, bquote(coef <- MVNsampler$start(p, fam[["scale.sigma"]])[[.(name)]]))
    } else {
      # account for scaling of covariates
      start <- add(start, bquote(coef <- Crnorm(.(q)) / colwise_maxabs(X)))
    }
    start <- add(start, quote(
      for (mc in mcs) {
        u <- coef[mc[["block.i"]]]
        if (mc[["type"]] == "gen" && mc[["gl"]]) {
          p[[mc$name]] <- check_and_get(p, mc[["name"]], mc[["q"]], \() u[mc[["i.v"]]])
          p[[mc$name_gl]] <- check_and_get(p, mc[["name_gl"]], mc$glp[["q"]], \() u[mc[["i.alpha"]]])
        } else if (mc[["type"]] == "s") {
          if (mc[["qf"]] > 0L) {
            p[[mc$name.f]] <- check_and_get(p, mc[["name.f"]], mc[["qf"]], \() u[mc[["i.f"]]])
            p[[mc$name.r]] <- check_and_get(p, mc[["name.r"]], mc[["qr"]], \() u[mc[["i.r"]]])
          } else {
            p[[mc$name.r]] <- check_and_get(p, mc[["name.r"]], mc[["qr"]], \() u)
          }
        } else {
          p[[mc$name]] <- check_and_get(p, mc[["name"]], mc[["q"]], \() u)
        }
      }
    ))
  }
  start <- add(start, quote(p))

  rm(mc)
  environment()
}
