
# mc: gen model component with AR1 component whose parameter should be inferred
# type: "TN" for independent truncated normal proposal
#       "RWTN" for truncated normal random walk proposal
#       "RWN" for normal random walk proposal
AR1_sampler <- function(mc) {
  opts <- mc$info$extra[[mc$AR1.inferred]]
  if (is.null(opts$control)) {
    # set default proposal
    if (is.null(mc$priorA)) {
      # independent truncated normal prior
      MH <- set_MH(type="TN")
    } else {
      # random walk truncated normal proposal
      MH <- set_MH(type="RWTN", scale=0.025, adaptive=TRUE, l=-1, u=1)
    }
  } else {
    MH <- opts$control
    if (!is.environment(MH)) stop("AR1: 'control' argument must be an environment created with function set_MH")
    MH$type <- match.arg(MH$type, c("TN", "RWTN", "RWN", "unif"))
    if (any(MH$type == c("RWTN", "unif"))) {
      if (is.null(MH$l)) MH$l <- -1
      if (is.null(MH$u)) MH$u <- 1
      if (MH$l > MH$u) stop("lower bound of MH truncated proposal exceeds upper bound")
    }
  }

  if (MH$type == "TN") {
    if (is.null(mc$priorA)) {
      Qv <- mc$generate_Qv()
      Q0.25 <- mc$kron_prod(mc$QA.template$update(0.25), Qv)
      Q0.5 <- mc$kron_prod(mc$QA.template$update(0.5), Qv)
    } else {
      stop("non-normal random effects, please choose random walk proposal 'RWTN' or 'RWN'")
    }
  }

  if (opts$prior$type == "unif") {
    tnprior <- FALSE
    l <- opts$prior$min
    u <- opts$prior$max
  } else if (opts$prior$type == "truncnormal") {
    tnprior <- TRUE
    mu0 <- opts$prior$mean
    prec0 <- opts$prior$precision
    l <- opts$prior$lower
    u <- opts$prior$upper
  }

  AR1.exp <- 0.5 * mc[["q"]] / mc$info$n[mc$AR1.inferred]

  name_v <- mc[["name"]]

  if (MH$type == "TN") {
    # independent truncated normal proposal, no adaptation in this case
    # here we need Q0.25 and Q0.5 to build up a template
    # row indices for -phi:
    ind1 <- which(abs(Q0.25@x - 0.5 * Q0.5@x) < sqrt(.Machine$double.eps))
    # row indices for 1 + phi^2:
    ind2 <- which(abs(Q0.25@x - ((1 + 0.25^2) / (1 + 0.5^2)) * Q0.5@x) < sqrt(.Machine$double.eps))
    col.ind <- rep(seq_len(ncol(Q0.5)), diff(Q0.5@p))  # col pointer --> col index
    # row, col indices of -phi elements
    i1 <- Q0.5@i[ind1] + 1L
    j1 <- col.ind[ind1]
    if (sum(i1 == j1) > 0) stop("unexpected: phi at diagonal")
    # row, col indices of (1 + phi^2) elements
    i2 <- Q0.5@i[ind2] + 1L
    j2 <- col.ind[ind2]
    diag.elements <- which(i2 == j2)
    i2.diag <- i2[diag.elements]
    j2.diag <- j2[diag.elements]
    i2 <- i2[-diag.elements]
    j2 <- j2[-diag.elements]
    ind2.diag <- ind2[diag.elements]
    ind2 <- ind2[-diag.elements]

    rm(Q0.25, Q0.5, col.ind, diag.elements, opts)

    draw <- function(phi, p) {
      v <- p[[name_v]]
      Qx <- p[[mc[["name_Q"]]]]
      # compute prec=alpha and beta in AR1 normal kernel MH proposal:
      # phi ~ p(phi) (1 - phi^2)^(q0 * l/l_AR1) exp(-0.5 * (alpha * phi^2 - 2 * beta * phi))
      prec <- (dotprodC(v[i2.diag], Qx[ind2.diag] * v[j2.diag]) +
                 2 * dotprodC(v[i2], Qx[ind2] * v[j2]) ) / (1 + phi^2)
      if (abs(phi) < sqrt(.Machine$double.eps))
        beta <- 0
      else
        beta <- - dotprodC(v[i1], Qx[ind1] * v[j1]) / phi
      mu <- beta/prec
      if (tnprior) {
        mu <- (prec0*mu0 + prec*mu)/(prec0 + prec)
        prec <- prec0 + prec
      }
      stdev <- 1/sqrt(prec)
      phi.star <- mu + stdev * Crtuvn((l - mu)/stdev, (u - mu)/stdev)
      log.ar <- AR1.exp * log((1 - phi.star^2) / (1 - phi^2))
      if (log(runif(1L)) < log.ar) phi.star else phi
    }
  } else {
    # random walk or uniform proposal

    lb <- max(-1, l)
    ub <- min(1, u)

    if (!is.null(mc$priorA)) name_omega <- mc[["name_omega"]]

    draw <- function(phi, p) {
      phi.star <- MH$propose(phi)
      if (phi.star < lb || phi.star > ub) return(phi)
      if (is.null(mc$priorA))
        QA.star <- mc$QA.template$update(phi.star)
      else
        QA.star <- crossprod_sym(mc$DA.template$update(phi.star), 1 / p[[name_omega]])
      Q.star <- mc$kron_prod(QA.star, p[[mc[["name_Qv"]]]])  # TODO kron_prod only available for components in a block?
      Q <- Q.star
      attr(Q, "x") <- p[[mc[["name_Q"]]]]
      v <- p[[name_v]]
      log.ar.post <- AR1.exp * log((1 - phi.star^2) / (1 - phi^2)) +
        0.5 * (dotprodC(v, Q %m*v% v) - dotprodC(v, Q.star %m*v% v))
      if (tnprior)
        log.ar.post <- log.ar.post + 0.5 * prec0 * ((phi - mu0)^2 - (phi.star - mu0)^2)
      if (MH$MH_accept(phi.star, phi, log.ar.post)) phi.star else phi
    }
  }

  environment()
}
