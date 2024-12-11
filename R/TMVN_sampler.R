#' Set up a sampler object for sampling from a possibly truncated and degenerate multivariate normal distribution
#'
#' This function sets up an object for multivariate normal sampling based on a specified precision matrix.
#' Linear equality and inequality restrictions are supported.
#' For sampling under inequality restrictions four algorithms are available. The default in that case is
#' an exact Hamiltonian Monte Carlo algorithm (Pakman and Paninski, 2014). A related algorithm is the zig-zag
#' Hamiltonian Monte Carlo method (Nishimura et al., 2021) in which momentum is sampled from a Laplace instead
#' of normal distribution. Alternatively, a Gibbs sampling algorithm can be used (Rodriguez-Yam et al., 2004).
#' The fourth option is a data augmentation method that samples from a smooth approximation to the truncated
#' multivariate normal distribution (Souris et al., 2018).
#'
#' The componentwise Gibbs sampler uses univariate truncated normal samplers as described
#' in Botev and L'Ecuyer (2016). These samplers are implemented in R package \pkg{TruncatedNormal},
#' but here translated to C++ for an additional speed-up.
#'
#' @examples
#' \donttest{
#' S <- cbind(diag(2), c(-1, 1), c(1.1, -1))  # inequality matrix
#' # S'x >= 0 represents the wedge x1 <= x2 <= 1.1 x1
#' # example taken from Pakman and Paninski (2014)
#' # 1. exact Hamiltonian Monte Carlo (Pakman and Paninski, 2014)
#' sampler <- create_TMVN_sampler(Q=diag(2), mu=c(4, 4), S=S, method="HMC")
#' sim <- MCMCsim(sampler, n.iter=600, verbose=FALSE)
#' summary(sim)
#' plot(as.matrix(sim$x), pch=".")
#' # 2. exact Hamiltonian Monte Carlo with Laplace momentum (Nishimura et al., 2021)
#' sampler <- create_TMVN_sampler(Q=diag(2), mu=c(4, 4), S=S, method="HMCZigZag")
#' sim <- MCMCsim(sampler, n.iter=600, verbose=FALSE)
#' summary(sim)
#' plot(as.matrix(sim$x), pch=".")
#' # 3. Gibbs sampling approach (Rodriguez-Yam et al., 2004)
#' sampler <- create_TMVN_sampler(Q=diag(2), mu=c(4, 4), S=S, method="Gibbs")
#' sim <- MCMCsim(sampler, burnin=500, n.iter=2000, verbose=FALSE)
#' summary(sim)
#' plot(as.matrix(sim$x), pch=".")
#' # 4. soft TMVN approximation (Souris et al., 2018)
#' sampler <- create_TMVN_sampler(Q=diag(2), mu=c(4, 4), S=S, method="softTMVN")
#' sim <- MCMCsim(sampler, n.iter=600, verbose=FALSE)
#' summary(sim)
#' plot(as.matrix(sim$x), pch=".")
#' }
#'
#' @export
#' @author Harm Jan Boonstra, with help from Grzegorz Baltissen
#' @param Q precision matrix of the (unconstrained) multivariate normal distribution.
#' @param mu mean of the (unconstrained) multivariate normal distribution.
#' @param Xy alternative to specifying mu; in this case \code{mu} is computed as \eqn{Q^{-1}\code{Xy}}.
#' @param update.Q whether \code{Q} is updated for each draw. Currently only supported
#'  by methods 'direct' and 'HMC'.
#' @param update.mu whether \code{mu} is updated for each draw. By default equal to \code{update.Q}.
#'  Currently only supported by methods 'direct' and 'HMC'.
#' @param name name of the TMVN vector parameter.
#' @param coef.names optional labels for the components of the vector parameter.
#' @param R equality restriction matrix.
#' @param r rhs vector for equality constraints \eqn{R'x = r}, where \eqn{R'} denotes the transpose of R.
#' @param S inequality restriction matrix.
#' @param s rhs vector for inequality constraints \eqn{S'x >= s}, where \eqn{S'} denotes the transpose of S.
#' @param lower alternative to \code{s} for two-sided inequality restrictions \eqn{\code{lower} <= S'x <= \code{upper}}.
#' @param upper alternative to \code{s} for two-sided inequality restrictions \eqn{\code{lower} <= S'x <= \code{upper}}.
#' @param check.constraints if \code{TRUE} check whether the starting values satisfy all constraints.
#' @param method sampling method. The options are "direct" for direct sampling from the
#'  unconstrained or equality constrained multivariate normal (MVN). For inequality constrained
#'  MVN sampling three methods are supported: "HMC" for (exact) Hamiltonian Monte Carlo,
#'  "HMCZigZag" for (exact) Hamiltonian Monte Carlo with Laplace momentum, "Gibbs" for a
#'  component-wise Gibbs sampling approach, and "softTMVN" for a data augmentation method that samples
#'  from a smooth approximation to the truncated MVN. Alternatively, the method setting
#'  functions \code{m_direct}, \code{m_HMC}, \code{m_HMC_ZigZag}, \code{m_Gibbs} or
#'  \code{m_softTMVN} can be used to select the method and possibly set some of its
#'  options to non-default values, see \code{\link{TMVN-methods}}.
#' @param reduce whether to a priori restrict the simulation to the subspace defined by the
#'  equality constraints.
#' @param chol.control options for Cholesky decomposition, see \code{\link{chol_control}}.
#' @param debug if \code{TRUE} a breakpoint is set at the beginning of the TMVN sampling
#'  function \code{draw}.
#' @returns An environment for sampling from a possibly degenerate and truncated multivariate normal
#'  distribution.
#' @references
#'  Z.I. Botev and P. L'Ecuyer (2016).
#'    Simulation from the Normal Distribution Truncated to an Interval in the Tail.
#'    in VALUETOOLS.
#'
#'  Y. Cong, B. Chen and M. Zhou (2017).
#'    Fast simulation of hyperplane-truncated multivariate normal distributions.
#'    Bayesian Analysis 12(4), 1017-1037.
#'
#'  Y. Li and S.K. Ghosh (2015). Efficient sampling methods for truncated multivariate normal
#'    and student-t distributions subject to linear inequality constraints.
#'    Journal of Statistical Theory and Practice 9(4), 712-732.
#'
#'  A. Nishimura, Z. Zhang and M.A. Suchard (2021). Hamiltonian zigzag sampler got more momentum
#'    than its Markovian counterpart: Equivalence of two zigzags under a momentum refreshment limit.
#'    arXiv:2104.07694.
#'
#'  A. Pakman and L. Paninski (2014).
#'    Exact Hamiltonian Monte Carlo for truncated multivariate gaussians.
#'    Journal of Computational and Graphical Statistics 23(2), 518-542.
#'
#'  G. Rodriguez-Yam, R.A. Davis and L.L. Scharf (2004).
#'    Efficient Gibbs sampling of truncated multivariate normal with application to constrained linear regression.
#'    Unpublished manuscript.
#'
#'  H. Rue and L. Held (2005).
#'    Gaussian Markov Random Fields.
#'    Chapman & Hall/CRC.
#'
#'  A. Souris, A. Bhattacharya and P. Debdeep (2018).
#'    The Soft Multivariate Truncated Normal Distribution.
#'    arXiv:1807.09155.
#'
#'  K.A. Valeriano, C.E. Galarza and L.A. Matos (2023).
#'    Moments and random number generation for the truncated elliptical family of distributions.
#'    Statistics and Computing 33(1), 1-20.
create_TMVN_sampler <- function(Q, mu=NULL, Xy=NULL, update.Q=FALSE, update.mu=update.Q,
                                name="x", coef.names=NULL,
                                R=NULL, r=NULL, S=NULL, s=NULL, lower=NULL, upper=NULL,
                                check.constraints = FALSE,
                                method=NULL, reduce=NULL,
                                chol.control=chol_control(), debug=FALSE) {

  if (name == "") stop("empty name")
  store_default <- function(prior.sampler=FALSE) name  # for direct use of create_TMVN_sampler

  if (!update.Q)
    Q <- economizeMatrix(Q, symmetric=TRUE, drop.zeros=TRUE, check=TRUE)
  n <- nrow(Q)

  if (is.null(Xy)) {
    mu <- if (is.null(mu)) numeric(n) else as.numeric(mu)
    if (length(mu) != n) stop("incompatible dimensions 'Q' and 'mu'")
  } else {
    if (!is.null(mu)) stop("only one of 'mu' and 'Xy' should be specified")
  }

  if (!is.null(coef.names)) {
    if (length(coef.names) != n) stop("incompatible length of 'coef.names'")
    coef.names <- list(coef.names)
    names(coef.names) <- name
  }

  if (is.null(R)) {
    reduce <- FALSE
    eq <- FALSE
  } else {
    if (check.constraints) Rnames <- colnames(R)
    R <- economizeMatrix(R, vec.as.diag=FALSE, check=TRUE)
    if (nrow(R) != n || ncol(R) > n) stop("incompatible constraint matrix 'R'")
    if (is.null(r)) {
      r <- numeric(ncol(R))
    } else {
      if (length(r) != ncol(R)) stop("length of 'r' should equal the number of columns of 'R'")
    }
    # TODO check that R has full column rank
    if (ncol(R) > 0L) {
      if (check.constraints) {
        R0 <- R
        r0 <- r
        check_equalities <- function(x, tol=sqrt(.Machine$double.eps)) {
          tol <- abs(tol)
          dif <- crossprod_mv(R0, x) - r0
          viol <- which(abs(dif) > tol)
          if (length(viol) >= 1L) {
            warn(length(viol), " equality restriction(s) violated")
            message("largest discrepancies:")
            o <- viol[order(abs(dif[viol]), decreasing = TRUE)[1:min(3L, length(viol))]]
            fdif <- format(dif[o], digits=3L)
            message(if (is.null(Rnames)) fdif else paste0(fdif, " (", Rnames[o], ")", collapse=", "))
          }
        }
      }
      eq <- TRUE
    } else {
      eq <- FALSE
    }
  }
  if (!eq) rm(R, r)

  if (is.null(S)) {
    rm(S, s)
    ineq <- FALSE
  } else {
    if (check.constraints) {
      Snames <- colnames(S)
      if (is.null(Snames)) Snames <- seq_len(ncol(S))
    }
    if (!is.null(lower)) {
      # use lower and upper bounds, and translate to Sx >= s
      if (is.null(upper)) stop("lower and upper bound should be specified together")
      if (!is.null(s)) warn("argument 's' ignored in combination with 'lower' and 'upper'")
      if (length(lower) != length(upper) || any(lower >= upper)) stop("'lower' and 'upper' incompatible")
      if (length(lower) != ncol(S)) stop("'S' incompatible with 'lower' and 'upper'")
      if (any(is.infinite(lower)) || any(is.infinite(upper))) {
        S <- cbind(S[, is.finite(lower), drop=FALSE], -S[, is.finite(upper), drop=FALSE])
        s <- c(lower[is.finite(lower)], upper[is.finite(upper)])
        if (!length(s)) stop("no effective inequalities")
        if (check.constraints) Snames <- c(paste(Snames[is.finite(lower)], "lower"), paste(Snames[is.finite(upper)], "upper"))
      } else {
        S <- cbind(S, -S)
        s <- c(lower, -upper)
        if (check.constraints) Snames <- paste(c(Snames, Snames), rep_each(c("lower", "upper"), ncol(S)))
      }
    }
    S <- economizeMatrix(S, vec.as.diag=FALSE, check=TRUE)
    if (nrow(S) != n) stop("incompatible constraint matrix 'S'")
    ncS <- ncol(S)
    if (is.null(s)) {
      s <- numeric(ncS)
    } else {
      s <- as.numeric(s)
      if (length(s) != ncS) stop("length of 's' should equal the number of columns of 'S'")
    }
    if (length(s)) {
      if (check.constraints) {
        S0 <- S
        s0 <- s
        check_inequalities <- function(x, tol=sqrt(.Machine$double.eps)) {
          dif <- crossprod_mv(S0, x) - s0
          viol <- which(dif < -abs(tol))
          if (length(viol) >= 1L) {
            warn(length(viol), " inequality restriction(s) violated")
            message("largest discrepancies:")
            o <- viol[order(dif[viol])[seq_len(min(3L, length(viol)))]]
            fdif <- format(dif[o], digits=3L)
            message(if (is.null(Snames)) fdif else paste0(fdif, " (", Snames[o], ")", collapse=", "))
          }
        }
      }
      ineq <- TRUE
    } else {
      ineq <- FALSE
    }
  }
  rm(lower, upper)

  if (!eq && !ineq) check.constraints <- FALSE

  if (check.constraints) {
    check_constraints <- function(x, tol=sqrt(.Machine$double.eps)) {
      tol <- abs(tol)
      if (eq) check_equalities(x, tol)
      if (ineq) check_inequalities(x, tol)
    }
  }

  if (is.null(method)) method <- if (ineq) "HMC" else "direct"
  if (is.character(method)) {
    method <- match.arg(method, c("direct", "Gibbs", "HMC", "HMCZigZag", "softTMVN"))
    method <- eval(call(paste0("m_", method)))
  }
  if (ineq) {
    if (method$name == "direct") stop("method 'direct' cannot be used for inequality constrained sampling")
  } else {
    if (method$name == "softTMVN") stop("method 'softTMVN' can only be used for inequality constrained sampling")
  }
  
  if (is.null(reduce))
    reduce <- eq && method$name == "Gibbs"
  else if (eq && method$name == "Gibbs" && !reduce)
    stop("for method 'Gibbs' equality restrictions are only supported with 'reduce=TRUE'")
  if (update.Q || update.mu) {
    if (all(method$name != c("direct", "HMC"))) stop("'update.Q=TRUE' or 'update.mu=TRUE' only supported for methods 'direct' and 'HMC'")
    if (reduce) stop("'update.Q=TRUE' or 'update.mu=TRUE' not supported in combination with 'reduce=TRUE'")
  }
  if (reduce || (ineq && method$name == "Gibbs")) name.tr <- paste0(name, "_tr_")

  if (method$name == "direct" && is.null(method$use.cholV)) {
    method$use.cholV <- !reduce && !update.Q && n <= 1000L
  }

  if (reduce) {

    if (!is.null(Xy)) mu <- build_chol(Q, control=chol.control)$solve(Xy)

    # transform to subspace defined by R
    if (method$name == "softTMVN" || method$name == "HMCZigZag") {
      reduced.frame <- eliminate_equalities(R, r, Q, mu, keep.Q=TRUE, chol.control)
      Q <- reduced.frame$Q
    } else {
      reduced.frame <- eliminate_equalities(R, r, Q, mu, chol.control=chol.control)
      rm(Q)
    }
    cholQ <- reduced.frame$cholQ
    n.frame <- reduced.frame$n
    mu <- numeric(n.frame)  # mean is zero in the new frame

    if (ineq) {
      # transform s and S to the reduced frame
      s <- reduced.frame$transform_restriction_rhs(S, s)
      S <- reduced.frame$transform_restriction_matrix(S)
    }

  } else {  # !reduce

    if (isTRUE(method$use.cholV)) {
      tryCatch(
        suppressWarnings({
          V <- economizeMatrix(solve(Q), symmetric=TRUE)
          cholV <- economizeMatrix(chol(V))
        }),
        error = function(err) {
          # TODO using determinant is not a very good way to check for non-pd
          d <- determinant(Q, logarithm=FALSE)
          if (d$sign * d$modulus < .Machine$double.eps) {
            stop("Non-positive-definite matrix in model component `", name, "`. Consider using 'remove.redundant=TRUE' or increasing the coefficients' prior precision.")
          } else {
            stop(err)
          }
        }
      )
      # remove any cached Cholesky factorizations
      if (class(V)[1L] == "dsCMatrix") attr(V, "factors") <- list()
    } else {
      cholQ <- build_chol(Q, control=chol.control)
    }

    if (all(method$name != c("softTMVN", "HMCZigZag"))) rm(Q)

    if (!is.null(Xy)) {
      mu <- if (isTRUE(method$use.cholV)) V %m*v% Xy else cholQ$solve(Xy)
    }

    if (ineq) {
      # currently tabMatrix disallowed because used as solve rhs below
      S <- economizeMatrix(S, sparse=if (class(cholQ$cholM)[1L] == "matrix") FALSE else NULL, allow.tabMatrix=FALSE)
    }
    if (method$name == "Gibbs" || method$name == "HMCZigZag") n.frame <- n

    if (eq && any(method$name == c("direct", "HMC"))) {
      # conditioning by Kriging to deal with equality constraints R'x = r
      # update.Q, dense cholQ --> faster to have dense R (and S) as well
      R <- economizeMatrix(R, sparse=if (update.Q && class(cholQ$cholM)[1L] == "matrix") FALSE else NULL, allow.tabMatrix=FALSE)
      if (isTRUE(method$use.cholV))
        projector <- create_projector(V=V, R=R)
      else
        projector <- create_projector(cholQ=cholQ, update.Q=update.Q, R=R)
      if (all(r == 0)) r <- NULL
    }  # END if (eq)
  }  # END !reduce
  if (isTRUE(method$use.cholV) && !update.mu) rm(V)
  rm(Xy)

  self <- environment()


  zero.mu <- !update.mu && all(mu == 0)
  if (any(method$name == c("HMC", "HMCZigZag"))) {
    if (eq && !reduce && any(r != 0)) zero.mu <- FALSE
  }
  if (zero.mu) mu <- numeric(0L)

  if (method$name == "Gibbs") {
    if (eq) {
      transform <- function(x) cholQ$Ltimes(reduced.frame$transform(x))
      untransform <- function(z) reduced.frame$untransform(cholQ$solve(z, system="Lt"))
    } else if (zero.mu) {
      transform <- function(x) cholQ$Ltimes(x)
      untransform <- function(z) cholQ$solve(z, system="Lt")
    } else {
      transform <- function(x) cholQ$Ltimes(x - mu)
      untransform <- function(z) mu + cholQ$solve(z, system="Lt")
    }
  } else if (reduce) {
    transform <- reduced.frame$transform
    untransform <- reduced.frame$untransform
  }


  # BEGIN draw function
  draw <- if (debug) function(p) {browser()} else function(p) {}
  if (method$name != "direct" && !(method$name == "Gibbs" && !ineq)) {
    if (reduce || (method$name == "Gibbs"))
      draw <- add(draw, bquote(x <- p[[.(name.tr)]]))
    else
      draw <- add(draw, bquote(x <- p[[.(name)]]))
  }

  if (method$name == "direct") {
    if (update.Q) draw <- add(draw, quote(cholQ$update(Q, Imult)))
    if (zero.mu)
      if (method$use.cholV)
        draw <- add(draw, bquote(x <- drawMVN_cholV(.(n), cholV, scale)))
      else
        draw <- add(draw, quote(x <- drawMVN_cholQ(cholQ, sd=scale)))
    else {
      if (update.mu)
        if (method$use.cholV)
          draw <- add(draw, bquote(x <- V %m*v% Xy + drawMVN_cholV(.(n), cholV, scale)))
        else
          draw <- add(draw, quote(x <- drawMVN_cholQ(cholQ, Xy, sd=scale)))
      else
        if (method$use.cholV)
          draw <- add(draw, bquote(x <- mu + drawMVN_cholV(.(n), cholV, scale)))
        else
          draw <- add(draw, quote(x <- mu + drawMVN_cholQ(cholQ, sd=scale)))
    }
    if (eq && !reduce) {
      if (update.Q)
        draw <- add(draw, quote(projector$signal_cholQ_change()))
      draw <- add(draw, quote(x <- projector$project(x, cholQ, r=r)))
    }
  }  # END direct

  if (method$name == "Gibbs") {
    if (ineq) {
      # transform to zero mean unit covariance matrix frame
      # inequalities U'v >= u
      u <- if (zero.mu) s else s - crossprod_mv(S, mu)
      # drop very small numbers in U as they give problems in the Gibbs sampler
      # and use transpose for faster col access in sparse case
      Ut <- economizeMatrix(
        t(drop0(cholQ$solve(S, system="L"), tol=1e-10)),
        vec.diag=FALSE, allow.tabMatrix=FALSE
      )
      if (all(class(Ut)[1L] != c("matrix", "dgCMatrix"))) Ut <- as(as(Ut, "CsparseMatrix"), "generalMatrix")
      
      eps <- method$eps
      # draw from truncated univariate normal full conditionals
      draw <- add(draw, quote(ustar <- u - Ut %m*v% x))
      if (method$slice) {
        if (is.matrix(Ut))
          draw <- add(draw, quote(x <- Crtmvn_slice_Gibbs_dense(x, Ut, ustar, eps)))
        else
          draw <- add(draw, quote(x <- Crtmvn_slice_Gibbs_sparse(x, Ut, ustar, eps)))
      } else {
        if (is.matrix(Ut))
          draw <- add(draw, quote(x <- Crtmvn_Gibbs_dense(x, Ut, ustar, eps)))
        else
          draw <- add(draw, quote(x <- Crtmvn_Gibbs_sparse(x, Ut, ustar, eps)))
      }
    } else {
      draw <- add(draw, bquote(x <- Crnorm(.(n.frame))))
    }
  }  # END Gibbs

  if (method$name == "softTMVN") {
    useV <- method$useV
    if (!is.null(method$CG) && eq) stop("conjugate gradients with equality constraints is not supported")
    if (method$PG.approx) {
      mPG <- as.integer(method$PG.approx.m)
      if (all(length(mPG) != c(1L, n))) stop("invalid value for option 'PG.approx.m'")
      rPolyaGamma <- function(b, c) CrPGapprox(ncS, b, c, mPG)
    } else {
      if (!requireNamespace("BayesLogit", quietly=TRUE)) stop("please install package 'BayesLogit' and try again")
      rpg <- BayesLogit::rpg
      rPolyaGamma <- function(b, c) rpg(ncS, b, c)
    }
    if (useV) {
      # 'dual' version of MVN sampling
      rm(Q)
      V <- cholQ$inverse()
      # define Diagonal matrix template to hold omega.tilde for use in sparse matrix templated sum
      D.omega.tilde.inv <- Cdiag(10*runif(ncS, 0.9, 1.1))  # choose large enough so that the matrix template below is pos-def
      matsum <- make_mat_sum(M0=economizeMatrix(crossprod_sym(S, V), symmetric=TRUE, drop.zeros=TRUE), M1=D.omega.tilde.inv)
      ch <- build_chol(matsum(D.omega.tilde.inv), control=chol.control)
      VS <- economizeMatrix(V %*% S, drop.zeros=TRUE)
      if (eq && !reduce) {
        VR <- economizeMatrix(V %*% R, drop.zeros=TRUE)
        RVR <- economizeMatrix(crossprod(R, VR), drop.zeros=TRUE, symmetric=TRUE)
        SVR <- economizeMatrix(crossprod(S, VR), drop.zeros=TRUE)
        temp <- crossprod_sym2(SVR, ch$solve(SVR))
        matsum_RVxR <- make_mat_sum(M0 = RVR, M1 = -temp)
        cholRVxR <- build_chol(matsum_RVxR(M1 = temp, w1 = -1), control=chol.control)
        rm(temp)
      }
    } else {
      St <- t(S)  # for use in crossprod_sym
      if (!zero.mu) Qmu <- Q %m*v% mu
      matsum <- make_mat_sum(
        M0 = Q, M1 = crossprod_sym(St, runif(ncS, 0.9, 1.1)),
        force.sparse = eq && inherits(R, "CsparseMatrix")
      )
      ch <- build_chol(matsum(crossprod_sym(St, runif(ncS, 0.9, 1.1))), control=chol.control)
      if (eq && !reduce)
        projector <- create_projector(cholQ=ch, update.Q=TRUE, R=R)
    }

    sharpness <- method$sharpness
    sharpness_squared <- sharpness^2
    sharpness_inv <- 1/sharpness
    # 1. draw Polya-Gamma mixture precisions
    draw <- draw |>
      add(quote(discr <- crossprod_mv(S, x) - s)) |>
      add(quote(omega <- rPolyaGamma(1, sharpness * discr))) |>
      # 2. draw vector of coefficients
      add(quote(omega.tilde <- sharpness_squared * omega))
    if (useV) {
      draw <- draw |>
        add(quote(attr(D.omega.tilde.inv, "x") <- 1 / omega.tilde)) |>
        add(quote(XX_V <- matsum(D.omega.tilde.inv))) |>
        add(quote(ch$update(XX_V))) |>
        add(quote(alpha <- sharpness * s + 0.5/omega))
      if (zero.mu)
        draw <- add(draw, quote(y1 <- drawMVN_cholQ(cholQ)))
      else
        draw <- add(draw, quote(y1 <- mu + drawMVN_cholQ(cholQ)))
      draw <- add(draw, bquote(y2 <- Crnorm(.(length(s))) * sqrt(1/omega)))
      draw <- add(draw, quote(x <- y1 + VS %m*v% ch$solve(sharpness_inv * (alpha - y2) - crossprod_mv(S, y1))))
      if (eq && !reduce) {
        draw <- draw |>
          add(quote(cholRVxR$update(matsum_RVxR(M1=crossprod_sym2(SVR, ch$solve(SVR)), w1=-1)))) |>
          add(quote(temp <- V %m*v% (R %m*v% cholRVxR$solve(r - crossprod_mv(R, x))))) |>
          add(quote(x <- x + temp - V %m*v% (S %m*v% ch$solve(crossprod_mv(S, temp)))))
      }
    } else {
      draw <- add(draw, quote(alpha <- 0.5 * sharpness + s * omega.tilde))
      if (zero.mu)
        draw <- add(draw, quote(Xy <- crossprod_mv(St, alpha)))
      else
        draw <- add(draw, quote(Xy <- crossprod_mv(St, alpha) + Qmu))
      if (is.null(method$CG)) {
        # Cholesky
        draw <- draw |>
          add(quote(XX <- crossprod_sym(St, omega.tilde))) |>
          add(quote(XX_Q <- matsum(XX))) |>
          add(quote(ch$update(XX_Q)))
        draw <- add(draw, quote(x <- drawMVN_cholQ(ch, Xy)))
        if (eq && !reduce) {
          draw <- draw |>
            add(quote(projector$signal_cholQ_change())) |>
            add(quote(x <- projector$project(x, ch, r)))
        }
      } else {
        # conjugate gradients
        draw <- add(draw, bquote(u <- Xy + crossprod_mv(St, sqrt(omega.tilde) * Crnorm(.(ncS))) + cholQ$Ltimes(Crnorm(.(n)), transpose=FALSE)))
        A_times <- function(x, omega.tilde) crossprod_mv(St, omega.tilde * (St %m*v% x)) + Q %m*v% x
        dQinv <- sqrt(1/diag(Q))  # seems better as a preconditioner than 1/diag(Q)
        M_solve <- function(x, omega.tilde) dQinv * x
        max.it <- method$CG$max.it
        stop.criterion <- method$CG$stop.criterion
        verbose <- method$CG$verbose
        draw <- add(draw, quote(x <- CG(u, self, x, max.it=max.it,
            e=stop.criterion, verbose=verbose, omega.tilde=omega.tilde)
        ))
      }
    }
  }  # END softTMVN

  if (method$name == "HMC") {
    # Hamiltonian H = 1/2 x'Mx - g'x + 1/2 p'M^{-1}p
    # precision matrix M, covariance matrix M^{-1}, mu=M^{-1}g must obey equality restrictions
    # Hamilton's equations: dx/dt = M^{-1}p, p = M dx/dt
    #                       dp/dt = - Mx + g
    if (eq && !reduce && !update.Q) {
      # project mean mu on eq constraint surface
      if (!zero.mu) mu <- projector$project(mu, cholQ, r=r)
    }
    if (ineq) {
      if (reduce) {
        # define VS as Q^-1 S = V S
        VS <- economizeMatrix(cholQ$solve(S), allow.tabMatrix=FALSE)
      } else if (!update.Q) {
        # define VS as Q^-1 S = V S; alternatively transform after projection
        VS <- economizeMatrix(cholQ$solve(S), sparse=if (is.matrix(cholQ$cholM)) FALSE else NULL, allow.tabMatrix=FALSE)
        if (eq) {
          # project inequality matrix on equality constraint surface, use Q as a metric
          VS <- economizeMatrix(projector$project(VS, cholQ), sparse=if (is.matrix(cholQ$cholM)) FALSE else NULL, allow.tabMatrix=FALSE)
        }
      } else {
        VS <- NULL
      }
      if (class(S)[1L] == "ddiMatrix") S <- as(as(S, "CsparseMatrix"), "generalMatrix")
      if (!is.null(VS) && class(VS)[1L] == "ddiMatrix") VS <- as(as(VS, "CsparseMatrix"), "generalMatrix")
      if (update.Q) {
        simplified <- FALSE
      } else {
        refl.fac <- 2 / colSums(S * VS)  # 2 / vector of normal' Q normal for all inequalities
        s.adj <- if (zero.mu) s else s - crossprod_mv(S, mu)  # if positive, mu violates inequalities
        abs.s.adj <- abs(s.adj)
        simplified <- !eq && identical(VS, S)  # simplified=TRUE is the case described in Pakman and Paninski: identity Q and no equality constraints
        if (simplified) VS <- NULL
      }
      max.events <- method$max.events
      diagnostic <- method$diagnostic
      if (diagnostic) {
        # set up a vector of wall bounce counts
        bounces <- setNames(integer(ncS), if (is.null(colnames(S))) seq_len(ncS) else colnames(S))
        display.n <- seq_len(min(10L, ncS))  # number of counts to display
      } else {
        bounces <- integer(0L)
      }
    }

    if (update.Q) {
      draw <- draw |>
        add(quote(cholQ$update(Q, Imult))) |>
        add(VS <- cholQ$solve(S))
      if (eq) {
        draw <- add(draw, quote(projector$signal_cholQ_change()))
        draw <- add(draw, quote(VS <- projector$project(VS, cholQ, r)))
      }
      draw <- add(draw, quote(refl.fac <- 2 / colSums(S * VS)))  # 2 / vector of normal' Q normal for all inequalities
    }
    if (update.mu) {
      draw <- add(draw, quote(mu <- cholQ$solve(Xy)))
      if (eq)
        draw <- add(draw, quote(mu <- projector$project(mu, cholQ, r)))
      draw <- draw |>
        add(quote(s.adj <- s - crossprod_mv(S, mu))) |>  # if positive, mu violates inequalities
        add(quote(abs.s.adj <- abs(s.adj)))
    }

    # 1. draw p from N(0, M^{-1}); instead draw v=dx/dt from N(0, M)
    #    this is the initial velocity (starting from final position x of previous draw)
    draw <- add(draw, quote(v <- drawMVN_cholQ(cholQ, sd=scale)))  # draw velocity at start time t=0
    if (!reduce) {
      if (eq) {
        # project velocity direction along constraint surface
        draw <- add(draw, quote(v <- projector$project(v, cholQ)))
      }
    }
    draw <- add(draw, quote(T0 <- method$Tsim()))
    # solution: x(t) = mu + v0*sin(t) + (x0 - mu)*cos(t)
    # inequalities: S'mu + S'v sin(t) + S'(x0 - mu) cos(t) >= s
    #           --> S'mu + u cos(t + phi) >= s
    if (ineq) {
      if (diagnostic) {
        draw <- draw |>
          add(quote(viol <- crossprod_mv(S, x) < s)) |>
          add(quote(if (any(viol)) cat("\nviolated constraints:", names(bounces)[viol])))
      }
      # NB bounces is updated by reference
      draw <- add(draw, quote(x <- TMVN_HMC_C(S, ncS, v, x, s.adj, refl.fac, zero.mu, mu, simplified, VS, diagnostic, bounces, T0, max.events)))
      if (diagnostic) {
        draw <- draw |>
          add(bquote(most_bounces <- sort(bounces, decreasing=TRUE)[.(display.n)])) |>
          add(quote(cat("\nmost bounces:", paste0(names(most_bounces), " (", most_bounces, ") "))))
      }
    } else {
      if (zero.mu)
        draw <- add(draw, quote(x <- sin(T0) * v + cos(T0) * x))
      else
        draw <- add(draw, quote(x <- mu + sin(T0) * v + cos(T0) * (x - mu)))
    }
  }  # END HMC

  if (method$name == "HMCZigZag") {
    # TODO: update.Q and update.mu cases

    base_which <- base::which  # faster, in case Matrix::which is loaded
    eps <- sqrt(.Machine$double.eps)
    negeps <- -eps

    if (class(Q)[1L] != "matrix") {
      iQ <- xQ <- list()
      for (j in seq_len(ncol(Q))) {
        temp <- Q[, j, drop=TRUE]
        iQ[[j]] <- base_which(temp != 0)
        xQ[[j]] <- temp[temp != 0]
      }
    }

    if (ineq) {
      if (class(S)[1L] == "ddiMatrix") S <- as(as(S, "CsparseMatrix"), "generalMatrix")
      if (class(S)[1L] == "dgCMatrix") {
        # store nonzero indices for each column
        inds <- list()
        len <- diff(S@p)
        for (j in seq_len(ncS)) {
          inds[[j]] <- S@i[(S@p[j] + 1L):(S@p[j] + len[j])] + 1L
        }
      }
    }

    if (eq && !reduce) {
      # TODO can we use the approximate precision matrix form here using matrix-inversion lemma?
      #      then we would nowhere need to factorize Q
      # project mean mu on eq constraint surface
      VR <- cholQ$solve(R)
      if (!zero.mu) {
        VR.RVRinv <- economizeMatrix(VR %*% solve(crossprod_sym2(R, VR)), drop.zeros=TRUE)
        mu <- mu + VR.RVRinv %m*v% (r - crossprod_mv(R, mu))
        rm(VR.RVRinv)
      }
      if (is.null(method$prec.eq))
        prec.eq <- 100 * diag(solve(crossprod_sym2(R, VR)))
      else
        prec.eq <- method$prec.eq
      rm(VR)
    }

    adapt <- method$adapt
    rate <- method$rate
    # NB if (reduce) then length of rate should be either 1 or the dimension of the reduced frame
    if (length(rate) == 1L) rate <- rep.int(rate, n.frame)
    if (length(rate) != n.frame) stop("Laplace distribution rate parameter has wrong length")
    unit.rate <- all(rate == 1) && !adapt
    diagnostic <- method$diagnostic
    if (diagnostic) {
      gradbounces <- setNames(integer(n), seq_len(n))
      display.n <- seq_len(min(10L, n))  # number of counts to display
      if (ineq) {
        # set up a vector of wall bounce counts
        bounces <- setNames(integer(ncS), if (is.null(colnames(S))) seq_len(ncS) else colnames(S))
        display.ncS <- seq_len(min(10L, ncS))  # number of counts to display
        flips <- setNames(integer(n), seq_len(n))
      }
    }
    if (adapt) {
      if (!ineq || !diagnostic) stop("adapt only in combination with diagnostic=TRUE")
    }

    draw <- draw |>
      add(quote(T0 <- method$Tsim())) |>
      add(quote(Tleft <- T0))
    if (diagnostic && ineq) {
      draw <- draw |>
        add(quote(viol <- crossprod_mv(S, x) < s)) |>
        add(quote(if (any(viol)) cat("\nviolated inequality constraints:", names(bounces)[viol])))  # print violated constraints
    }
    # draw pM from Laplace distribution, with parameters w; use pM for momentum to distinguish from state p
    if (unit.rate) {
      # this seems a fast way to generate Laplace variates:
      draw <- draw |>
        add(bquote(pM <- rexp(.(n.frame)) - rexp(.(n.frame)))) |>
        add(quote(v <- sign(pM)))
    } else {
      draw <- draw |>
        add(bquote(pM <- rexp(.(n.frame), rate=rate) - rexp(.(n.frame), rate=rate))) |>
        add(quote(v <- rate * sign(pM)))
    }
    draw <- draw |>
      add(bquote(Qx <- .(if (zero.mu) quote(Q %m*v% x) else quote(Q %m*v% (x - mu))))) |>
      add(quote(Qv <- Q %m*v% v))
    if (eq && !reduce) {
      if (zero.mu)
        draw <- add(draw, quote(Qx <- Qx + R %m*v% (prec.eq * crossprod_mv(R, x))))
      else
        draw <- add(draw, quote(Qx <- Qx + R %m*v% (prec.eq * crossprod_mv(R, x - mu))))
      draw <- add(draw, quote(Qv <- Qv + R %m*v% (prec.eq * crossprod_mv(R, v))))
    }

    draw <- add(draw, quote(event <- 0L))
    draw <- add(draw, quote(
      repeat {
        # compute gradient event time dt.gr, and possibly boundary event times
        discr <- Qx*Qx + 2*Qv*pM
        iposD <- base_which(discr > 0)
        dt.gr <- Inf
        gradient.event <- TRUE
        if (length(iposD)) {
          sqrtD <- sqrt(discr[iposD])
          t1 <- (-Qx[iposD] - sqrtD) / Qv[iposD]
          ind1 <- base_which(t1 > 0)
          t2 <- (-Qx[iposD] + sqrtD) / Qv[iposD]
          ind2 <- base_which(t2 > 0)
          if (length(ind1) && length(ind2)) {
            i1 <- ind1[which.min(t1[ind1])]
            i2 <- ind2[which.min(t2[ind2])]
            dt.gr <- min(t1[i1], t2[i2])
            istar <- if (t1[i1] < t2[i2]) iposD[i1] else iposD[i2]
          } else {
            if (length(ind1)) {
              i1 <- ind1[which.min(t1[ind1])]
              dt.gr <- t1[i1]
              istar <- iposD[i1]
            } else if (length(ind2)) {
              i2 <- ind2[which.min(t2[ind2])]
              dt.gr <- t2[i2]
              istar <- iposD[i2]
            }
          }
        }

        dt.1st <- dt.gr

        if (ineq) {
          # boundary events
          # TODO can we update Sv and Sx more efficiently then simply recompute them?
          Sv <- crossprod_mv(S, v)
          dtS <- (s - crossprod_mv(S, x)) / Sv

          # ignore inequality boundaries at which the particle currently is (including previously hit walls) while moving to the interior
          # use negeps here, as small round-off like violations of inequalities are bound to occur
          iS <- base_which(dtS > negeps & Sv < 0)  # allow passing to the right side
          if (length(iS)) {
            jstar <- iS[which.min(dtS[iS])]
            dt.bd <- dtS[jstar]  # first boundary event time
            if (dt.bd < dt.gr) {
              gradient.event <- FALSE
              dt.1st <- dt.bd
            }
          }
        }

        if (dt.1st >= Tleft) {
          x <- x + Tleft * v
          # no need to update pM, as it will be refreshed in the next iteration
          break
        }

        # update x and pM
        x <- x + dt.1st * v
        pM <- pM - dt.1st * Qx - 0.5 * dt.1st * dt.1st * Qv
        # pM[istar] will be 0 after a gradient event,
        # but for numerical robustness (important!!) set it to 0 below
        # update Qx
        if (!eq) Qx <- Qx + dt.1st * Qv
        # all.equal(Qx, crossprod_mv(Q, x - mu))
        # update v
        if (gradient.event) {
          pM[istar] <- 0  # for numerical robustness (important!)
          v[istar] <- -v[istar]
          if (!eq || reduce) {
            # update Qv (NB v[istar] has already changed sign)
            #Qv <- Qv + 2 * v[istar] * get_col(Q, istar)
            if (class(Q)[1L] == "matrix") {
              Qv <- Qv + 2 * v[istar] * Q[, istar, drop=TRUE]
            } else {
              Qv[iQ[[istar]]] <- Qv[iQ[[istar]]] + 2 * v[istar] * xQ[[istar]]
            }
            # all.equal(Qv, crossprod_mv(Q, v))
          }
          if (diagnostic) {
            gradbounces[istar] <- gradbounces[istar] + 1L
            if (ineq) flips[istar] <- flips[istar] + 1L
          }
        } else {
          # group of L1 isometries: all signed permutations
          # use a simple version: flip all pi's with si != 0
          # this guarantees that sum(Sj * v) > 0 after the bounce
          # TODO check whether other options e.g. including permutations are valid
          if (class(S)[1L] == "matrix")
            ind <- base_which(S[, jstar, drop=TRUE] != 0)
          else
            ind <- inds[[jstar]]
          pM[ind] <- -pM[ind]
          v[ind] <- -v[ind]
          # update Qv
          if (!eq || reduce) Qv <- crossprod_mv(Q, v)  # TODO update using only changed v components
          if (diagnostic) {
            bounces[jstar] <- bounces[jstar] + 1L
            flips[ind] <- flips[ind] + 1L
          }
        }

        if (eq && !reduce) {
          # recompute Qx and Qv
          # TODO make this more efficient?
          Qx <- if (zero.mu) Q %m*v% x else Q %m*v% (x - mu)
          if (zero.mu)
            Qx <- Qx + R %m*v% (prec.eq * crossprod_mv(R, x))
          else
            Qx <- Qx + R %m*v% (prec.eq * crossprod_mv(R, x - mu))
          Qv <- Q %m*v% v + R %m*v% (prec.eq * crossprod_mv(R, v))
        }

        Tleft <- Tleft - dt.1st
      }  # END repeat
    ))

    if (diagnostic) {
      if (ineq) {
        draw <- draw |>
          add(bquote(most_bounces <- sort(flips, decreasing=TRUE)[.(display.n)])) |>
          add(quote(cat("\nmost flips:", paste0(names(most_bounces), " (", most_bounces, ") ")))) |>
          add(bquote(most_bounces <- sort(bounces, decreasing=TRUE)[.(display.ncS)])) |>
          add(quote(cat("\nmost bounces:", paste0(names(most_bounces), " (", most_bounces, ") "))))
      }
      draw <- draw |>
        add(bquote(most_bounces <- sort(gradbounces, decreasing=TRUE)[.(display.n)])) |>
        add(bquote(cat("\nmost grad bounces:", paste0(names(most_bounces), " (", most_bounces, ") "))))
      if (adapt) {
        # TODO account for the number of inequalities each variable is involved in
        draw <- add(draw, quote(rate[flips > 100 * T0] <<- rate[flips > 100 * T0] / 1.1))
      }
    }
  }  # END HMCZigZag

  if (reduce || method$name == "Gibbs") {
    if (method$name != "direct" && !(method$name == "Gibbs" && !ineq))
      draw <- add(draw, bquote(p[[.(name.tr)]] <- x))
    draw <- add(draw, bquote(p[[.(name)]] <- untransform(x)))
  } else {
    draw <- add(draw, bquote(p[[.(name)]] <- x))
  }
  draw <- add(draw, quote(p))
  # set function signature (avoiding check NOTE)
  if (update.Q) {
    # total precision matrix to use in chol is Q + Imult*I
    formals(draw) <- alist(p=, scale=1, Q=, Imult=0, Xy=)
  } else {
    if (update.mu)
      formals(draw) <- alist(p=, scale=1, Xy=)
    else
      formals(draw) <- alist(p=, scale=1)
  }
  # END draw function

  if (ineq && method$name == "Gibbs") {
    start <- function(p=list(), scale=1) {
      if (!is.null(p[[name.tr]])) {
        if (length(p[[name.tr]]) != n.frame) stop("wrong length of start value for '", name.tr, "'")
        p[[name]] <- untransform(p[[name.tr]])
      } else if (!is.null(p[[name]])) {
        if (length(p[[name]]) != n) stop("wrong length of start value for '", name, "'")
        p[[name.tr]] <- transform(p[[name]])
      } else {
        if (!requireNamespace("lintools", quietly=TRUE)) stop("please install package lintools and try again")
        temp <- -as(as(Ut, "TsparseMatrix"), "generalMatrix")
        # sparse version
        A <- data.frame(
          row = temp@i,  # constraint index
          col = temp@j,  # variable index
          coef = temp@x
        )
        x0 <- Crnorm(n.frame, sd=scale)
        for (i in 1:10) {  # at most 10 trials
          epsilon <- scale * rexp(ncS)
          res <- lintools::sparse_project(x=x0, A, b=-(u + epsilon), neq=0L, base=0L, sorted=FALSE)
          scale <- 0.3 * scale
          if (anyNA(res$x)) {
            x0 <- Crnorm(n.frame, sd=scale)
          } else {
            if (all(Ut %m*v% res$x >= u)) break
            x0 <- res$x + Crnorm(n.frame, sd=scale)
          }
        }
        p[[name.tr]] <- res$x
        p[[name]] <- untransform(p[[name.tr]])
      }
      if (check.constraints) check_constraints(p[[name]])
      p
    }
  } else {
    # TODO here too use lintools to find a better starting value
    if (reduce) {
      start <- function(p=list(), scale=1) {
        if (!is.null(p[[name.tr]])) {
          return(p)
        } else {
          if (!is.null(p[[name]])) {
            p[[name.tr]] <- transform(p[[name]])
            return(p)
          }
        }
      }
    } else {
      start <- function(p=list(), scale=1) {
        if (!is.null(p[[name]])) {
          p[[name]] <- as.numeric(p[[name]])
          return(p)
        }
      }
    }
    if (isTRUE(method$use.cholV)) {
      if (zero.mu) {
        start <- add(start, bquote(x <- drawMVN_cholV(.(n), cholV, scale)))
      } else {
        start <- add(start, bquote(x <- mu + drawMVN_cholV(.(n), cholV, scale)))
      }
    } else {
      if (zero.mu || update.Q) {
        start <- add(start, quote(x <- drawMVN_cholQ(cholQ, sd=scale)))
      } else {
        start <- add(start, quote(x <- mu + drawMVN_cholQ(cholQ, sd=scale)))
      }
    }
    if (eq && !reduce) {
      if (method$name == "softTMVN") {
        if (useV) {
          start <- start |>
            add(quote(cholRVxR$update(matsum_RVxR(M1=crossprod_sym2(SVR, ch$solve(SVR)), w1=-1)))) |>
            add(quote(temp <- V %m*v% (R %m*v% cholRVxR$solve(r - crossprod_mv(R, x))))) |>
            add(bquote(p[[.(name)]] <- x + temp - V %m*v% (S %m*v% ch$solve(crossprod_mv(S, temp)))))
        } else
          start <- add(start, bquote(p[[.(name)]] <- projector$project(x, ch, r)))
      } else {
        if (method$name == "HMCZigZag")
          start <- add(start, bquote(p[[.(name)]] <- x))  # maybe also create and use projector?
        else
          start <- add(start, bquote(p[[.(name)]] <- projector$project(x, cholQ, r)))
      }
    } else {
      if (reduce) {
        start <- start |>
          add(bquote(p[[.(name.tr)]] <- x)) |>
          add(bquote(p[[.(name)]] <- untransform(x)))
      } else {
        start <- add(start, bquote(p[[.(name)]] <- x))
      }
    }
    if (check.constraints) start <- add(start, bquote(check_constraints(p[[.(name)]])))
    start <- add(start, quote(p))
  }

  rm(chol.control)

  self
}
