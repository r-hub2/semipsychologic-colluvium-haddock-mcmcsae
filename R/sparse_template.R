# mc: model component
# update.XX: whether XX is updated in each draw
# control: sampler options
# this function is called for its side effect of assigning a MVN sampler object
#   to environment mc, as well as an update method and possibly mat_sum and kron_prod methods
sparse_template <- function(mc, update.XX=FALSE, control=NULL) {

  if (mc[["type"]] == "block") {
    Q <- mc[["QT"]]
    keep.kp <- FALSE
  } else {
    # build a Kronecker product template
    Qv <- mc$generate_Qv()
    if (mc[["gl"]]) {
      QA <- mc$glp[["QA.ext"]]
    } else if (mc$strucA[["update.Q"]]) {
      QA <- mc$strucA$update_Q(mc[["QA"]], runif(1L, 0.25, 0.75))
    } else if (is.null(mc[["priorA"]])) {
      QA <- if (is.null(mc[["AR1.inferred"]])) mc[["QA"]] else mc$QA.template$QA0.5
    } else {
      if (is.null(mc[["AR1.inferred"]]))
        QA <- crossprod_sym(mc[["DA"]], runif(mc[["lD"]], 0.75, 1.25))
      else
        QA <- crossprod_sym(mc$DA.template[["DA0.5"]], runif(mc[["lD"]], 0.75, 1.25))
    }

    kron_prod <- build_kron(QA, Qv, q2=mc[["q0"]],
      M1.fixed = is.null(mc[["priorA"]]) && !mc$strucA[["update.Q"]] && is.null(mc[["AR1.inferred"]])
    )
    Q <- kron_prod(QA, Qv)
    if (mc[["in_block"]]) {
      if (is.matrix(Q)) Q <- as(as(Q, "CsparseMatrix"), "symmetricMatrix")
      assign("Q", Q, envir=mc)  # only for single-time use in create_mc_block
    }
    keep.kp <- TRUE
  }

  if (mc[["type"]] == "block" || !mc[["in_block"]]) {
    if (mc[["type"]] == "gen" && mc[["gl"]]) {
      XX <- mc$glp[["XX.ext"]]
      R <- mc$glp[["R"]]
    } else {
      XX <- mc[["XX"]]
      R <- mc[["R"]]
    }

    add.eps.I <- if (is.null(control)) FALSE else control[["add.eps.I"]]
    add.outer.R <- if (add.eps.I || is.null(control)) FALSE else control[["add.outer.R"]]
    cMVN.sampler <- if (is.null(control)) FALSE else control[["cMVN.sampler"]]
    if (cMVN.sampler) {
      if (update.XX) {
        mat_sum <- make_mat_sum(M1=XX, M2=Q, force.sparse=TRUE)
        XX_Q <- mat_sum(XX, Q)
      } else {
        mat_sum <- make_mat_sum(M0=XX, M1=Q, force.sparse=TRUE)
        XX_Q <- mat_sum(Q)
      }
      mc$MVNsampler <- create_block_cMVN_sampler(
        mbs=mc[["mcs"]], X=mc[["X"]], Q=XX_Q,
        R=R, r=mc[["r"]],
        sampler=mc[["e"]],
        name=mc[["name"]], chol.control=mc[["e"]]$control[["chol.control"]]
      )
    } else {
      mc$MVNsampler <- NULL
      R0 <- R
      if (is.null(add.outer.R)) {
        # unresolved singularities usually only occur when at least two components are improper
        add.outer.R <- mc[["type"]] == "block" && sum(!b_apply(mc[["mcs"]], `[[`, "is.proper")) >= 2L
        if (add.outer.R && !is.null(R)) {
          # it should not always be necessary to use the full matrix R
          # leave out the largest of the individual model components' contribution:
          cols <- i_apply(mc[["mcs"]], \(x) if (is.null(x[["R"]])) 0L else ncol(x[["R"]]))
          cols <- cols[cols != 0L]
          if (length(cols) > 1L) {
            inds <- c(0L, cumsum(cols))
            maxR <- which.max(cols)
            R0 <- R[, -((inds[maxR] + 1L):inds[maxR + 1L]), drop=FALSE]
          }
        }
      }

      # possibly affecting constraints: IGMRF factors, user-defined constraints, bym2, gl
      constraints <- set_constraints(R=R, r=mc[["r"]], S=mc[["S"]], s=mc[["s"]])

      if (!add.outer.R || is.null(R)) {
        if (update.XX) {
          mat_sum <- make_mat_sum(M1=XX, M2=Q)
          XX_Q <- mat_sum(XX, Q)
        } else {
          mat_sum <- make_mat_sum(M0=XX, M1=Q)
          XX_Q <- mat_sum(Q)
        }
        tryCatch(
          suppressWarnings(
            mc$MVNsampler <- create_TMVN_sampler(
              Q=XX_Q, update.Q=TRUE, name=mc[["name"]],
              constraints=constraints,
              chol.control=mc[["e"]]$control[["chol.control"]]
            )
          ),
          error = function(e) {
            #if (grepl("Cholesky", e$message)) {  # TODO are there other possible messages for singular matrix conditions?
            # the problem is most likely a non-positive-definite XX + Q because of unresolved
            #   singularities in XX (due to levels without observations, or blocking and noninformative regression prior)
            # try to correct this by adding a multiple of tcrossprod(R) to XX (Rue and Held, 2005)
            #   the effect of which is cancelled by imposing these restrictions (R) during sampling
            # unfortunately, this typically leads to less sparse XX, XX+Q, and therefore slower Cholesky, or even out-of-memory
            # TODO add tcrossprod(R0) where R0 is the minimal selection of columns s.t. XX + R0 R0' is positive definite
            #   i.e. restrict R0 to a selection of constraints that involve all non-observed levels
            #   and perhaps use weights to further reduce the condition number
            if (!add.eps.I && is.null(R)) {
              if (mc[["type"]] == "block") {
                noninf <- b_apply(mc[["mcs"]], \(x) isFALSE(x[["informative.prior"]]))
                if (any(noninf))
                  message("Error in setting up a multivariate normal sampler, possibly due to ",
                       "a singular precision matrix. You may try to replace the non-informative prior ",
                       "of model component(s) '", paste0(names(mc[["mcs"]])[noninf], collapse="', '"),
                       "' by a (weakly) informative prior. Alternatively set 'add.eps.I' to TRUE in create_sampler control argument.")
              }
              stop(e)
            } else {
              if (!add.eps.I && class(XX_Q)[1L] == "dsCMatrix" && XX_Q@Dim[1L] >= 500L) {
                message("adding outer product of constraint matrix to the precision matrix ",
                  "to improve its condition number\n",
                  "note that this may be inefficient if the resulting matrix is much less sparse"
                )
              }
            }
          }
        )
      }
      if (is.null(mc[["MVNsampler"]])) {  # 2nd attempt or add.outer.R
        if (is.null(R)) stop("singular precision matrix and no constraints")
        if (update.XX) {
          if (add.eps.I)
            mat_sum <- make_mat_sum(M0=control[["eps"]] * CdiagU(ncol(XX)), M1=XX, M2=Q)
          else
            mat_sum <- make_mat_sum(M0=tcrossprod(R0), M1=XX, M2=Q)
          XX_Q <- mat_sum(XX, Q)
        } else {
          if (add.eps.I)
            mat_sum <- make_mat_sum(M0=economizeMatrix(XX + control[["eps"]] * CdiagU(ncol(XX)), symmetric=TRUE), M1=Q)
          else
            mat_sum <- make_mat_sum(M0=economizeMatrix(XX + tcrossprod(R0), symmetric=TRUE), M1=Q)
          XX_Q <- mat_sum(Q)
        }
        tryCatch(
          suppressWarnings(
            mc$MVNsampler <- create_TMVN_sampler(
              Q=XX_Q, update.Q=TRUE, name=mc[["name"]],
              constraints=constraints,
              chol.control=mc[["e"]]$control[["chol.control"]]
            )
          ),
          error = function(e) {
            if (mc[["type"]] == "block") {
              noninf <- b_apply(mc[["mcs"]], \(x)
                isFALSE(x[["informative.prior"]]) || (x[["type"]] == "gen" && x[["gl"]] && isFALSE(x$glp$informative.prior))
              )
              if (any(noninf))
                warn("Error in setting up a multivariate normal sampler, possibly due to ",
                     "a singular precision matrix. You may try to replace the non-informative prior ",
                     "of model component(s) '", paste0(names()[noninf], collapse="', '"),
                     "' by a (weakly) informative prior.")
            }
            stop(e)
          }
        )
      }
    }

    keep.ms <- TRUE
    if (mc[["type"]] == "gen" && mc[["unit_Q"]]) {  # unit QA, scalar Qv (with unit Q0), and no constraints
      if (is_unit_diag(XX)) {
        # TODO cleaner and more efficient handling of unit XX case (unit_Q + unit XX --> scalar update)
        XX.expanded <- if (is_unit_ddi(XX)) expand_unit_ddi(XX) else XX
        mc$update <- function(XX, QA, Qv, w) list(Q=XX.expanded, Imult=Qv*w)
        keep.kp <- keep.ms <- FALSE
      } else {
        if (class(XX)[1L] == class(if (update.XX) mat_sum(XX, Q) else mat_sum(Q))[1L]) {
          keep.kp <- keep.ms <- FALSE
          mc$update <- function(XX, QA, Qv, w) list(Q=XX, Imult=Qv*w)
        } else {  # keep mat_sum, kron_prod because for some reason XX has type incompatible with ch
          if (update.XX)
            mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(XX, kron_prod(QA, Qv * w)), Imult=0)
          else
            mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(kron_prod(QA, Qv * w)), Imult=0)
        }
      }
    } else {
      if (mc[["type"]] == "block" || mc[["q0"]] == 1L) {  # scalar Qv, or block
        if (mc[["type"]] != "block" && is.null(mc[["AR1.inferred"]])) keep.kp <- FALSE
        if (update.XX)
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(XX, QA, w2=Qv*w), Imult=0)
        else
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(QA, w1=Qv*w), Imult=0)
      } else {  # q0 > 1, no block
        if (update.XX)
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(XX, kron_prod(QA, Qv), w2=w), Imult=0)
        else
          mc$update <- function(XX, QA, Qv, w) list(Q=mat_sum(kron_prod(QA, Qv), w1=w), Imult=0)
      }
    }
    if (keep.ms) mc$mat_sum <- mat_sum
  }

  if (mc[["type"]] != "block" && mc[["in_block"]]) keep.kp <- TRUE
  if (keep.kp) mc$kron_prod <- kron_prod

}
