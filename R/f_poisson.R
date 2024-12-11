
#' @export
#' @rdname mcmcsae-family
f_poisson <- function(link="log", size=100) {
  family <- "poisson"
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv
  size <- as.numeric(size)
  if (length(size) != 1L || is.na(size) || size <= 0) stop("'size' must be a positive scalar")
  # NB size is only used as an internal negbinomial parameter, and only for model fitting
  log.size <- log(size)
  make_llh <- function(y) {
    llh_0 <- -sum(lgamma(y + 1))
    # Poisson log-likelihood, first remove internal offset from linear predictor p[["e_"]]
    function(p) {
      eta <- p[["e_"]] + log.size
      llh_0 + sum(y * eta - exp(eta))
    }
  }
  make_llh_i <- function(y) {
    llh_0_i <- -lgamma(y + 1)
    # NB e_i must be the linear predictor excluding the internal intercept
    function(draws, i, e_i) {
      nr <- dim(e_i)[1L]
      rep_each(llh_0_i[i], nr) + rep_each(y[i], nr) * e_i - exp(e_i)
    }
  }
  make_rpredictive <- function(newdata) {
    if (is.null(newdata)) {
      # in-sample prediction/replication, linear predictor,
      # or custom X case, see prediction.R
      rpredictive <- function(p, lp) rpois(length(lp), lambda=exp(lp))
    } else {
      newn <- nrow(newdata)
      rpredictive <- function(p, lp) rpois(newn, lambda=exp(lp))
    }
  }
  environment()
}
