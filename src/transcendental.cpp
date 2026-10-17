
#include <RcppEigen.h>
using namespace Rcpp;


//' Compute sqrt(x) or sqrt(1/x) using Eigen
//'
//' @noRd
//' @param x input vector.
//' @param inv if \code{true} sqrt(1/x) is computed instead of sqrt(x).
//' @returns The computed vector sqrt(x).
// [[Rcpp::export(rng=false)]]
SEXP sqrtC(SEXP x, const bool inv = false) {
  if (!Rf_isReal(x)) Rcpp::stop("x must be a numeric (double) vector");

  const R_xlen_t n = Rf_xlength(x);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));

  // make it a matrix/array if the input was one
  Rf_setAttrib(out, R_DimSymbol, Rf_getAttrib(x, R_DimSymbol));

  Eigen::Map<Eigen::VectorXd> map_x(REAL(x), n);
  Eigen::Map<Eigen::VectorXd> map_out(REAL(out), n);

  if (inv) {
    map_out.array() = map_x.array().rsqrt();
  } else {
    map_out.array() = map_x.array().sqrt();
  }

  UNPROTECT(1);
  return out;
}

//' Compute log(x) using Eigen
//'
//' @noRd
//' @param x input vector.
//' @returns The computed vector log(x).
// [[Rcpp::export(rng=false)]]
SEXP logC(SEXP x) {
  if (!Rf_isReal(x)) Rcpp::stop("x must be a numeric (double) vector");
  
  const R_xlen_t n = Rf_xlength(x);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));
  
  // Copy dimensions (makes it a matrix/array if the input was one)
  Rf_setAttrib(out, R_DimSymbol, Rf_getAttrib(x, R_DimSymbol));
  
  Eigen::Map<Eigen::VectorXd> map_x(REAL(x), n);
  Eigen::Map<Eigen::VectorXd> map_out(REAL(out), n);

  map_out.array() = map_x.array().log();

  UNPROTECT(1);
  return out;
}

//' Compute exp(x) or exp(-x) using Eigen
//'
//' @noRd
//' @param x input vector.
//' @param neg if false (default) compute exp(x), otherwise exp(-x).
//' @returns The computed vector exp(x) or exp(-x).
// [[Rcpp::export(rng=false)]]
SEXP expC(SEXP x, const bool neg = false) {
  if (!Rf_isReal(x)) Rcpp::stop("x must be a numeric (double) vector");

  const R_xlen_t n = Rf_xlength(x);
  SEXP out = PROTECT(Rf_allocVector(REALSXP, n));

  // make it a matrix/array if the input was one
  Rf_setAttrib(out, R_DimSymbol, Rf_getAttrib(x, R_DimSymbol));

  Eigen::Map<Eigen::VectorXd> map_x(REAL(x), n);
  Eigen::Map<Eigen::VectorXd> map_out(REAL(out), n);

  if (neg) {
    map_out.array() = (-map_x.array()).exp();
  } else {
    map_out.array() = map_x.array().exp();
  }

  UNPROTECT(1);
  return out;
}

//' Compute log(1 + exp(x)) in a numerically robust way
//'
//' @noRd
//' @param x input vector.
//' @returns The computed vector log(1 + exp(x)).
//' @references
//'   M. Maechler (2012)
//'     Accurately Computing log(1 - exp(-|a|)) Assessed by the Rmpfr package.
//'     URL: \url{https://cran.r-project.org/package=Rmpfr/vignettes/log1mexp-note.pdf}.
// [[Rcpp::export(rng=false)]]
NumericVector log1pexpC(const NumericVector & x) {
  const int n = x.size();
  NumericVector out = no_init(n);
  for (int i = 0; i < n; i++) {
    if (x[i] <= -37) {
      out[i] = std::exp(x[i]);
    } else if (x[i] <= 18) {
      out[i] = std::log1p(std::exp(x[i]));
    } else if (x[i] <= 33.3) {
      out[i] = x[i] + std::exp(-x[i]);
    } else {
      out[i] = x[i];
    }
  }
  return out;
}

//' Compute column-wise log of sum of exponentials for a dense matrix
//'
//' @noRd
//' @param M numeric dense matrix.
//' @returns The vector of column-wise log-sum-exp of M.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd logSumExpColwiseC(const Eigen::MatrixXd & M) {
  Eigen::RowVectorXd maxCoeff = M.colwise().maxCoeff();
  Eigen::MatrixXd shifted = M.rowwise() - maxCoeff;
  Eigen::RowVectorXd sumExp = shifted.array().exp().colwise().sum();
  Eigen::RowVectorXd logSumExp = maxCoeff + sumExp.array().log().matrix();
  return logSumExp.transpose();
}
