
#include <Rcpp.h>


extern "C" SEXP CHM_dsC_Cholesky(SEXP a, SEXP perm, SEXP super, SEXP Imult, SEXP m, SEXP LDL);
extern "C" SEXP CHMf_solve(SEXP a, SEXP b, SEXP system);
extern "C" SEXP CHMf_solve_matrix(SEXP a, SEXP b, SEXP system);
extern "C" SEXP CHMf_spsolve(SEXP a, SEXP b, SEXP system);
extern "C" SEXP CHM_update_inplace(SEXP object, SEXP parent, SEXP mult);


// [[Rcpp::export]]
SEXP cCHM_dsC_Cholesky(SEXP a, SEXP perm, SEXP super, SEXP Imult, SEXP m, SEXP LDL) {
  return CHM_dsC_Cholesky(a, perm, super, Imult, m, LDL);
}

// [[Rcpp::export]]
SEXP cCHMf_solve(SEXP a, SEXP b, SEXP system) {
  return CHMf_solve(a, b, system);
}

// [[Rcpp::export]]
SEXP cCHMf_solve_matrix(SEXP a, SEXP b, SEXP system) {
  return CHMf_solve_matrix(a, b, system);
}

// [[Rcpp::export]]
SEXP cCHMf_spsolve(SEXP a, SEXP b, SEXP system) {
  return CHMf_spsolve(a, b, system);
}

// [[Rcpp::export]]
SEXP cCHM_update_inplace(SEXP object, SEXP parent, SEXP mult) {
  return CHM_update_inplace(object, parent, mult);
}
