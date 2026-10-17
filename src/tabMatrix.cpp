
#include <RcppEigen.h>
using namespace Rcpp;


//' Create a tabMatrix
//'
//' @noRd
//' @param Dim integer vector of size two containing the number of rows and columns.
//' @param reduced whether zero columns are allowed.
//' @param perm integer vector of column indices corresponding to the unique (if any)
//'  non-zero entry per row.
//' @param num whether the non-zeros can be different from 1.
//' @param x in case num is true, the numeric vector of non-zero values.
//' @returns A matrix of class \code{tabMatrix}.
// [[Rcpp::export(rng=false)]]
SEXP Ctab(const IntegerVector & Dim, const bool reduced, const IntegerVector & perm, const bool num, const NumericVector & x) {
  if (reduced && num) stop("'reduced' and 'num' should not both be true");
  S4 out("tabMatrix");
  out.slot("Dim") = clone(Dim);
  out.slot("reduced") = reduced;
  if (perm.size() != Dim[0]) stop("length of perm slot must equal row dimension of tabMatrix");
  out.slot("perm") = clone(perm);
  out.slot("num") = num;
  if (num && x.size() != Dim[0]) stop("length of x slot (if nonempty) must equal row dimension of tabMatrix");
  out.slot("x") = clone(x);
  return out;
}

//' Matrix product of a matrix and the transpose of a tabMatrix
//'
//' @noRd
//' @param y a matrix.
//' @param A a tabMatrix.
//' @returns The matrix product \code{y t(A)}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_tab_tcrossprod(const Eigen::Map<Eigen::MatrixXd> & y, const SEXP A) {
  if (!Rf_isS4(A) || !Rf_inherits(A, "tabMatrix")) stop("A is not a tabMatrix");
  const S4 A_S4(A);
  const IntegerVector perm(A_S4.slot("perm"));
  const IntegerVector Dim(A_S4.slot("Dim"));
  if (Dim[1] != y.cols()) stop("incompatible dimensions");
  const int n = perm.size();
  const bool reduced(::Rf_asLogical(A_S4.slot("reduced")));
  const bool num(::Rf_asLogical(A_S4.slot("num")));
  Eigen::MatrixXd out(y.rows(), n);
  if (reduced) {
    for (int i = 0; i < n; i++) {
      if (perm[i] < 0) {
        out.col(i).setZero();
      } else {
        out.col(i) = y.col(perm[i]);
      }
    }
  } else if (num) {
    const NumericVector x(A_S4.slot("x"));
    for (int i = 0; i < n; i++) {
      out.col(i) = x[i] * y.col(perm[i]);
    }
  } else {
    for (int i = 0; i < n; i++) {
      out.col(i) = y.col(perm[i]);
    }
  }
  return out;
}

//' Values of unary crossprod of tabMatrix
//'
//' @noRd
//' @param A a tabMatrix.
//' @returns The entries of the diagonal matrix representing the unary cross-product of A.
// [[Rcpp::export(rng=false)]]
NumericVector Ctab_unary_crossprod(const SEXP A) {
  if (!Rf_isS4(A) || !Rf_inherits(A, "tabMatrix")) stop("A is not a tabMatrix");
  const S4 A_S4(A);
  const IntegerVector perm(A_S4.slot("perm"));
  const IntegerVector Dim(A_S4.slot("Dim"));
  const int n = Dim[0];
  NumericVector diag(Dim[1]);
  const bool reduced(::Rf_asLogical(A_S4.slot("reduced")));
  const bool num(::Rf_asLogical(A_S4.slot("num")));
  if (reduced) {
    for (int i = 0; i < n; i++) {
      if (perm[i] >= 0) {
        diag[perm[i]]++;
      }
    }
  } else if (num) {
    const NumericVector x(A_S4.slot("x"));
    for (int i = 0; i < n; i++) {
      diag[perm[i]] += x[i]*x[i];
    }
  } else {
    for (int i = 0; i < n; i++) {
      diag[perm[i]]++;
    }
  }
  return diag;
}

//' Coerce a \code{tabMatrix} to a \code{dgCMatrix}
//'
//' @noRd
//' @param M a \code{tabMatrix} object.
//' @returns The same matrix as a \code{dgCMatrix}.
// [[Rcpp::export(rng=false)]]
SEXP Ctab2dgC(const SEXP M) {
  if (!Rf_isS4(M) || !Rf_inherits(M, "tabMatrix")) stop("M is not a tabMatrix");
  const S4 M_S4(M);
  const IntegerVector Dim(M_S4.slot("Dim"));
  const IntegerVector perm(M_S4.slot("perm"));
  const NumericVector x(M_S4.slot("x"));
  const bool reduced(::Rf_asLogical(M_S4.slot("reduced")));
  const bool num(::Rf_asLogical(M_S4.slot("num")));
  S4 out("dgCMatrix");
  out.slot("Dim") = clone(Dim);
  IntegerVector tab(Dim[1]);
  if (reduced) {
    for (int k=0; k < perm.size(); k++) {
      if (perm[k] >= 0) tab[perm[k]]++;
    }
  } else if (num) {
    for (int k=0; k < perm.size(); k++) {
      if (x[k] != 0) tab[perm[k]]++;
    }
  } else {
    for (int k=0; k < perm.size(); k++) {
      tab[perm[k]]++;
    }
  }
  int csum = 0;
  IntegerVector colpointers(Dim[1] + 1);
  for (int k=0; k < tab.size(); k++) {
    csum += tab[k];
    colpointers[k + 1] = csum;
  }
  out.slot("p") = colpointers;
  IntegerVector ind = clone(colpointers);
  IntegerVector rowpointers(csum);
  if (num) {
    NumericVector xdgC = no_init(csum);
    for (int k=0; k < perm.size(); k++) {
      if (x[k] != 0) {
        rowpointers[ind[perm[k]]] = k;
        xdgC[ind[perm[k]]] = x[k];
        ind[perm[k]]++;
      }
    }
    out.slot("x") = xdgC;
  } else {
    out.slot("x") = rep(1., csum);
    if (reduced) {
      for (int k=0; k < perm.size(); k++) {
        if (perm[k] >= 0) {
          rowpointers[ind[perm[k]]] = k;
          ind[perm[k]]++;
        }
      }
    } else {
      for (int k=0; k < perm.size(); k++) {
        rowpointers[ind[perm[k]]] = k;
        ind[perm[k]]++;
      }
    }
  }
  out.slot("i") = rowpointers;
  return out;
}

//' Coerce a \code{tabMatrix} to a \code{matrix}
//'
//' @noRd
//' @param M a \code{tabMatrix} object.
//' @returns The same matrix as an ordinary dense matrix.
// [[Rcpp::export(rng=false)]]
NumericMatrix Ctab2mat(const SEXP M) {
  if (!Rf_isS4(M) || !Rf_inherits(M, "tabMatrix")) stop("M is not a tabMatrix");
  const S4 M_S4(M);
  const IntegerVector Dim(M_S4.slot("Dim"));
  const IntegerVector perm(M_S4.slot("perm"));
  const bool reduced(::Rf_asLogical(M_S4.slot("reduced")));
  const bool num(::Rf_asLogical(M_S4.slot("num")));
  NumericMatrix out(Dim[0], Dim[1]);
  if (num) {
    const NumericVector x(M_S4.slot("x"));
    for (int i = 0; i < Dim[0]; i++) {
      out(i, perm[i]) = x[i];
    }
  } else if (reduced) {
    for (int i = 0; i < Dim[0]; i++) {
      if (perm[i] >= 0) {
        out(i, perm[i]) = 1;
      }
    }
  } else {
    for (int i = 0; i < Dim[0]; i++) {
      out(i, perm[i]) = 1;
    }
  }
  return out;
}
