
#include <RcppEigen.h>
using namespace Rcpp;


//’ Create a tabMatrix
//’
//’ @param Dim integer vector of size two containing the number of rows and columns.
//’ @param reduced whether zero columns are allowed.
//’ @param perm integer vector of column indices corresponding to the unique (if any)
//’  non-zero entry per row.
//’ @param num whether the non-zeros can be different from 1.
//’ @param x in case num is true, the numeric vector of nonzero values.
//’ @returns A matrix of class \code{tabMatrix}.
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

//’ Matrix product of a sparse tabMatrix object with a vector
//’
//’ @param A a tabMatrix.
//’ @param y a numeric vector.
//’ @param ignore_x whether to use only the indicator part of the tabMatrix (for expansion).
//’ @returns The vector \code{Ay}.
// [[Rcpp::export(rng=false)]]
Rcpp::NumericVector Ctab_numeric_prod(const SEXP A, const NumericVector & y, const bool ignore_x = false) {
  if (!Rf_isS4(A) || !Rf_inherits(A, "tabMatrix")) stop("A is not a tabMatrix");
  const S4 A_s4(A);  // Only do this once
  const IntegerVector perm = A_s4.slot("perm");
  const IntegerVector Dim = A_s4.slot("Dim");
  const int n = perm.size();
  if (Dim[1] != y.size()) stop("incompatible dimensions");
  NumericVector out(no_init(n));
  const bool reduced = as<bool>(A_s4.slot("reduced"));
  const bool num = as<bool>(A_s4.slot("num"));
  const int* perm_ptr = perm.begin();
  const double* y_ptr = y.begin();
  double* out_ptr = out.begin();
  if (reduced) {
    for (int i = 0; i < n; ++i) {
      out_ptr[i] = (perm_ptr[i] >= 0) ? y_ptr[perm_ptr[i]] : 0.0;
    }
  } else if (num && !ignore_x) {
    const Rcpp::NumericVector x = A_s4.slot("x");
    const double* x_ptr = x.begin();
    for (int i = 0; i < n; ++i) {
      out_ptr[i] = x_ptr[i] * y_ptr[perm_ptr[i]];
    }
  } else {
    for (int i = 0; i < n; ++i) {
      out_ptr[i] = y_ptr[perm_ptr[i]];
    }
  }
  return out;
}

//’ Matrix product of a sparse tabMatrix with a matrix
//’
//’ @param A a tabMatrix.
//’ @param y a matrix.
//’ @returns The matrix product \code{Ay}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Ctab_dense_prod(const SEXP A, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (!Rf_isS4(A) || !Rf_inherits(A, "tabMatrix")) stop("A is not a tabMatrix");
  const S4 A_S4(A);
  const IntegerVector perm(A_S4.slot("perm"));
  const IntegerVector Dim(A_S4.slot("Dim"));
  if (Dim[1] != y.rows()) stop("incompatible dimensions");
  const int n = perm.size();
  const bool reduced(::Rf_asLogical(A_S4.slot("reduced")));
  const bool num(::Rf_asLogical(A_S4.slot("num")));
  Eigen::MatrixXd out(n, y.cols());
  if (reduced) {
    for (int i = 0; i < n; i++) {
      if (perm[i] < 0) {
        out.row(i).setZero();
      } else {
        out.row(i) = y.row(perm[i]);
      }
    }
  } else if (num) {
    const NumericVector x(A_S4.slot("x"));
    for (int i = 0; i < n; i++) {
      out.row(i) = x[i] * y.row(perm[i]);
    }
  } else {
    for (int i = 0; i < n; i++) {
      out.row(i) = y.row(perm[i]);
    }
  }
  return out;
}

//’ Matrix product of a matrix and the transpose of a tabMatrix
//’
//’ @param y a matrix.
//’ @param A a tabMatrix.
//’ @returns The matrix product \code{yA'}.
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

//’ Crossproduct of a tabMatrix with a vector
//’
//’ @param A a tabMatrix.
//’ @param y a numeric vector.
//’ @returns The vector \code{A'y}.
// [[Rcpp::export(rng=false)]]
NumericVector Ctab_numeric_crossprod(const SEXP A, const NumericVector & y) {
  if (!Rf_isS4(A) || !Rf_inherits(A, "tabMatrix")) stop("A is not a tabMatrix");
  const S4 A_S4(A);
  const IntegerVector perm(A_S4.slot("perm"));
  const IntegerVector Dim(A_S4.slot("Dim"));
  const int n = y.size();
  if (Dim[0] != n) stop("incompatible dimensions");
  NumericVector out(Dim[1]);
  const bool reduced(::Rf_asLogical(A_S4.slot("reduced")));
  const bool num(::Rf_asLogical(A_S4.slot("num")));
  if (reduced) {
    for (int i = 0; i < n; i++) {
      if (perm[i] >= 0) {
        out[perm[i]] += y[i];
      }
    }
  } else if (num) {
    const NumericVector x(A_S4.slot("x"));
    for (int i = 0; i < n; i++) {
      out[perm[i]] += x[i]*y[i];
    }
  } else {
    for (int i = 0; i < n; i++) {
      out[perm[i]] += y[i];
    }
  }
  return out;
}

//’ Crossproduct of a tabMatrix with a matrix
//’
//’ @param A a tabMatrix.
//’ @param y a numeric matrix.
//’ @returns The vector \code{A'y}.
// [[Rcpp::export(rng=false)]]
NumericMatrix Ctab_dense_crossprod(const SEXP A, const NumericMatrix & y) {
  if (!Rf_isS4(A) || !Rf_inherits(A, "tabMatrix")) stop("A is not a tabMatrix");
  const S4 A_S4(A);
  const IntegerVector perm(A_S4.slot("perm"));
  const IntegerVector Dim(A_S4.slot("Dim"));
  const int yrows = y.rows();
  if (Dim[0] != yrows) stop("incompatible dimensions");
  const int ycols = y.cols();
  NumericMatrix out(Dim[1], ycols);
  const bool reduced(::Rf_asLogical(A_S4.slot("reduced")));
  const bool num(::Rf_asLogical(A_S4.slot("num")));
  if (reduced) {
    for (int i = 0; i < yrows; i++) {
      if (perm[i] >= 0) {
        for (int j = 0; j < ycols; j++) {
          out(perm[i], j) += y(i, j);
        }
      }
    }
  } else if (num) {
    const NumericVector x(A_S4.slot("x"));
    for (int j = 0; j < ycols; j++) {
      for (int i = 0; i < yrows; i++) {
        out(perm[i], j) += x[i] * y(i, j);
      }
    }
  } else {
    for (int j = 0; j < ycols; j++) {
      for (int i = 0; i < yrows; i++) {
        out(perm[i], j) += y(i, j);
      }
    }
  }
  return out;
}

//’ Values of unary crossprod of tabMatrix
//’
//’ @param A a tabMatrix.
//’ @returns The entries of the diagonal matrix representing the unary cross-product of A.
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

//’ Coerce a \code{tabMatrix} to a \code{dgCMatrix}
//’
//’ @param M a \code{tabMatrix} object.
//’ @returns The same matrix as a \code{dgCMatrix}.
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

//’ Coerce a \code{tabMatrix} to a \code{matrix}
//’
//’ @param M a \code{tabMatrix} object.
//’ @returns The same matrix as an ordinary dense matrix.
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
