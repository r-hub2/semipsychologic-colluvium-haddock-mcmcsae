
// RcppEigen.h also pulls in Rcpp.h
#include <RcppEigen.h>
using namespace Rcpp;


//' Copy an object in R
//'
//' @noRd
//' @param x an R object, e.g. a vector, to be copied.
//' @returns A newly allocated copy of the object.
// [[Rcpp::export(rng=false)]]
SEXP copy_obj(SEXP x) {
  return Rf_duplicate(x);
}

// set x[i:(i + length(z) - 1)] to z in-place
// [[Rcpp::export(rng=false)]]
void set_in_place(SEXP xSEXP, const int i, const SEXP zSEXP) {
  const int xsize = Rf_length(xSEXP);
  const int zsize = Rf_length(zSEXP);
  if (i < 0 || i + zsize > xsize) stop("index out of bounds");
  Eigen::Map<Eigen::VectorXd> x(REAL(xSEXP), xsize);
  Eigen::Map<Eigen::VectorXd> z(REAL(zSEXP), zsize);
  x.segment(i, z.size()) = z;
}

//' Inverse of a symmetric positive definite dense matrix
//'
//' @noRd
//' @param M a symmetric positive-definite matrix.
//' @returns The inverse of \code{M}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd inverseSPD(const Eigen::Map<Eigen::MatrixXd> & M) {
  const int d = M.rows();
  if (M.cols() != d) stop("not a square matrix");
  return M.selfadjointView<Eigen::Upper>().ldlt().solve(Eigen::MatrixXd::Identity(d,d));
}

//' Cholesky factor of a symmetric positive definite dense matrix
//'
//' @noRd
//' @param M a symmetric positive definite matrix.
//' @returns The upper triangular Cholesky factor of \code{M}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Ccholesky(const Eigen::Map<Eigen::MatrixXd> M) {
  return M.llt().matrixU();
}

//' Dense triangular solve, in particular of Lt system, with vector right hand side
//'
//' @noRd
//' @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//' @param y a numeric vector.
//' @returns The solution of \code{Mx=y}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Cbacksolve(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & y) {
  if (M.cols() != y.size()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().solve(y);
}

//' Dense triangular solve, in particular of L system, with vector right hand side.
//'
//' @noRd
//' @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//' @param y a numeric vector.
//' @returns The solution of \code{t(M)x=y}.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Cforwardsolve(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & y) {
  if (M.cols() != y.size()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().transpose().solve(y);
}

//' Dense triangular solve, in particular of Lt system, with matrix right hand side
//'
//' @noRd
//' @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//' @param y a numeric matrix.
//' @returns The solution of \code{Mx=y}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd CbacksolveM(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (M.cols() != y.rows()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().solve(y);
}

//' Dense triangular solve, in particular of L system, with matrix right hand side
//'
//' @noRd
//' @param M an upper triangular matrix, typically the Lt factor of a LLt decomposition.
//' @param y a numeric matrix.
//' @returns The solution of \code{t(M)x=y}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd CforwardsolveM(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::MatrixXd> & y) {
  if (M.cols() != y.rows()) stop("incompatible dimensions");
  return M.triangularView<Eigen::Upper>().transpose().solve(y);
}

//' Inner product of two vectors
//'
//' @noRd
//' @param x a numeric vector.
//' @param y a numeric vector.
//' @returns The inner product \code{t(x)y}.
// [[Rcpp::export(rng=false)]]
double dotprodC(const Eigen::Map<Eigen::VectorXd> & x, const Eigen::Map<Eigen::VectorXd> & y) {
  if (x.size() != y.size()) stop("incompatible dimensions");
  return x.dot(y);
}

//' Fast summation of a vector by group
//'
//' @noRd
//' @param x a numeric vector.
//' @param group an integer vector of the same length as \code{x} defining the groups.
//' @param n the number of groups, assumed to be labelled \code{1:n}.
//' @returns A vector with totals of \code{x} by group.
// [[Rcpp::export(rng=false)]]
NumericVector fast_aggrC(const NumericVector & x, const IntegerVector & group, const int n) {
  const int xsize = x.size();
  if (xsize != group.size()) stop("incompatible dimensions");
  NumericVector out(n);
  for (int i = 0; i < xsize; i++) {
    out[group[i] - 1] += x[i];
  }
  return out;
}

//' Matrix product of a matrix with a sparse matrix
//'
//' @noRd
//' @param A a matrix.
//' @param B a numeric compressed, sparse, column-oriented matrix.
//' @returns The matrix product \code{AB}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_sparse_prod(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::SparseMatrix<double>> & B) {
  if (A.cols() != B.rows()) stop("incompatible dimensions");
  return A * B;
}

//' Matrix product of a matrix with a symmetric sparse (dsC) matrix
//'
//' @noRd
//' @param M a matrix.
//' @param Q a numeric symmetric compressed, sparse, column-oriented matrix.
//'  It is assumed that the upper triangular part is stored.
//' @returns The matrix product \code{MQ}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_sparseS_prod(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::SparseMatrix<double>> & Q) {
  if (M.cols() != Q.rows()) stop("incompatible dimensions");
  return M * Q.selfadjointView<Eigen::Upper>();
}

//' Matrix product of a matrix with a diagonal matrix represented by a vector
//'
//' @noRd
//' @param M a dense matrix.
//' @param d a numeric vector representing a diagonal matrix.
//' @returns The matrix product \code{M diag(d)}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_diag_prod(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & d) {
  if (M.cols() != d.size()) stop("incompatible dimensions");
  return M * d.asDiagonal();
}

//' Matrix product of the transpose of a matrix with a sparse matrix
//'
//' @noRd
//' @param A a matrix.
//' @param B a numeric compressed, sparse, column-oriented matrix.
//' @returns The matrix product \code{t(A)B}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_sparse_crossprod(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::SparseMatrix<double>> & B) {
  if (A.rows() != B.rows()) stop("incompatible dimensions");
  return A.transpose() * B;
}

//' Matrix product of the transpose of a matrix with a diagonal matrix represented by a vector
//'
//' @noRd
//' @param M a dense matrix.
//' @param d a numeric vector representing a diagonal matrix.
//' @returns The matrix product \code{M diag(d)}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_diag_crossprod(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & d) {
  if (M.rows() != d.size()) stop("incompatible dimensions");
  return M.transpose() * d.asDiagonal();
}

//' Product of a matrix with the transpose of a sparse matrix
//'
//' @noRd
//' @param y a (dense) matrix.
//' @param A a numeric compressed, sparse, column-oriented matrix.
//' @returns The tcrossproduct \code{y t(A)}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_sparse_tcrossprod(const Eigen::Map<Eigen::MatrixXd> & y, const Eigen::Map<Eigen::SparseMatrix<double>> & A) {
  if (y.cols() != A.cols()) stop("incompatible dimensions");
  return y * A.transpose();
}

//' Product of a diagonal and sparse matrix
//'
//' @noRd
//' @param x a numeric vector, representing a diagonal matrix.
//' @param A a numeric compressed, sparse, column-oriented matrix.
//' @returns The product \code{diag(x) A}.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Cdiag_sparse_prod(const Eigen::Map<Eigen::VectorXd> & x, const Eigen::Map<Eigen::SparseMatrix<double>> & A) {
  if (x.size() != A.rows()) stop("incompatible dimensions");
  return x.asDiagonal() * A;
}

// Used in closure for sparse matrix addition to compute x-slot.
// [[Rcpp::export(rng=false)]]
NumericVector sparse_sum_x(const int n,
                           const IntegerVector ind1, const IntegerVector ind2,
                           const NumericVector M1x, const NumericVector M2x,
                           const bool UD1, const bool UD2,
                           const double w1, const double w2) {
  NumericVector out(n);
  double* p_out = REAL(out);
  const int n1 = ind1.size();
  const int* p_ind1 = INTEGER(ind1);
  const double* p_M1x = REAL(M1x);

  if (UD1) {
    for (int i = 0; i < n1; i++) {
      p_out[p_ind1[i]] = w1;
    }
  } else {
    for (int i = 0; i < n1; i++) {
      p_out[p_ind1[i]] = w1 * p_M1x[i];
    }
  }

  const int n2 = ind2.size();
  if (n2 > 0) {
    const int* p_ind2 = INTEGER(ind2);
    const double* p_M2x = REAL(M2x);
    if (UD2) {
      for (int i = 0; i < n2; i++) {
        p_out[p_ind2[i]] += w2;
      }
    } else {
      for (int i = 0; i < n2; i++) {
        p_out[p_ind2[i]] += w2 * p_M2x[i];
      }
    }
  }

  return out;
}

//' Extract the diagonal of a dense matrix
//'
//' @noRd
//' @param A a numeric dense matrix.
//' @returns The diagonal of A as a vector.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd diagC(const Eigen::Map<Eigen::MatrixXd> & A) {
  return A.diagonal();
}

//' Add a vector to the diagonal of a dense matrix after copying it
//'
//' @noRd
//' @param A a numeric dense matrix.
//' @param d a numeric vector.
//' @returns Matrix A with d added to its diagonal.
// [[Rcpp::export(rng=false)]]
SEXP add_diagC(SEXP A, SEXP d) {
  SEXP A_ = Rf_duplicate(A);
  NumericMatrix M(A_);
  NumericVector d_(d);
  const int n = M.nrow();
  if (d_.size() == 1) {
    const double v = d_[0];
    for (int i = 0; i < n; ++i)
      M(i, i) += v;
  } else {
    for (int i = 0; i < n; ++i)
      M(i, i) += d_[i];
  }
  return A_;
}

//' Compute the symmetric crossprod \code{t(M)diag(q)M} for dense matrix \code{M}
//'
//' @noRd
//' @param M a numeric dense matrix.
//' @param q a numeric vector.
//' @returns The product \code{t(M)diag(q)M}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_crossprod_sym(const Eigen::Map<Eigen::MatrixXd> & M, const Eigen::Map<Eigen::VectorXd> & q) {
  if (M.rows() != q.size()) stop("incompatible input");
  const int n = M.cols();
  Eigen::VectorXd sq = q.array().sqrt();
  Eigen::MatrixXd sqM = sq.asDiagonal() * M;
  return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(sqM.transpose());
}

//' Compute the symmetric crossprod \code{t(P)QP} for permutation matrix \code{P} and symmetric sparse Q
//'
//' @noRd
//' @param Q a symmetric sparse matrix, represented by its upper triangle.
//' @param p an integer vector representing a permutation matrix P.
//' @returns The product \code{t(P)QP}.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_sym_twist(const Eigen::Map<Eigen::SparseMatrix<double>> Q, const Eigen::Map<Eigen::VectorXi> & p) {
  if (Q.rows() != p.size()) stop("incompatible dimensions");
  Eigen::SparseMatrix<double> out;
  out = Q.selfadjointView<Eigen::Upper>().twistedBy( p.asPermutation() );
  return out.triangularView<Eigen::Upper>();
}

//' Create a unit Diagonal Matrix
//'
//' @noRd
//' @param n an integer representing the dimension.
//' @returns The unit diagonal matrix of class \code{ddiMatrix} and dimension \code{n}.
// [[Rcpp::export(rng=false)]]
SEXP CdiagU(const int n) {
  S4 out("ddiMatrix");
  out.slot("Dim") = IntegerVector::create(n, n);
  out.slot("diag") = "U";
  return out;
}

//' Create a Diagonal Matrix
//'
//' @noRd
//' @param x a numeric vector representing the matrix diagonal.
//' @returns The diagonal matrix of class \code{ddiMatrix} with diagonal \code{x}.
// [[Rcpp::export(rng=false)]]
SEXP Cdiag(const NumericVector x) {
  S4 out("ddiMatrix");
  const int n = x.size();
  out.slot("Dim") = IntegerVector::create(n, n);
  out.slot("x") = x;  // NB points to x, no copy
  return out;
}

//' Compute scaled dense matrix
//'
//' @noRd
//' @param A a (dense) square matrix.
//' @param d a vector of size the number of rows or columns of \code{A}, or a scalar.
//' @returns The sparse matrix \code{DAD} with \code{D=diag(d)}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cscale_dense(const Eigen::Map<Eigen::MatrixXd> & A, const Eigen::Map<Eigen::VectorXd> & d) {
  if (d.size() == 1) {
    return (d[0] * d[0]) * A;
  } else {
    return d.asDiagonal() * A * d.asDiagonal();
  }
}

//' Compute scaled sparse matrix
//'
//' @noRd
//' @param A a sparse matrix.
//' @param d a vector of size the number of rows or columns of \code{A}, or a scalar.
//' @returns The sparse matrix \code{t(D)AD} with \code{D=diag(d)}.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Cscale_sparse(const Eigen::Map<Eigen::SparseMatrix<double>> & A, const Eigen::Map<Eigen::VectorXd> & d) {
  if (d.size() == 1) {
    return (d[0] * d[0]) * A;
  } else {
    return d.asDiagonal() * A * d.asDiagonal();
  }
}

//' Compute symmetric crossproduct of sparse and dense matrices
//'
//' @noRd
//' @param A a sparse matrix.
//' @param Q a dense matrix, assumed to be symmetric.
//' @returns The (dense) matrix \code{t(A)QA}.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Csparse_dense_crossprod_sym(const Eigen::Map<Eigen::SparseMatrix<double>> & A, const Eigen::Map<Eigen::MatrixXd> & Q) {
  if (A.rows() != Q.rows()) stop("incompatible dimensions");
  return (A.adjoint() * Q * A).selfadjointView<Eigen::Upper>();
}

//' Compute symmetric crossproduct of sparse matrices
//'
//' @noRd
//' @param A a sparse matrix.
//' @param Q a sparse matrix, assumed to be symmetric.
//' @returns The sparse matrix \code{t(A)QA}, in upper triangular view.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_crossprod_sym(const Eigen::Map<Eigen::SparseMatrix<double>> & A, const Eigen::Map<Eigen::SparseMatrix<double>> & Q) {
  if (A.rows() != Q.rows()) stop("incompatible dimensions");
  return (A.adjoint() * Q.selfadjointView<Eigen::Upper>() * A).triangularView<Eigen::Upper>();
}

//' Compute symmetric weighted crossproduct of sparse matrices
//'
//' @noRd
//' @param A a sparse matrix.
//' @param Q a vector, representing a diagonal matrix.
//' @returns The sparse matrix \code{t(A)diag(Q)A}, in upper triangular view.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_diag_crossprod_sym(const Eigen::Map<Eigen::SparseMatrix<double>> & A, const Eigen::Map<Eigen::VectorXd> & Q) {
  if (Q.size() != A.rows()) stop("incompatible dimensions");
  // it seems faster to first assign QA to a temporary sparse matrix
  Eigen::SparseMatrix<double> temp = Q.asDiagonal() * A;
  return (A.adjoint() * temp).triangularView<Eigen::Upper>();
}

//' Compute crossproduct of two sparse matrices, known to yield a symmetric (sparse) matrix
//'
//' @noRd
//' @param A a sparse matrix.
//' @param B a sparse matrix, assumed to be the product of a symmetric matrix with \code{A}.
//' @returns The sparse matrix \code{t(A)B}, in upper triangular view.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Csparse_crossprod_sym2(const Eigen::Map<Eigen::SparseMatrix<double>> & A, const Eigen::Map<Eigen::SparseMatrix<double>> & B) {
  if (A.rows() != B.rows()) stop("incompatible dimensions");
  return (A.adjoint() * B).triangularView<Eigen::Upper>();
}

// [[Rcpp::export(rng=false)]]
Rcpp::List prec2se_cor(const Eigen::Map<Eigen::MatrixXd> & Q) {
  Eigen::MatrixXd V = inverseSPD(Q);
  Eigen::VectorXd se = V.diagonal().array().sqrt();
  Eigen::VectorXd se_inv = se.array().inverse();
  V = se_inv.asDiagonal() * V * se_inv.asDiagonal(); // overwrite V by the correlation matrix
  const int d = se.size();
  Eigen::VectorXd cor(d*(d-1)/2);
  int ind = 0;
  for (int i = 1; i < d; i++) {
    cor.segment(ind, i) = V.col(i).head(i);
    ind += i;
  }
  return Rcpp::List::create(Rcpp::Named("se") = se,
                            Rcpp::Named("cor") = cor);
}

//' Compute Kronecker product of two dense matrices
//'
//' @noRd
//' @param M1 numeric dense matrix.
//' @param M2 numeric dense matrix.
//' @returns The Kronecker product of M1 and M2 as a dense matrix.
// [[Rcpp::export(rng=false)]]
Eigen::MatrixXd Cdense_kron(const Eigen::Map<Eigen::MatrixXd> & M1, const Eigen::Map<Eigen::MatrixXd> & M2) {
  return kroneckerProduct(M1, M2);  //.selfadjointView<Eigen::Lower>() for symmetric matrices --> slower
}

//' Generalisation of rep function useful for sparse kronecker product templates
//'
//' @noRd
//' @param v numeric vector.
//' @param n integer vector.
//' @param M2 numeric vector.
//' @returns A vector with replicated multiples of the elements of v.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd Crepgen(const Eigen::Map<Eigen::VectorXd> & v, const Eigen::Map<Eigen::VectorXi> & n, const Eigen::Map<Eigen::VectorXd> & M2) {
  const int q = M2.size();
  Eigen::VectorXd out(q * v.size());
  Eigen::VectorXd x;
  int ind = 0;
  int ind_out = 0;
  int ni;
  for (int i = 0; i < n.size(); i++) {
    ni = n.coeff(i);
    x = v.segment(ind, ni);
    for (int j = 0; j < q; j++) {
      out.segment(ind_out, ni) = M2.coeff(j) * x;
      ind_out += ni;
    }
    ind += ni;
  }
  return out;
}

//' Compute the number of non-zeros per column of a sparse template matrix
//'
//' This function computes the number of non-zeros per column of the sparse
//' template matrix for weighted sparse symmetric crossproducts, as computed
//' by function Ccreate_sparse_crossprod_sym.
//'
//' @noRd
//' @param X a sparse (dgC) matrix.
//' @param j1_ind integer vector containing row indices of non-zero elements in t(X)QX.
//' @param j2_ind integer vector containing column indices of non-zero elements in t(X)QX.
//' @returns An integer vector containing the number of non-zeros per column of a
//'  sparse template matrix for fast updating of the weighted symmetric crossproduct
//'  of X.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXi Cnnz_per_col_scps_template(
    const Eigen::Map<Eigen::SparseMatrix<double>> & X,
    const Eigen::VectorXi & j1_ind,
    const Eigen::VectorXi & j2_ind) {

  const int nnz = j1_ind.size();
  if (j2_ind.size() != nnz) stop("'j1_ind' and 'j2_ind' should have the same length");

  Eigen::VectorXi nnz_per_col(nnz);

  // loop over nonzero elements of dsCMatrix X'QX in the order of its i and x slots
  for (int j = 0; j < nnz; j++) {
    int nnz_in_col_prod = 0;
    Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator i1_(X, j1_ind[j]);
    Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator i2_(X, j2_ind[j]);
    while (i1_) {
      while (i2_ && i2_.index() < i1_.index()) ++i2_;
      if (i2_ && i2_.index() == i1_.index()) nnz_in_col_prod++;
      ++i1_;
    }
    nnz_per_col[j] = nnz_in_col_prod;
  }
  return nnz_per_col;
}

//' Compute a sparse template matrix for the sparse weighted symmetric crossproduct
//'
//' The sparse template matrix can be used for fast updating of the sparse weighted
//' symmetric crossproduct of X, when the diagonal matrix of weights W is subject
//' to change
//'
//' @noRd
//' @param X a sparse (dgC) matrix.
//' @param j1_ind integer vector containing row indices of nonzero elements in t(X)QX.
//' @param j2_ind integer vector containing column indices of nonzero elements in t(X)QX.
//' @param nnz_per_col integer vector containing the number of non-zeros per column
//'  of the sparse template matrix.
//' @returns A sparse matrix containing the nonzero products of columns of X.
// [[Rcpp::export(rng=false)]]
Eigen::SparseMatrix<double> Ccreate_sparse_crossprod_sym_template(
    const Eigen::Map<Eigen::SparseMatrix<double>> & X,
    const Eigen::VectorXi & j1_ind,
    const Eigen::VectorXi & j2_ind,
    const Eigen::Map<Eigen::VectorXi> & nnz_per_col) {

  const int n = X.rows();
  const int nnz = j1_ind.size();
  if (j2_ind.size() != nnz) stop("'j1_ind' and 'j2_ind' should have the same length");

  // reserve space for the sparse matrix using the computed nnz's per column
  Eigen::SparseMatrix<double> out(n, nnz);
  out.reserve(nnz_per_col);

  // fill the columns of the sparse matrix with X[,j1] * X[,j2]
  for (int j = 0; j < nnz; j++) {
    Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator i1_(X, j1_ind[j]);
    Eigen::Map<Eigen::SparseMatrix<double>>::InnerIterator i2_(X, j2_ind[j]);
    while (i1_) {
      while (i2_ && i2_.index() < i1_.index()) ++i2_;
      if (i2_ && i2_.index() == i1_.index()) {
        out.insert(i1_.index(), j) = i1_.value() * i2_.value();
      }
      ++i1_;
    }
  }
  return out;
}

//' Compute variances of each row of a dense matrix
//'
//' @noRd
//' @param M a (dense) matrix.
//' @returns A vector containing the row variance.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd rowVarsC(const Eigen::Map<Eigen::MatrixXd> & M) {
  const int n_rows = M.rows();
  const int n_cols = M.cols();
  Eigen::VectorXd vars(n_rows);
  for (int i = 0; i < n_rows; ++i) {
    double mean = M.row(i).mean();
    double var = 0.0;
    for (int j = 0; j < n_cols; ++j) {
      double diff = M(i, j) - mean;
      var += diff * diff;
    }
    vars[i] = var / (n_cols - 1); // sample variance (not population variance)
  }
  return vars;
}

// helper function for BART prediction: recursive traversal
int predict_tree_rec(
    const IntegerVector& var,
    const NumericVector& value,
    const NumericMatrix& X,
    NumericVector& pred,
    const std::vector<int>& idx,
    int node
) {
  // leaf
  if (var[node] == -1) {
    double v = value[node];
    for (int i : idx)
      pred[i] += v;
    return 1;
  }

  int split_var = var[node] - 1;
  double split_val = value[node];

  std::vector<int> left;
  std::vector<int> right;
  left.reserve(idx.size());
  right.reserve(idx.size());

  for (int i : idx) {
    if (X(i, split_var) <= split_val)
      left.push_back(i);
    else
      right.push_back(i);
  }

  int n_left = predict_tree_rec(
    var, value, X, pred, left, node + 1
  );

  int n_right = predict_tree_rec(
    var, value, X, pred, right, node + 1 + n_left
  );

  return 1 + n_left + n_right;
}

//' Compute predictions for a BART component
//'
//' @noRd
//' @param tree tree column of dbarts trees data frame.
//' @param var var column of dbarts trees data frame.
//' @param value value column of dbarts trees data frame.
//' @param Xnew matrix of BART covariates to predict for.
//' @param n_trees number of trees.
//' @returns A vector with raw BART predictions corresponding.
// [[Rcpp::export(rng=false)]]
NumericVector tree_predictC(
    IntegerVector tree, IntegerVector var,
    NumericVector value, NumericMatrix Xnew,
    const int n_trees) {

  const int n = Xnew.nrow();
  NumericVector pred(n);

  std::vector<int> all_idx(n);
  for (int i = 0; i < n; ++i) all_idx[i] = i;

  for (int t = 1; t <= n_trees; ++t) {

    // collect rows for this tree
    std::vector<int> rows;
    for (int i = 0; i < tree.size(); ++i)
      if (tree[i] == t) rows.push_back(i);
    if (rows.empty()) continue;

    // create compact views
    IntegerVector var_t(rows.size());
    NumericVector value_t(rows.size());

    for (size_t i = 0; i < rows.size(); ++i) {
      var_t[i]   = var[rows[i]];
      value_t[i] = value[rows[i]];
    }

    predict_tree_rec(
      var_t, value_t, Xnew, pred, all_idx, 0
    );
  }

  return pred;
}
