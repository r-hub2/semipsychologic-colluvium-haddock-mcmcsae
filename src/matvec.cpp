#ifndef PTR_RESTRICT
#  if defined(__clang__) || defined(__GNUC__)
#    define PTR_RESTRICT __restrict__
#  elif defined(_MSC_VER)
#    define PTR_RESTRICT __restrict
#  else
#    define PTR_RESTRICT
#  endif
#endif

#include <RcppEigen.h>
#include <vector>
#include <utility>
using namespace Rcpp;
using namespace Eigen;


// generic update_op supporting Eigen Block references and Maps
template <typename DestExpr, typename SrcExpr>
inline void update_op(DestExpr&& y_dest, const SrcExpr& expr, const int op) {
  switch(op) {
  case 0: y_dest.noalias() = expr; break;
  case 1: y_dest.noalias() += expr; break;
  case 2: y_dest.noalias() -= expr; break;
  default: stop("invalid op code; use 0 (set), 1 (add), or 2 (subtract)");
  }
}

// extract dimensions across supported types
inline void get_matrix_dims(SEXP A_in, int& rows, int& cols) {
  if (Rf_isS4(A_in)) {
    static SEXP sym_Dim = Rf_install("Dim");
    const int* dim = INTEGER(R_do_slot(A_in, sym_Dim));
    rows = dim[0]; cols = dim[1];
  } else if (Rf_isMatrix(A_in)) {
    const int* dim = INTEGER(Rf_getAttrib(A_in, R_DimSymbol));
    rows = dim[0]; cols = dim[1];
  } else if (TYPEOF(A_in) == REALSXP) {
    rows = Rf_length(A_in);
    cols = rows;
  } else {
    stop("unsupported matrix type in list");
  }
}

inline void get_rhs_dims(SEXP x_in, int& x_rows, int& K, bool& is_matrix) {
  is_matrix = Rf_isMatrix(x_in);
  if (is_matrix) {
    const int* dim = INTEGER(Rf_getAttrib(x_in, R_DimSymbol));
    x_rows = dim[0];
    K = dim[1];
  } else if (Rf_isReal(x_in)) {
    x_rows = Rf_length(x_in);
    K = 1;
  } else {
    stop("x must be a numeric vector or matrix");
  }
}

inline std::vector<SEXP> get_matrix_list(SEXP A) {
  if (TYPEOF(A) == VECSXP) {
    List L(A);
    return std::vector<SEXP>(L.begin(), L.end());
  }
  return std::vector<SEXP>{A};
}


void dense_mvprod_op(SEXP A_in, const int rows, const int cols, Map<const MatrixXd>& X, Map<MatrixXd>& Y,
                     const int op, const int x_offset, const int K) {
  Map<const MatrixXd> A(REAL(A_in), rows, cols);
  update_op(Y.block(0, 0, rows, K), A * X.block(x_offset, 0, cols, K), op);
}

void sparse_mvprod_op(SEXP A_in, const int rows, const int cols, Map<const MatrixXd>& X, Map<MatrixXd>& Y,
                      const int op, const int x_offset, const int K) {
  static SEXP sym_p = Rf_install("p");
  static SEXP sym_i = Rf_install("i");
  static SEXP sym_x = Rf_install("x");

  SEXP p_sexp = R_do_slot(A_in, sym_p);
  SEXP i_sexp = R_do_slot(A_in, sym_i);
  SEXP x_sexp = R_do_slot(A_in, sym_x);

  const int nnz = Rf_length(i_sexp);
  const Map<SparseMatrix<double>> A(rows, cols, nnz, INTEGER(p_sexp), INTEGER(i_sexp), REAL(x_sexp));

  auto x_block = X.block(x_offset, 0, cols, K);
  auto y_block = Y.block(0, 0, rows, K);

  if (Rf_inherits(A_in, "dsCMatrix")) {
    update_op(y_block, A.selfadjointView<Upper>() * x_block, op);
  } else if (Rf_inherits(A_in, "dtCMatrix")) {
    update_op(y_block, A.triangularView<Upper>() * x_block, op);
  } else {
    update_op(y_block, A * x_block, op);
  }
}

void diag_mvprod_op(SEXP A_in, const int n, Map<const MatrixXd>& X, Map<MatrixXd>& Y,
                     const int op, const int x_offset, const int y_offset, const int K) {
  static SEXP sym_x    = Rf_install("x");
  static SEXP sym_diag = Rf_install("diag");

  bool is_unit = false;
  double* d_ptr = nullptr;

  if (Rf_isS4(A_in) && Rf_inherits(A_in, "ddiMatrix")) {
    SEXP diag_slot = R_do_slot(A_in, sym_diag);
    const char* diag_str = CHAR(STRING_ELT(diag_slot, 0));
  
    if (diag_str[0] == 'U') {
      is_unit = true;
    } else {
      d_ptr = REAL(R_do_slot(A_in, sym_x));
    }
  } else if (TYPEOF(A_in) == REALSXP) {
    d_ptr = REAL(A_in);
  } else {
    stop("input must be a 'ddiMatrix' or a numeric vector");
  }

  auto x_block = X.block(x_offset, 0, n, K);
  auto y_block = Y.block(y_offset, 0, n, K);

  if (is_unit) {
    update_op(y_block, x_block, op);
  } else {
    Map<const VectorXd> d(d_ptr, n);
    update_op(y_block, d.asDiagonal() * x_block, op);
  }
}

void tab_mvprod_op(SEXP A_in, const int rows, SEXP x_in, SEXP y,
                   const int op = 0, const int x_offset = 0,
                   const int K = 1, const int X_rows = 0, const int Y_rows = 0) {

  static SEXP sym_perm    = Rf_install("perm");
  static SEXP sym_reduced = Rf_install("reduced");
  static SEXP sym_num     = Rf_install("num");
  static SEXP sym_x       = Rf_install("x");

  const bool reduced  = LOGICAL(R_do_slot(A_in, sym_reduced))[0];
  const bool num      = LOGICAL(R_do_slot(A_in, sym_num))[0];
  const int* perm_ptr = INTEGER(R_do_slot(A_in, sym_perm));

  const double* PTR_RESTRICT x_slot_ptr = num ? REAL(R_do_slot(A_in, sym_x)) : nullptr;
  const double* PTR_RESTRICT x_ptr = REAL(x_in) + x_offset;
  double* PTR_RESTRICT y_ptr = REAL(y);

  for (int k = 0; k < K; ++k) {
    const double* PTR_RESTRICT x_col = x_ptr + k * X_rows;
    double* PTR_RESTRICT y_col = y_ptr + k * Y_rows;

  #define RUN_LOOP(VAL_EXPR)                                      \
    if (op == 0) {                                                \
      for (int i = 0; i < rows; ++i) y_col[i] = (VAL_EXPR);       \
    } else if (op == 1) {                                         \
      for (int i = 0; i < rows; ++i) y_col[i] += (VAL_EXPR);      \
    } else if (op == 2) {                                         \
      for (int i = 0; i < rows; ++i) y_col[i] -= (VAL_EXPR);      \
    } else {                                                      \
      stop("invalid op code");                                    \
    }

    if (reduced) {
      RUN_LOOP((perm_ptr[i] >= 0) ? x_col[perm_ptr[i]] : 0.0)
    } else if (num) {
      RUN_LOOP(x_slot_ptr[i] * x_col[perm_ptr[i]])
    } else {
      RUN_LOOP(x_col[perm_ptr[i]])
    }
  #undef RUN_LOOP
  }
}

//' Compute the product of a (sparse) matrix with a vector or dense matrix
//'  and possibly add to or subtract (in place) from an existing vector or
//'  dense matrix
//'
//' Note that it is assumed that \code{x} and \code{y} do not alias.
//'
//' @noRd
//' @param A a matrix object. Supported types are matrix,
//'  sparse (dgC/dsC/dtC) matrices, tabMatrix and diagonal
//'  matrices (ddi or numeric vector interpreted as such),
//'  and a list of any combination of such matrices with equal
//'  numbers of rows, interpreted as the column-binded matrix.
//'  For symmetric and triangular sparse matrices currently only
//'  upper dsC/dtC, non-unit-diagonal dtC matrices are supported.
//' @param x a numeric vector or dense matrix with length/number of rows equal
//'  to the (total) number of columns of A.
//' @param y an optional numeric vector or dense matrix with length/number of
//'  rows equal to the number of rows of A, which will be set or updated
//'  in place.
//' @param op an integer specifying the operation to apply on y.
//'  If 0 y is set to Ax, 1 adds Ax to y and 2 subtracts Ax from y.
//'  Defaults to 0, and can only be 0 if no y is passed.
//' @returns A numeric vector or dense matrix containing the result
//'  \code{Ax (+/- y)}. The vector is newly allocated only if no y is passed,
//'  otherwise y will be updated in place.
//' @noRd
// [[Rcpp::export(rng=false)]]
SEXP m_op(SEXP A, SEXP x, SEXP y = R_NilValue, const int op = 0) {

  int x_rows, K;
  bool is_x_matrix;
  get_rhs_dims(x, x_rows, K, is_x_matrix);

  std::vector<SEXP> A_list = get_matrix_list(A);
  const int n_matrices = A_list.size();
  if (n_matrices == 0) stop("matrix input cannot be empty");

  std::vector<std::pair<int,int>> dims(n_matrices);
  int rows = 0, total_cols = 0;
  for (int k = 0; k < n_matrices; ++k) {
    int r, c;
    get_matrix_dims(A_list[k], r, c);
    dims[k] = {r, c};
    if (k == 0) rows = r; else if (r != rows) stop("row dimension mismatch across matrices");
    total_cols += c;
  }

  if (x_rows != total_cols) {
    stop("dimension mismatch: row size of x must equal total columns of A");
  }

  Rcpp::RObject y_out;

  if (Rf_isNull(y)) {
    if (op != 0) stop("cannot add or subtract when y is NULL");
    if (is_x_matrix) {
      y_out = Rcpp::NumericMatrix(rows, K);
    } else {
      y_out = Rcpp::NumericVector(rows);
    }
  } else {
    int y_rows, y_cols;
    bool y_is_matrix;
    get_rhs_dims(y, y_rows, y_cols, y_is_matrix);
    if (y_rows != rows || y_cols != K) stop("dimension mismatch between y and expected output shape");
    y_out = y;
  }

  Map<const MatrixXd> X_map(REAL(x), x_rows, K);
  Map<MatrixXd> Y_map(REAL(SEXP(y_out)), rows, K);

  int x_offset = 0;
  for (int k = 0; k < n_matrices; ++k) {
    SEXP curr_A = A_list[k];
    const int r = dims[k].first;
    const int c = dims[k].second;

    int curr_op = (op == 0) ? ((k == 0) ? 0 : 1) : op;

    if (Rf_isS4(curr_A)) {
      if (Rf_inherits(curr_A, "tabMatrix")) {
        tab_mvprod_op(curr_A, r, x, y_out, curr_op, x_offset, K, x_rows, rows);
      } else if (Rf_inherits(curr_A, "ddiMatrix")) {
        diag_mvprod_op(curr_A, c, X_map, Y_map, curr_op, x_offset, 0, K);
      } else {
        sparse_mvprod_op(curr_A, r, c, X_map, Y_map, curr_op, x_offset, K);
      }
    } else if (Rf_isMatrix(curr_A)) {
      dense_mvprod_op(curr_A, r, c, X_map, Y_map, curr_op, x_offset, K);
    } else if (TYPEOF(curr_A) == REALSXP) {
      diag_mvprod_op(curr_A, c, X_map, Y_map, curr_op, x_offset, 0, K);
    }

    x_offset += c;
  }

  return y_out;
}


void dense_crossprod_op(SEXP A_in, const int rows, const int cols, Map<const MatrixXd>& X, Map<MatrixXd>& Y,
                        const int op, const int y_offset, const int K) {
  Map<const MatrixXd> A(REAL(A_in), rows, cols);
  update_op(Y.block(y_offset, 0, cols, K), A.transpose() * X.block(0, 0, rows, K), op);
}

void sparse_crossprod_op(SEXP A_in, const int rows, const int cols, Map<const MatrixXd>& X, Map<MatrixXd>& Y,
                         const int op, const int y_offset, const int K) {
  static SEXP sym_p = Rf_install("p");
  static SEXP sym_i = Rf_install("i");
  static SEXP sym_x = Rf_install("x");

  SEXP p_sexp = R_do_slot(A_in, sym_p);
  SEXP i_sexp = R_do_slot(A_in, sym_i);
  SEXP x_sexp = R_do_slot(A_in, sym_x);

  const int nnz = Rf_length(i_sexp);
  const Map<SparseMatrix<double>> A(rows, cols, nnz, INTEGER(p_sexp), INTEGER(i_sexp), REAL(x_sexp));

  auto x_block = X.block(0, 0, rows, K);
  auto y_block = Y.block(y_offset, 0, cols, K);

  if (Rf_inherits(A_in, "dsCMatrix")) {
    update_op(y_block, A.selfadjointView<Upper>() * x_block, op);
  } else if (Rf_inherits(A_in, "dtCMatrix")) {
    update_op(y_block, A.triangularView<Upper>().transpose() * x_block, op);
  } else {
    update_op(y_block, A.transpose() * x_block, op);
  }
}

void diag_crossprod_op(SEXP A_in, const int n, Map<const MatrixXd>& X, Map<MatrixXd>& Y,
                       const int op, const int y_offset, const int K) {
  diag_mvprod_op(A_in, n, X, Y, op, 0, y_offset, K);
}

void tab_crossprod_op(SEXP A_in, const int rows, const int cols, SEXP x_in, SEXP y,
                      const int op = 0, const int y_offset = 0,
                      const int K = 1, const int X_rows = 0, const int Y_rows = 0) {

  static SEXP sym_perm    = Rf_install("perm");
  static SEXP sym_reduced = Rf_install("reduced");
  static SEXP sym_num     = Rf_install("num");
  static SEXP sym_x       = Rf_install("x");

  const bool reduced  = LOGICAL(R_do_slot(A_in, sym_reduced))[0];
  const bool num      = LOGICAL(R_do_slot(A_in, sym_num))[0];
  const int* perm_ptr = INTEGER(R_do_slot(A_in, sym_perm));

  const double* PTR_RESTRICT x_slot_ptr = num ? REAL(R_do_slot(A_in, sym_x)) : nullptr;
  const double* PTR_RESTRICT x_ptr = REAL(x_in);
  double* PTR_RESTRICT y_ptr = REAL(y) + y_offset;

  for (int k = 0; k < K; ++k) {
    const double* PTR_RESTRICT x_col = x_ptr + k * X_rows;
    double* PTR_RESTRICT y_col = y_ptr + k * Y_rows;

    if (op == 0) {
      std::fill_n(y_col, cols, 0.0);
    }

    if (reduced) {
      if (op == 0 || op == 1) {
        for (int i = 0; i < rows; ++i) {
          int p = perm_ptr[i];
          if (p >= 0) y_col[p] += x_col[i];
        }
      } else if (op == 2) {
        for (int i = 0; i < rows; ++i) {
          int p = perm_ptr[i];
          if (p >= 0) y_col[p] -= x_col[i];
        }
      }
    } else if (num) {
      if (op == 0 || op == 1) {
        for (int i = 0; i < rows; ++i) {
          y_col[perm_ptr[i]] += x_slot_ptr[i] * x_col[i];
        }
      } else if (op == 2) {
        for (int i = 0; i < rows; ++i) {
          y_col[perm_ptr[i]] -= x_slot_ptr[i] * x_col[i];
        }
      }
    } else {
      if (op == 0 || op == 1) {
        for (int i = 0; i < rows; ++i) {
          y_col[perm_ptr[i]] += x_col[i];
        }
      } else if (op == 2) {
        for (int i = 0; i < rows; ++i) {
          y_col[perm_ptr[i]] -= x_col[i];
        }
      }
    }
  }
}

//' Compute transposed-matrix-vector/matrix product and possibly add to or
//' subtract (in place) from an existing vector/matrix
//' Compute the transposed product of a (sparse) matrix with a vector or dense
//' matrix and possibly add to or subtract (in place) from an existing vector or
//' dense matrix
//'
//' Note that it is assumed that \code{x} and \code{y} do not alias.
//'
//' @noRd
//' @param A a matrix object. Supported types are matrix,
//'  sparse (dgC/dsC/dtC) matrices, tabMatrix and diagonal
//'  matrices (ddi or numeric vector interpreted as such),
//'  and a list of any combination of such matrices with equal
//'  numbers of rows, interpreted as the column-binded matrix.
//'  For symmetric and triangular sparse matrices currently only
//'  upper dsC/dtC, non-unit-diagonal dtC matrices are supported.
//' @param x a numeric vector of length the number of rows of A, or
//'  a dense matrix with the same number of rows.
//' @param y an optional numeric vector or dense matrix with length/number of
//'  columns equal to the number of rows of A, which will be set or updated
//'  in place.
//' @param op an integer specifying the operation to apply on y.
//'  If 0 y is set to Ax, 1 adds Ax to y and 2 subtracts Ax from y.
//'  Defaults to 0, and can only be 0 if no y is passed.
//' @returns A numeric vector or dense matrix containing the result
//'  \code{t(A)x (+/- y)}. The vector is newly allocated only if no y is passed,
//'  otherwise y will be updated in place.
//' @noRd
// [[Rcpp::export(rng=false)]]
SEXP tm_op(SEXP A, SEXP x, SEXP y = R_NilValue, const int op = 0) {

  int x_rows, K;
  bool is_x_matrix;
  get_rhs_dims(x, x_rows, K, is_x_matrix);

  std::vector<SEXP> A_list = get_matrix_list(A);
  const int n_matrices = A_list.size();
  if (n_matrices == 0) stop("matrix input cannot be empty");

  std::vector<std::pair<int,int>> dims(n_matrices);
  int total_rows = 0, total_cols = 0;
  for (int k = 0; k < n_matrices; ++k) {
    int r, c;
    get_matrix_dims(A_list[k], r, c);
    dims[k] = {r, c};
    if (k == 0) total_rows = r;
    else if (r != total_rows) stop("row dimension mismatch across matrices");
    total_cols += c;
  }

  if (x_rows != total_rows) {
    stop("dimension mismatch: row size of x must equal total rows of A");
  }

  Rcpp::RObject y_out;

  if (Rf_isNull(y)) {
    if (op != 0) stop("cannot add or subtract when y is NULL");
    if (is_x_matrix) {
      y_out = Rcpp::NumericMatrix(total_cols, K);
    } else {
      y_out = Rcpp::NumericVector(total_cols);
    }
  } else {
    int y_rows, y_cols;
    bool y_is_matrix;
    get_rhs_dims(y, y_rows, y_cols, y_is_matrix);
    if (y_rows != total_cols || y_cols != K) stop("dimension mismatch between y and expected output shape");
    y_out = y;
  }

  Map<const MatrixXd> X_map(REAL(x), x_rows, K);
  Map<MatrixXd> Y_map(REAL(SEXP(y_out)), total_cols, K);

  int y_offset = 0;
  for (int k = 0; k < n_matrices; ++k) {
    SEXP curr_A = A_list[k];
    const int r = dims[k].first;
    const int c = dims[k].second;

    if (Rf_isS4(curr_A)) {
      if (Rf_inherits(curr_A, "tabMatrix")) {
        tab_crossprod_op(curr_A, r, c, x, y_out, op, y_offset, K, x_rows, total_cols);
      } else if (Rf_inherits(curr_A, "ddiMatrix")) {
        diag_crossprod_op(curr_A, c, X_map, Y_map, op, y_offset, K);
      } else {
        sparse_crossprod_op(curr_A, r, c, X_map, Y_map, op, y_offset, K);
      }
    } else if (Rf_isMatrix(curr_A)) {
      dense_crossprod_op(curr_A, r, c, X_map, Y_map, op, y_offset, K);
    } else if (TYPEOF(curr_A) == REALSXP) {
      diag_crossprod_op(curr_A, c, X_map, Y_map, op, y_offset, K);
    }

    y_offset += c;
  }

  return y_out;
}


//' Update a vector in place
//'
//' @noRd
//' @param y a numeric vector to be updated in place.
//' @param op an integer specifying the operation to apply: 0 sets y to x,
//'  1 adds x to y, 2 subtracts x from y, 3 multiplies y by x, and
//'  4 divides y by x.
//' @param x a numeric vector of the same length as y, or a numeric scalar.
//' @returns The updated vector y.
// [[Rcpp::export(rng = false)]]
SEXP v_op(SEXP y, const int op, SEXP x) {

  if (!Rf_isReal(y) || !Rf_isReal(x)) {
    stop("both inputs must be numeric vectors");
  }

  const int n = Rf_length(y);
  const int n_x = Rf_length(x);

  Map<ArrayXd> vecY(REAL(y), n);
  
  if (n_x == 1) {  // scalar broadcasting
    const double val = REAL(x)[0];
    if      (op == 0) vecY = val;
    else if (op == 1) vecY += val;
    else if (op == 2) vecY -= val;
    else if (op == 3) vecY *= val;
    else if (op == 4) vecY /= val;
    else Rcpp::stop("invalid op code");
  } else if (n_x == n) {
    Map<const ArrayXd> vecX(REAL(x), n);
    if      (op == 0) vecY = vecX;         // SET
    else if (op == 1) vecY += vecX;        // ADD
    else if (op == 2) vecY -= vecX;        // SUB
    else if (op == 3) vecY *= vecX;        // MUL
    else if (op == 4) vecY /= vecX;        // DIV
    else stop("invalid op code; use 0 (SET), 1 (ADD), 2 (SUB), 3 (MUL) or 4 (DIV)");
  } else {
    stop("length of x (%d) must be 1 or equal to length of y (%d)", n_x, n);
  }

  return y;
}
