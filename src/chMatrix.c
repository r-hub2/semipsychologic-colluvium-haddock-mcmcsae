
// need to start cholmod, in addition to one started by Matrix
// --> probably two cholmod workspaces; how (memory-)inefficient is this?

#include <R_ext/RS.h>
#include <Matrix/alloca.h>
#include <Matrix/cholmod.h>
#include <Matrix/stubs.c>


extern cholmod_common c;  // see mcmcsae_init.c


// make sure mcmcsae_init.c includes this:
/*
#include "Matrix/Matrix.h"

cholmod_common c;

// and at the end:
void R_init_mcmcsae(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);

  M_R_cholmod_start(&c);
}

void R_unload_mcmcsae(DllInfo *dll) {
  M_cholmod_finish(&c);
}
*/


// m method, integer
void chm_set_ordering(const int m) {
  if (m == -1) {
    // natural ordering, i.e. no permutation
    c.nmethods = 1; c.method[0].ordering = CHOLMOD_NATURAL; c.postorder = FALSE;
  } else if (m == 0) {
    c.default_nesdis = TRUE;
    c.nmethods = 0;  // the default, but without METIS since that does not seem to be available in Matrix
  } else if (m == 1) {
    // only AMD
    c.nmethods = 1; c.method[0].ordering = CHOLMOD_AMD; c.postorder = TRUE;
  } else if (m == 2) {
    // natural ordering, but with postordering
    c.nmethods = 1; c.method[0].ordering = CHOLMOD_NATURAL; c.postorder = TRUE;
  } else if (m == 3) {
    // most extensive search
    c.nmethods = 9;
  }
}

// Cholesky of dsCMatrix
// added argument m: ordering method (integer)
SEXP CHM_dsC_Cholesky(SEXP a, SEXP perm, SEXP super, SEXP Imult, SEXP m, SEXP LDL) {
  CHM_FR L;
  CHM_SP A = AS_CHM_SP__(a);
  double beta[2] = {0, 0};
  beta[0] = asReal(Imult);

  int iSuper = asLogical(super),
      iPerm  = asLogical(perm),
      iLDL   = asLogical(LDL);
  int im     = asInteger(m);
  if ((im < -1) || (im > 3)) error("Cholesky ordering method must be an integer between -1 and 3");

  if (iSuper == NA_LOGICAL)	iSuper = -1;  // NA --> let CHOLMOD choose
  if (iLDL > 0) iSuper = 0;

  c.final_ll = (iLDL == 0) ? 1 : 0;
  c.supernodal = (iSuper > 0) ? CHOLMOD_SUPERNODAL :
    ((iSuper < 0) ? CHOLMOD_AUTO : CHOLMOD_SIMPLICIAL);

  if (iPerm) {
    chm_set_ordering(im);
  } else {  // no permutation, m ignored in this case
    chm_set_ordering(-1);
  }

  L = M_cholmod_analyze(A, &c);

  M_cholmod_factorize_p(A, beta, (int*)NULL, 0 /*fsize*/, L, &c);
  int ok = (L->minor == L->n);
  if (ok) {
    //Rprintf("ok, c.final_ll = %d \n", c.final_ll);
    SEXP out = PROTECT(M_cholmod_factor_as_sexp(L, 0 /* do not free */));
    M_cholmod_free_factor(&L, &c);
    UNPROTECT(1);
    return out;
  } else {
    M_cholmod_free_factor(&L, &c);
    error("Cholesky factorization failed");
    return R_NilValue;
  }
}

// there is no M_chm_dense_to_SEXP in Matrix_stubs so write it here
// based on Matrix package's chm_dense_to_matrix in chm_common.c
SEXP chm_dense_to_matrixSEXP(CHM_DN a) {
  if (a->xtype != CHOLMOD_REAL) error("not a real type cholmod object");
  SEXP ans = PROTECT(allocMatrix(REALSXP, a->nrow, a->ncol));
  Memcpy(REAL(ans), (double *) a->x, a->nrow * a->ncol);
  M_cholmod_free_dense(&a, &c);
  UNPROTECT(1);
  return ans;
}

// basically chm_dense_to_vector, see chm_common.c of Matrix package
SEXP chm_dense_to_vectorSEXP(CHM_DN a) {
  if (a->xtype != CHOLMOD_REAL) error("not a real type cholmod object");
  SEXP ans = PROTECT(allocVector(REALSXP, a->nrow * a->ncol));
  Memcpy(REAL(ans), (double *) a->x, a->nrow * a->ncol);
  M_cholmod_free_dense(&a, &c);
  UNPROTECT(1);
  return ans;
}

// dense vector solve
SEXP CHMf_solve(SEXP a, SEXP b, SEXP system) {
  CHM_FR L = AS_CHM_FR(a);
  int n = LENGTH(b);
  CHM_DN B = N_AS_CHM_DN(REAL(b), n, 1);

  int sys = asInteger(system);

  if (!(sys--)) error("invalid system argument");

  SEXP ans = chm_dense_to_vectorSEXP(M_cholmod_solve(sys, L, B, &c));
  return ans;
}

// dense matrix solve
SEXP CHMf_solve_matrix(SEXP a, SEXP b, SEXP system) {
  CHM_FR L = AS_CHM_FR(a);
  int* dim = INTEGER(getAttrib(b, R_DimSymbol));
  CHM_DN B = N_AS_CHM_DN(REAL(b), dim[0], dim[1]);

  int sys = asInteger(system);

  if (!(sys--)) error("invalid system argument");

  SEXP ans = chm_dense_to_matrixSEXP(M_cholmod_solve(sys, L, B, &c));
  return ans;
}

SEXP CHMf_spsolve(SEXP a, SEXP b, SEXP system) {
  CHM_FR L = AS_CHM_FR(a);
  CHM_SP B = AS_CHM_SP__(b);
  int sys = asInteger(system);

  if (!(sys--)) error("invalid system argument");

  SEXP ans = M_cholmod_sparse_as_sexp(
    M_cholmod_spsolve(sys, L, B, &c), 1, 0, 0, "", R_NilValue);
  return ans;
}

// destructive Cholesky
SEXP CHM_update_inplace(SEXP object, SEXP parent, SEXP mult) {
  CHM_FR L = AS_CHM_FR(object);
  CHM_SP A = AS_CHM_SP__(parent);

  M_cholmod_factor_update(L, A, asReal(mult));
  return R_NilValue;
}
