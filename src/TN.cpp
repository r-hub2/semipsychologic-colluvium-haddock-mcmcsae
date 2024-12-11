
// RcppEigen.h also pulls in Rcpp.h
#include <RcppEigen.h>
using namespace Rcpp;


// the code for univariate truncated standard normal sampling is based on
// MatLab code by Z.I. Botev, see https://web.maths.unsw.edu.au/~zdravkobotev/
// and R package TruncatedNormal


static const double A = 0.4;
static const double minusA = -0.4;
static const double B = 2.05;
static const double epsHMC = std::sqrt(std::numeric_limits<double>::epsilon());


// Rayleigh rejection sampling
double nt(const double l, const double u) {
  double x;
  const double c = 0.5*l*l;
  const double f = std::expm1(c - 0.5*u*u);
  do {
    x = c - std::log1p(f * R::runif(0, 1));
  } while (x * std::pow(R::runif(0, 1), 2) > c);
  return std::sqrt(2*x);
}

// simple rejection sampling
double trnd(const double l, const double u) {
  double x;
  do {
    x = R::rnorm(0, 1);
  } while (x < l || x > u);
  return x;
}

//’ Generate a random value from a standardized univariate truncated normal distribution
//’
//’ @param l lower truncation bound.
//’ @param u upper truncation bound.
//’ @returns A single draw from the standard univariate truncated normal distribution.
// [[Rcpp::export(rng=true)]]
double Crtuvn(const double l, const double u) {
  double out;
  if (l > A) {
    out = nt(l, u);
  } else if (u < minusA) {
    out = -nt(-u, -l);
  } else {
    if (std::abs(u - l) > B) {
      out = trnd(l, u);
    } else {
      out = R::qnorm(
        R::pnorm(l, 0, 1, true, false) + 
          (R::pnorm(u, 0, 1, true, false) - R::pnorm(l, 0, 1, true, false)) * R::runif(0, 1),
        0, 1, true, false
      );
    }
  }
  return out;
}

//’ Generate a Gibbs cycle for a standardized multivariate truncated normal distribution with dense constraint matrix
//’
//’ @param v start value.
//’ @param Ut dense matrix encoding the inequality constrained linear combinations.
//’ @param ustar vector encoding the (initial) lower bounds corresponding to Ut.
//’ @param eps small positive value to control the numerical stability of the Gibbs sampler.
//’ @returns A single draw from the standardized multivariate truncated normal distribution.
//’ @references
//’  Y. Li and S.K. Ghosh (2015). Efficient sampling methods for truncated multivariate normal
//’    and student-t distributions subject to linear inequality constraints.
//’    Journal of Statistical Theory and Practice 9(4), 712-732.
// [[Rcpp::export(rng=true)]]
Eigen::VectorXd Crtmvn_Gibbs_dense(const Eigen::Map<Eigen::VectorXd> & v, const Eigen::Map<Eigen::MatrixXd> Ut,
    const Eigen::Map<Eigen::VectorXd> & ustar, const double eps) {
  double a, b, vi, x, temp;
  Eigen::VectorXd u(ustar);
  int n = v.size();
  int m = u.size();
  Eigen::VectorXd out(n);
  // loop over variables
  for (int i = 0; i < n; i++) {
    a = R_NegInf;
    b = R_PosInf;
    vi = v[i];
    // loop over constraints that variable i is involved in
    for (int j = 0; j < m; j++) {
      x = Ut(j,i);
      u[j] += x * vi;
      if (x > eps) {
        temp = u[j]/x;
        if (temp > a) a = temp;
      } else if (x < -eps) {
        temp = u[j]/x;
        if (temp < b) b = temp;
      }
    }
    if (a == R_NegInf && b == R_PosInf) {
      out[i] = R::rnorm(0, 1);
    } else if (a == b) {
      out[i] = a;
    } else if (a < b) {
      out[i] = Crtuvn(a, b);
    } else {
      // this seems a numerically stable way to deal with numerical inaccuracy:
      if (a < vi) {
        out[i] = a;
      } else if (b > vi) {
        out[i] = b;
      } else {
        out[i] = vi;
      }
    }
    // kan in de laatste iteratie overslaan:
    u -= out[i] * Ut.col(i);
  }
  return out;
}

//’ Generate a Gibbs cycle for a standardized multivariate truncated normal distribution with sparse constraint matrix
//’
//’ @param v start value.
//’ @param Ut sparse matrix encoding the inequality constrained linear combinations.
//’ @param ustar vector encoding the (initial) lower bounds corresponding to Ut.
//’ @param eps small numerical value to stabilize the Gibbs sampler.
//’ @returns A single draw from the standardized multivariate truncated normal distribution.
// [[Rcpp::export(rng=true)]]
NumericVector Crtmvn_Gibbs_sparse(const NumericVector & v, const SEXP Ut, const NumericVector & ustar, const double eps) {
  double a, b, vi, x, temp;
  if (!Rf_isS4(Ut) || !Rf_inherits(Ut, "dgCMatrix")) stop("Ut is not a dgCMatrix");
  const IntegerVector Utp(as<S4>(Ut).slot("p"));
  const IntegerVector Uti(as<S4>(Ut).slot("i"));
  const NumericVector Utx(as<S4>(Ut).slot("x"));
  NumericVector u(clone(ustar));
  int n = v.size();
  NumericVector out = no_init(n);
  // loop over variables
  for (int i = 0; i < n; i++) {
    a = R_NegInf;
    b = R_PosInf;
    vi = v[i];
	for (int j = Utp[i]; j < Utp[i + 1]; j++) {
      x = Utx[j];
      u[Uti[j]] += x * vi;
      if (x > eps) {
        temp = u[Uti[j]]/x;
        if (temp > a) a = temp;
      } else if (x < -eps) {
        temp = u[Uti[j]]/x;
        if (temp < b) b = temp;
      }
    }
    if (a == R_NegInf && b == R_PosInf) {
      out[i] = R::rnorm(0, 1);
    } else if (a == b) {
      out[i] = a;
    } else if (a < b) {
      out[i] = Crtuvn(a, b);
    } else {
      // this seems a numerically stable way to deal with numerical inaccuracy:
      if (a < vi) {
        out[i] = a;
      } else if (b > vi) {
        out[i] = b;
      } else {
        out[i] = vi;
      }
    }
    for (int j = Utp[i]; j < Utp[i + 1]; j++) {
      u[Uti[j]] -= Utx[j] * out[i];
    }
  }
  return out;
}  

//’ Use a Gibbs within slice method to generate a next draw from a standardized multivariate truncated normal distribution with dense constraint matrix
//’
//’ @param v start value.
//’ @param Ut dense matrix encoding the inequality constrained linear combinations.
//’ @param ustar vector encoding the (initial) lower bounds corresponding to Ut.
//’ @param eps small numerical value to stabilize the Gibbs sampler.
//’ @returns A single draw from the standardized multivariate truncated normal distribution.
//’ @references
//’  Y. Li and S.K. Ghosh (2015). Efficient sampling methods for truncated multivariate normal
//’    and student-t distributions subject to linear inequality constraints.
//’    Journal of Statistical Theory and Practice 9(4), 712-732.
//’  K.A. Valeriano, C.E. Galarza and L.A. Matos (2023). Moments and random number generation
//’    for the truncated elliptical family of distributions.
//’    Statistics and Computing 33(1), 1-20.
// [[Rcpp::export(rng=true)]]
Eigen::VectorXd Crtmvn_slice_Gibbs_dense(const Eigen::Map<Eigen::VectorXd> & v, const Eigen::Map<Eigen::MatrixXd> Ut,
                                   const Eigen::Map<Eigen::VectorXd> & ustar, const double eps) {
  double a, b, vi, x, temp;
  Eigen::VectorXd u(ustar);
  int n = v.size();
  int m = u.size();
  double xx = v.dot(v);  // not necessary to recompute in every iteration, but may be more numerically stable
  double y = R::runif(0, std::exp(-0.5*xx));
  double kappa_y = -2*std::log(y);
  Eigen::VectorXd out(n);
  // loop over variables
  for (int i = 0; i < n; i++) {
    vi = v[i];
    xx -= vi*vi;
    b = std::sqrt(kappa_y - xx);
    a = -b;
    // loop over constraints that variable i is involved in
    for (int j = 0; j < m; j++) {
      x = Ut(j,i);
      u[j] += x * vi;
      if (x > eps) {
        temp = u[j]/x;
        if (temp > a) a = temp;
      } else if (x < -eps) {
        temp = u[j]/x;
        if (temp < b) b = temp;
      }
    }
    if (a == b) {
      out[i] = a;
    } else if (a < b) {
      out[i] = R::runif(a, b);
    } else {
      // this seems a numerically stable way to deal with numerical inaccuracy:
      if (a < vi) {
        out[i] = a;
      } else if (b > vi) {
        out[i] = b;
      } else {
        out[i] = vi;
      }
    }
    // kan in de laatste iteratie overslaan:
    u -= out[i] * Ut.col(i);
    xx += out[i] * out[i];
  }
  return out;
}

//’ Use a Gibbs within slice method to generate a next draw from a standardized multivariate truncated normal distribution with sparse constraint matrix
//’
//’ @param v start value.
//’ @param Ut sparse matrix encoding the inequality constrained linear combinations.
//’ @param ustar vector encoding the (initial) lower bounds corresponding to Ut.
//’ @param eps small numerical value to stabilize the Gibbs sampler.
//’ @returns A single draw from the standardized multivariate truncated normal distribution.
// [[Rcpp::export(rng=true)]]
NumericVector Crtmvn_slice_Gibbs_sparse(const NumericVector & v, const SEXP Ut, const NumericVector & ustar, const double eps) {
  double a, b, vi, x, temp;
  if (!Rf_isS4(Ut) || !Rf_inherits(Ut, "dgCMatrix")) stop("Ut is not a dgCMatrix");
  const IntegerVector Utp(as<S4>(Ut).slot("p"));
  const IntegerVector Uti(as<S4>(Ut).slot("i"));
  const NumericVector Utx(as<S4>(Ut).slot("x"));
  NumericVector u(clone(ustar));
  int n = v.size();
  // not necessary to recompute xx in every iteration, but may be more numerically stable
  double xx = 0;
  for (int i = 0; i < v.size(); i++) {
    xx += v[i] * v[i];
  }
  double y = R::runif(0, std::exp(-0.5*xx));
  double kappa_y = -2*std::log(y);
  NumericVector out = no_init(n);
  // loop over variables
  for (int i = 0; i < n; i++) {
    vi = v[i];
    xx -= vi*vi;
    b = std::sqrt(kappa_y - xx);
    a = -b;
    for (int j = Utp[i]; j < Utp[i + 1]; j++) {
      x = Utx[j];
      u[Uti[j]] += x * vi;
      if (x > eps) {
        temp = u[Uti[j]]/x;
        if (temp > a) a = temp;
      } else if (x < -eps) {
        temp = u[Uti[j]]/x;
        if (temp < b) b = temp;
      }
    }
    if (a == b) {
      out[i] = a;
    } else if (a < b) {
      out[i] = R::runif(a, b);
    } else {
      // this seems a numerically stable way to deal with numerical inaccuracy:
      if (a < vi) {
        out[i] = a;
      } else if (b > vi) {
        out[i] = b;
      } else {
        out[i] = vi;
      }
    }
    for (int j = Utp[i]; j < Utp[i + 1]; j++) {
      u[Uti[j]] -= Utx[j] * out[i];
    }
    xx += out[i] * out[i];
  }
  return out;
}

//’ Generate a vector of truncated normal variates for use in probit binomial model
//’
//’ @param mu vector of normal means.
//’ @param y response vector.
//’ @returns A vector of truncated normal variates.
// [[Rcpp::export(rng=true)]]
NumericVector CrTNprobit(const NumericVector & mu, const NumericVector & y) {
  double l, u;
  const int n = mu.size();
  NumericVector out = no_init(n);
  for (int i = 0; i < n; i++) {
    l = y[i] == 1 ? -mu[i] : R_NegInf;
    u = y[i] == 0 ? -mu[i] : R_PosInf;
    out[i] = mu[i] + Crtuvn(l, u);
  }
  return out;
}

//’ Simulate Hamiltonian dynamics for a multivariate truncated normal distribution
//’
//’ @param S inequality constraint matrix, referencing either matrix or dgCMatrix
//’ @param S_cols number of coulmns of S, i.e. number of inequalities
//’ @param v0 generated starting velocity.
//’ @param x0 current state.
//’ @param s_adj adjusted rhs vector corresponding to S.
//’ @param refl_fac precomputed vector for reflection.
//’ @param zero_mu whether mean vector mu is all zeros.
//’ @param mu mean vector; ignored if zero_mu is true.
//’ @param simplified true in case of identity covariance matrix and no equality constraints.
//’ @param VS standardized inequality matrix S; ignored if simplified is true.
//’ @param diagnostic whether to keep track of the number of bounces.
//’ @param bounces number of bounces off each of the inequality walls.
//’ @param t_sim simulation time.
//’ @param max_refl break off the simulation if there are more than so many bounces.
//’ @returns New state x. If diagnostic, then bounces are updated.
//’ @references
//’  A. Pakman and L. Paninski (2014).
//’    Exact Hamiltonian Monte Carlo for truncated multivariate gaussians.
//’    Journal of Computational and Graphical Statistics 23(2), 518-542.
// [[Rcpp::export(rng=false)]]
Eigen::VectorXd TMVN_HMC_C(
    const SEXP S,
    const int S_cols,  // number of columns of S, i.e. number of inequalities
    const Eigen::Map<Eigen::VectorXd> & v0,
    const Eigen::Map<Eigen::VectorXd> & x0,
    const Eigen::Map<Eigen::VectorXd> & s_adj,
    const Eigen::Map<Eigen::VectorXd> & refl_fac,
    const bool zero_mu,
    const Eigen::Map<Eigen::VectorXd> & mu,
    const bool simplified,
    const SEXP VS,
    const bool diagnostic,
    Eigen::Map<Eigen::VectorXi> & bounces,
    double t_sim,
    int max_refl) {

  double cost, sint, vproj;

  Eigen::VectorXd vtemp(v0.size());

  Eigen::VectorXd v = v0;
  Eigen::VectorXd x = x0;

  Eigen::VectorXd Sv(S_cols);
  Eigen::VectorXd Sx(S_cols);

  int n_refl = 0;

  while (true) {
    if (Rf_isMatrix(S)) {
      Eigen::Map<Eigen::MatrixXd> Sm = as<Eigen::Map<Eigen::MatrixXd> >(S);
      Sv.noalias() = Sm.transpose() * v;
      if (zero_mu) {
        Sx.noalias() = Sm.transpose() * x;
      } else {
        Sx.noalias() = Sm.transpose() * (x - mu);
      }
    } else if (Rf_inherits(S, "dgCMatrix")) {
      Eigen::MappedSparseMatrix<double> Sm = as<Eigen::MappedSparseMatrix<double> >(S);
      Sv.noalias() = Sm.transpose() * v;
      if (zero_mu) {
        Sx.noalias() = Sm.transpose() * x;
      } else {
        Sx.noalias() = Sm.transpose() * (x - mu);
      }
    } else {
      stop("unexpected matrix type");
    }
    Eigen::VectorXd u = sqrt(Sv.array().square() + Sx.array().square());

    int ind_first_hit = -1;
    double t_first_hit = std::numeric_limits<double>::infinity();

    int n_walls = u.size();
    double u_inv, phi, t_hit1, t_hit2;
    for (int i = 0; i < n_walls; i++) {
      if (u[i] <= abs(s_adj[i])) continue;
      u_inv = 1/u[i];
      phi = ((Sv[i] <= 0) - (Sv[i] >= 0)) * std::acos(Sx[i] * u_inv);
      t_hit1 = std::acos(s_adj[i] * u_inv) - phi;
      t_hit2 = -(t_hit1 + 2 * phi);
      if (t_hit1 < -epsHMC) t_hit1 += 2 * M_PI;
      if (t_hit1 < epsHMC) t_hit1 = 0;
      if (t_hit1 < t_first_hit) {
        if (t_hit1 > 0 || Sv[i] <= 0) {
          t_first_hit = t_hit1;
          ind_first_hit = i;
        }
      }
      if (t_hit2 < -epsHMC) t_hit2 += 2 * M_PI;
      if (t_hit2 < epsHMC) t_hit2 = 0;
      if (t_hit2 < t_first_hit) {
        if (t_hit2 > 0 || Sv[i] <= 0) {
          t_first_hit = t_hit2;
          ind_first_hit = i;
        }
      }
    }
    
    if (t_first_hit >= t_sim) break;
    
    // next wall hit happens within simulation period
    if (t_first_hit > 0) {
      if (n_refl > max_refl) {
        t_first_hit = R::runif(0, t_first_hit);  // break somewhere before next bounce
      }
      cost = cos(t_first_hit);
      sint = sin(t_first_hit);
      vtemp = v;
      if (zero_mu) {
        v = cost * vtemp - sint * x;
        x = sint * vtemp + cost * x;
      } else {
        v = cost * vtemp - sint * (x - mu);
        x = mu + sint * vtemp + cost * (x - mu);
      }
      t_sim -= t_first_hit;
    }
    if (n_refl > max_refl) break;
    
    // reflect momentum
    if (Rf_isMatrix(S)) {
      Eigen::Map<Eigen::MatrixXd> Sm = as<Eigen::Map<Eigen::MatrixXd> >(S);
      vproj = Sm.col(ind_first_hit).dot(v);
    } else {
      Eigen::MappedSparseMatrix<double> Sm = as<Eigen::MappedSparseMatrix<double> >(S);
      vproj = Sm.col(ind_first_hit).dot(v);
    }
    if (vproj < 0) {  // a bounce
      double alpha = vproj * refl_fac[ind_first_hit];
      if (simplified) {
        if (Rf_isMatrix(S)) {
          Eigen::Map<Eigen::MatrixXd> Sm = as<Eigen::Map<Eigen::MatrixXd> >(S);
          v -= alpha * Sm.col(ind_first_hit);
        } else if (Rf_inherits(S, "dgCMatrix")) {
          Eigen::MappedSparseMatrix<double> Sm = as<Eigen::MappedSparseMatrix<double> >(S);
          v -= alpha * Sm.col(ind_first_hit);
        }
      } else {
        if (Rf_isMatrix(VS)) {
          Eigen::Map<Eigen::MatrixXd> VSm = as<Eigen::Map<Eigen::MatrixXd> >(VS);
          v -= alpha * VSm.col(ind_first_hit);
        } else if (Rf_inherits(VS, "dgCMatrix")) {
          Eigen::MappedSparseMatrix<double> VSm = as<Eigen::MappedSparseMatrix<double> >(VS);
          v -= alpha * VSm.col(ind_first_hit);
        }
      }
      if (diagnostic) {
        bounces[ind_first_hit]++;
      }
      n_refl++;
    }  // else do nothing: passage to the feasible side of the inequality wall
  }  // END while (true)

  if (t_sim < epsHMC) return x;
  if (zero_mu) {
    return sin(t_sim) * v + cos(t_sim) * x;
  } else {
    return mu + sin(t_sim) * v + cos(t_sim) * (x - mu);
  }
}
