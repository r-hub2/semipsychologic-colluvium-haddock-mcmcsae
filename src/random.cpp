
#include <Rcpp.h>
#include <GIGrvg.h>
using namespace Rcpp;


// constants used in approximate Polya-Gamma sampling
static const double PI2 = M_PI * M_PI;
static const double EPS = 10.0 * std::numeric_limits<double>::epsilon();


//’ Draw a vector of (approximate) Polya-Gamma variates
//’
//’ @param n the size of the vector.
//’ @param b shape parameter. Either a scalar or vector of nonnegative doubles of size n.
//’ @param z exponential tilting parameter. Either a scalar or vector of doubles of size n.
//’ @param m integer scalar or vector of size n specifying the number of explicit gamma draws
//’  used in the approximation. A -1 value indicates that a normal moment matching approximation
//’  is to be used. For a value less than -1 a default choice for the approximation will be used.
//’ @returns A vector of size n with (approximate) Polya-Gamma draws.
// [[Rcpp::export(rng=true)]]
NumericVector CrPGapprox(const int n, const NumericVector & b, const NumericVector & z, const IntegerVector & m) {
  double bi, hzi, th, mu, Sigma, dninv, rgs;
  int mi;
  const int nb = b.size();
  const int nz = z.size();
  const int nm = m.size();
  NumericVector out = no_init(n);
  for (int i = 0; i < n; ++i) {
    bi = nb == 1 ? b[0] : b[i];
    if (bi < EPS) {
      out[i] = 0;
    } else {
      hzi = nz == 1 ? 0.5*z[0] : 0.5*z[i];
      // compute mean and variance of PG(bi, zi), up to factor bi
      if (std::abs(hzi) < 0.01) {
        mu = 0.25 * (1 - hzi*hzi / 3);
        Sigma = (1 - 0.8 * hzi*hzi) / 24;
      } else {
        th = tanh(hzi);
        mu = 0.25 * th / hzi;
        Sigma = 0.0625 * (th - hzi / std::pow(cosh(hzi), 2)) / std::pow(hzi, 3);
      }
      mi = nm == 1 ? m[0] : m[i];
      if (mi < -1) {
        // default
        if (bi > 200) {
          mi = -1;
        } else if (bi > 20) {
          mi = 0;
        } else if (bi > 2) {
          mi = 1;
        } else if (bi > 0.5) {
          mi = 2;
        } else {
          mi = 4;
        }
      }
      // draw from approximation to PG(b, z)
      switch(mi) {
      case -1:
        out[i] = R::rnorm(bi*mu, std::sqrt(bi*Sigma));
        break;
      case 0:
        out[i] = R::rgamma(bi*mu*mu/Sigma, Sigma/mu);
        break;
      case 1:
        dninv = 2 / (0.25*PI2 + hzi*hzi);
        mu -= 0.25 * dninv;
        Sigma -= 0.0625 * dninv*dninv;
        out[i] = 0.25 * dninv * R::rgamma(bi, 1) + R::rgamma(bi*mu*mu/Sigma, Sigma/mu);
        break;
      default:
        rgs = 0;
        for (int j = 0; j < mi; ++j) {
          dninv = 2 / (PI2 * std::pow(j + 0.5, 2) + hzi*hzi);
          rgs += 0.25 * dninv * R::rgamma(bi, 1);
          mu -= 0.25 * dninv;
          Sigma -= 0.0625 * dninv*dninv;
        }
        out[i] = rgs + R::rgamma(bi*mu*mu/Sigma, Sigma/mu);
      }
    }
  }
  return out;
}


//’ Draw a vector of normal variates
//’
//’ @param n the size of the vector.
//’ @param mean scalar mean.
//’ @param sd scalar standard deviation.
//’ @returns A vector of size n with draws from a normal distribution.
// [[Rcpp::export(rng=true)]]
NumericVector Crnorm(const int n, const double mean = 0, const double sd = 1) {
  return rnorm(n, mean, sd);
}


// (internal) wrapper function for do_rgig in GIGrvg package
double do_rgig1(double lambda, double chi, double psi) {
  SEXP (*fun)(int, double, double, double) = NULL;
  if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
  return as<double>(fun(1, lambda, chi, psi));
}

//’ Draw a vector of generalized inverse gaussian (GiG) variates
//’
//’ @param n the size of the vector.
//’ @param p (vector of) shape parameters.
//’ @param a (vector of) shape/scale parameters.
//’ @param b (vector of) shape/scale parameters.
//’ @returns A vector of size n with draws from a GiG distribution.
// [[Rcpp::export(rng=true)]]
NumericVector Crgig(const int n, const NumericVector & p, const NumericVector & a, const NumericVector & b) {
  NumericVector out = no_init(n);
  const int np=p.size();
  const int na=a.size();
  const int nb=b.size();
  //double pi,ai,bi;
  for (int i = 0; i < n; ++i) {
    /*
    NB issue has been solved in GiGrvg 0.7
    // edge case issue in GIGrvg; for now we deal with these (gamma/invgamma) cases ourselves
    pi = np == 1 ? p[0] : p[i];
    ai = na == 1 ? a[0] : a[i];
    bi = nb == 1 ? b[0] : b[i];
    if (ai < EPS || bi < EPS) {
      if (pi > 0.0) {
        out[i] = R::rgamma(pi, 2.0/ai);
      } else {
        out[i] = 1.0/R::rgamma(-pi, 2.0/bi);
      }
    } else {
      out[i] = do_rgig1(pi, bi, ai);
    }
    */
    // parameter translation: lambda=p, chi=b, psi=a
    out[i] = do_rgig1(
      np == 1 ? p[0] : p[i],
      nb == 1 ? b[0] : b[i],
      na == 1 ? a[0] : a[i]
    );
  }
  return out;
}


//’ Draw a vector of (approximate) Chinese Restaurant Table (CRT) variates
// Used in a Gibbs sampler for negative binomial model with modeled shape parameter.
// The approximation is based on Le Cam's theorem, i.e. the approximation of a convolution
// of Bernoulli random variables by a Poisson distribution. The sampling is exact for all
// values of \code{y} less than or equal to \code{2*m}.
//’
//’ @param y data vector.
//’ @param r (inverse) dispersion or shape parameter, can be scalar or vector.
//’ @param m positive integer; larger values give more accuracy but slower performance.
//’ @returns A vector of (approximate) CRT variates.
// [[Rcpp::export(rng=true)]]
IntegerVector CrCRT(const NumericVector & y, const NumericVector & r, const int m=20) {
  double prob, lambda;
  int m_expl;
  const int two_m = 2 * m;
  const int n = y.size();
  const int nr = r.size();
  double ri = r[0];
  IntegerVector out(n);
  for (int i = 0; i < n; i++) {
    if (nr > 1) ri = r[i];
    if (y[i] <= two_m) {
      // exact CRT sampling
      for (int j = 0; j < y[i]; j++) {
        prob = ri / (ri + j);
        if (R::runif(0, 1) < prob) out[i]++;
      }
    } else {
      m_expl = std::min(m, (int)ri);
      // first m_expl Bernoulli draws
      for (int j = 0; j < m_expl; j++) {
        prob = ri / (ri + j);
        if (R::runif(0, 1) < prob) out[i]++;
      }
      // then approximate remaining y[i] - m_expl draws
      lambda = ri * (R::digamma(y[i] + ri) - R::digamma(m_expl + ri));
      out[i] += R::rpois(lambda);
    }
  }
  return out;
}
