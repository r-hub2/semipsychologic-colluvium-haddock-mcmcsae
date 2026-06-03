
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
  double bi, hzi, hzi2, th, mu, Sigma, dninv, rgs;
  int mi;
  const int nb = b.size();
  const int nz = z.size();
  const int nm = m.size();
  NumericVector out = no_init(n);
  for (int i = 0; i < n; ++i) {
    bi = nb == 1 ? b[0] : b[i];
    if (bi < EPS) {
      out[i] = 0.0;
    } else {
      hzi = nz == 1 ? 0.5*z[0] : 0.5*z[i];
      hzi2 = hzi * hzi;
      // compute mean and variance of PG(bi, zi), up to factor bi
      if (std::abs(hzi) < 0.01) {
        mu = 0.25 * (1.0 - hzi2 / 3.0);
        Sigma = (1.0 - 0.8 * hzi2) / 24.0;
      } else {
        th = tanh(hzi);
        mu = 0.25 * th / hzi;
        Sigma = 0.0625 * (th - hzi * (1 - th * th)) / (hzi2 * hzi);
      }
      mi = nm == 1 ? m[0] : m[i];
      if (mi < -1) {
        // default
        if (bi > 200.0) {
          mi = -1;
        } else if (bi > 20.0) {
          mi = 0;
        } else if (bi > 2.0) {
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
          dninv = 2.0 / (0.25*PI2 + hzi2);
          mu -= 0.25 * dninv;
          Sigma -= 0.0625 * dninv*dninv;
          out[i] = 0.25 * dninv * R::rgamma(bi, 1) + R::rgamma(bi*mu*mu/Sigma, Sigma/mu);
          break;
        default:
          double rgs = 0.0;
          double jhalf;
          for (int j = 0; j < mi; ++j) {
            jhalf = j + 0.5;
            dninv = 2.0 / (PI2 * jhalf * jhalf + hzi2);
            rgs += 0.25 * dninv * R::rgamma(bi, 1);
            mu -= 0.25 * dninv;
            Sigma -= 0.0625 * dninv*dninv;
          }
          out[i] = rgs + R::rgamma(bi*mu*mu/Sigma, Sigma/mu);
      }  // END switch(mi)
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
  const int n = y.size();
  const int nr = r.size();
  const int two_m = 2 * m;
  IntegerVector out(n);
  const double* p_y = y.begin();
  const double* p_r = r.begin();
  int* p_out = out.begin();
  if (nr == 1) {  // scalar r
    const double ri = p_r[0];
    const int m_expl = std::min(m, (int)ri);
    const double digamma_m_ri = R::digamma(m_expl + ri);
    std::vector<double> prob_cache(two_m + 1);
    for (int j = 0; j <= two_m; j++) {
      prob_cache[j] = ri / (ri + j);
    }
    for (int i = 0; i < n; i++) {
      const double yi = p_y[i];
      int count = 0;
      if (yi <= two_m) {
        // exact CRT sampling
        for (int j = 0; j < yi; j++) {
          if (R::runif(0, 1) < prob_cache[j]) count++;
        }
      } else {
        // first m_expl Bernoulli draws
        for (int j = 0; j < m_expl; j++) {
          if (R::runif(0, 1) < prob_cache[j]) count++;
        }
        // then approximate remaining y[i] - m_expl draws
        double lambda = ri * (R::digamma(yi + ri) - digamma_m_ri);
        count += R::rpois(lambda);
      }
      p_out[i] = count;
    }
  } else {  // vector r
    for (int i = 0; i < n; i++) {
      const double yi = p_y[i];
      const double ri = p_r[i];
      int count = 0;
      if (yi <= two_m) {
        for (int j = 0; j < yi; j++) {
          if (R::runif(0, 1) < (ri / (ri + j))) count++;
        }
      } else {
        const int m_expl = std::min(m, (int)ri);
        for (int j = 0; j < m_expl; j++) {
          if (R::runif(0, 1) < (ri / (ri + j))) count++;
        }
        double lambda = ri * (R::digamma(yi + ri) - R::digamma(m_expl + ri));
        count += R::rpois(lambda);
      }
      p_out[i] = count;
    }
  }
  return out;
}
