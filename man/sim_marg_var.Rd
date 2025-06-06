% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cMVN_sampler.R
\name{sim_marg_var}
\alias{sim_marg_var}
\title{Compute a Monte Carlo estimate of the marginal variances of a (I)GMRF}
\usage{
sim_marg_var(
  D,
  Q = NULL,
  R = NULL,
  r = NULL,
  eps1 = 1e-09,
  eps2 = 1e-09,
  nSim = 100L
)
}
\arguments{
\item{D}{factor of precision matrix Q such that Q=D'D.}

\item{Q}{precision matrix.}

\item{R}{equality restriction matrix.}

\item{r}{rhs vector for equality constraints \eqn{R'x = r}, where \eqn{R'} denotes the transpose of R.}

\item{eps1}{passed to \code{\link{create_cMVN_sampler}}.}

\item{eps2}{passed to \code{\link{create_cMVN_sampler}}.}

\item{nSim}{number of Monte Carlo samples used to estimate the marginal variances.}
}
\value{
A vector of Monte Carlo estimates of the marginal variances.
}
\description{
Estimate marginal variances of a (I)GMRF prior defined in terms
of a sparse precision matrix and possibly a set of equality constraints.
The marginal variances might be used to rescale the precision matrix
such that a default prior for a corresponding variance component is
more appropriate.
}
\references{
S.H. Sorbye and H. Rue (2014).
   Scaling intrinsic Gaussian Markov random field priors in spatial modelling.
   Spatial Statistics, 8, 39-51.
}
