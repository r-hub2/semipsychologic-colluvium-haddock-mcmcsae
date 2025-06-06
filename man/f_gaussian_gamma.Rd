% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_gamma.R
\name{f_gaussian_gamma}
\alias{f_gaussian_gamma}
\title{Specify a Gaussian-Gamma sampling distribution}
\usage{
f_gaussian_gamma(link = "identity", var.model, ...)
}
\arguments{
\item{link}{the name of a link function. Currently the only allowed link function
for this distribution family is \code{"identity"}.}

\item{var.model}{a formula specifying the terms of a variance model.
The left-hand side of the formula should specify the observed variances,
unless the family object is used for data generation only.
Several types of model terms on the right-hand side of the formula are supported:
a regression term for the log-variance specified with \code{\link{vreg}(...)},
and a term \code{\link{vfac}(...)} for multiplicative modelled factors
at a certain level specified by a factor variable. In addition, \code{\link{reg}} and \code{\link{gen}}
can be used to specify regression or random effect terms. In that case the prior distribution
of the coefficients is not exactly normal, but instead Multivariate Log inverse Gamma (MLiG),
see also \code{\link{pr_MLiG}}.}

\item{...}{further arguments passed to \code{\link{f_gamma}}.}
}
\value{
A family object.
}
\description{
This function can be used in the \code{family} argument of
\code{\link{create_sampler}} or \code{\link{generate_data}} to specify a
Gaussian-Gamma sampling distribution, i.e., a Gaussian sampling distribution
whose variances are observed subject to error according to a Gamma
distribution.
}
