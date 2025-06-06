% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models.R
\name{correlation}
\alias{correlation}
\alias{iid}
\alias{RW1}
\alias{RW2}
\alias{AR1}
\alias{season}
\alias{spatial}
\alias{splines}
\alias{custom}
\title{Correlation factor structures in generic model components}
\usage{
iid(name)

RW1(name, circular = FALSE, w = NULL)

RW2(name)

AR1(name, phi, w = NULL, control = NULL)

season(name, period)

spatial(
  name,
  graph = NULL,
  snap = sqrt(.Machine$double.eps),
  queen = TRUE,
  poly.df = NULL,
  derive.constraints = FALSE
)

splines(name, knots, degree)

custom(name, D = NULL, Q = NULL, R = NULL, derive.constraints = NULL)
}
\arguments{
\item{name}{name of a variable, unquoted.}

\item{circular}{whether the random walk is circular.}

\item{w}{a vector of weights.}

\item{phi}{prior distribution, or fixed value, for an
autoregressive parameter. The default is a uniform prior over
the interval [-1, 1]. A single numeric value is interpreted as
a fixed value, corresponding to a degenerate prior, which can also
be specified as \code{pr_fixed(value)}. Alternatively,
\code{link{pr_truncnormal}} can be used to specify a truncated
normal prior.}

\item{control}{options for Metropolis-Hastings sampling from the
conditional posterior for an autoregressive parameter. These options
can be set using function \code{\link{set_MH}}. Supported proposal
types are "TN" and "RWTN". By default an independence truncated
normal proposal (type="TN"), or a random walk truncated normal
proposal (type="RWTN") with adaptive scale initialised at 0.025
is used, depending on whether the specified random effects'
distribution is Gaussian or non-Gaussian.}

\item{period}{a positive integer specifying the seasonal period.}

\item{graph}{either a spatial object of class \code{SpatialPolygons},
\code{sf}, \code{sfc}, or a neighbours list of class \code{nb}.}

\item{snap}{passed to \code{\link[spdep]{poly2nb}}.
Ignored if \code{graph} is a neighbours list.}

\item{queen}{passed to \code{\link[spdep]{poly2nb}}.
Ignored if \code{graph} is a neighbours list.}

\item{poly.df}{a spatial data frame. DEPRECATED, use argument
\code{graph} instead.}

\item{derive.constraints}{whether to derive the constraint matrix for an
IGMRF model component numerically from the precision matrix.
The use of \code{derive.constraints} in function \code{spatial}
is DEPRECATED, as it is no longer needed.}

\item{knots}{passed to \code{\link[splines]{splineDesign}}.}

\item{degree}{passed to \code{\link[splines]{splineDesign}}.}

\item{D}{custom incidence matrix.}

\item{Q}{custom precision matrix.}

\item{R}{custom restriction matrix.}
}
\description{
Element 'factor' of a model component created using function
\code{\link{gen}} is a formula composed of several possible terms described
below. It is used to derive a (typically sparse) precision matrix for a set of
coefficients, and possibly a matrix representing a set of linear constraints
to be imposed on the coefficient vector.
\describe{
  \item{iid(f)}{Independent effects corresponding to the levels of factor \code{f}.}
  \item{RW1(f, circular=FALSE, w=NULL)}{First-order random walk over the levels of factor \code{f}.
    The random walk can be made circular and different (fixed) weights can be attached to the innovations.
    If specified, \code{w} must be a positive numeric vector of length one less than the number of
    factor levels. For example, if the levels correspond to different times, it would often be
    reasonable to choose \code{w} proportional to the reciprocal time differences. For equidistant
    times there is generally no need to specify \code{w}.}
  \item{RW2(f)}{Second-order random walk.}
  \item{AR1(f, phi, w=NULL, control=NULL)}{First-order autoregressive correlation structure among
    the levels of \code{f}. Argument \code{phi} can be a single numerical value of the
    autoregressive parameter, or an appropriate prior specification if phi should be inferred.
    If not supplied, a uniform prior on (-1, 1] is assumed.
    For irregularly spaced AR(1) processes weights can be specified, in the same way as for
    \code{RW1}.}
  \item{season(f, period)}{Dummy seasonal with period \code{period}.}
  \item{spatial(f, graph, snap, queen)}{CAR spatial correlation.
    Argument \code{graph} can either be an object of (S4) class \code{SpatialPolygonsDataFrame}
    or an object of (S3) class \code{sf}. The latter can be obtained, e.g., by reading in a
    shape file using function \code{\link[sf]{st_read}}. Arguments \code{snap} and \code{queen}
    are passed to \code{\link[spdep]{poly2nb}}, which computes a neighbours list.
    Alternatively, a neighbours list object of class \code{nb} can be passed directly
    as argument \code{graph}.}
  \item{splines(f, knots, degree)}{P-splines, i.e. penalised B-splines structure over
    the domain of a quantitative variable f. Arguments knots and degree are passed to
    \code{\link[splines]{splineDesign}}. If \code{knots} is a single value it is interpreted as
    the number of knots, otherwise as a vector of knot positions. By default 40 equally spaced
    knots are used, and a degree of 3.}
  \item{custom(f, D=NULL, Q=NULL, R=NULL, derive.constraints=NULL)}{Either a custom precision or incidence
    matrix associated with factor f can be passed to argument \code{Q} or \code{D}. Optionally a
    constraint matrix can be supplied as \code{R}, or constraints can be derived from the null space
    of the precision matrix by setting \code{derive.constraints=TRUE}.}
}
}
\examples{
\donttest{
# example of CAR spatial random effects
if (requireNamespace("sf")) {
  # 1. load a shape file of counties in North Carolina
  nc <- sf::st_read(system.file("shape/nc.shp", package="sf"))
  # 2. generate some data according to a model with a few regression
  # effects, as well as spatial random effects
  gd <- generate_data(
    ~ reg(~ AREA + BIR74, prior=pr_normal(precision=1), name="beta") +
      gen(factor = ~ spatial(NAME, graph=nc), name="vs"),
    family=f_gaussian(var.prior = pr_invchisq(df=10, scale=1)),
    data = nc
  )
  # add the generated target variable and the spatial random effects to the
  # spatial dataframe object
  nc$y <- gd$y
  nc$vs_true <- gd$pars$vs
  # 3. fit a model to the generated data, and see to what extent the
  #    parameters used to generate the data, gd$pars, are reproduced
  sampler <- create_sampler(
    y ~ reg(~ AREA + BIR74, prior=pr_normal(precision=1), name="beta") +
    gen(factor = ~ spatial(NAME, graph=nc), name="vs"),
    data=nc
  )
  # increase burnin and n.iter below to improve MCMC convergence
  sim <- MCMCsim(sampler, store.all=TRUE, burnin=100, n.iter=250, n.chain=2, verbose=FALSE)
  (summ <- summary(sim))
  nc$vs <- summ$vs[, "Mean"]
  plot(nc[c("vs_true", "vs")])
  plot(gd$pars$vs, summ$vs[, "Mean"]); abline(0, 1, col="red")
}
}

}
\references{
B. Allevius (2018).
   On the precision matrix of an irregularly sampled AR(1) process.
   arXiv:1801.03791v2.

 H. Rue and L. Held (2005).
   Gaussian Markov Random Fields.
   Chapman & Hall/CRC.
}
