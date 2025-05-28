#' Markov Chain Monte Carlo Small Area Estimation
#'
#' Fit multi-level models with possibly correlated random effects using MCMC.
#'
#' Functions to fit multi-level models with Gaussian, binomial, multinomial,
#' negative binomial or Poisson likelihoods using MCMC. Models with a linear predictor
#' consisting of various possibly correlated random effects are supported, allowing
#' flexible modelling of temporal, spatial or other kinds of dependence structures.
#' For Gaussian models the variance can be modelled too. By modelling variances
#' at the unit level the marginal distribution can be changed to a Student-t or Laplace
#' distribution, which may account better for outliers.
#' The package has been developed with applications to small area estimation
#' in official statistics in mind. The posterior samples for the model
#' parameters can be passed to a prediction function to generate samples from
#' the posterior predictive distribution for user-defined quantities such as
#' finite population domain means. For model assessment, posterior predictive
#' checks and DIC/WAIC criteria can easily be computed.
#'
#' @name mcmcsae-package
#' @aliases mcmcsae
NULL

#' @importFrom Rcpp evalCpp
#' @useDynLib mcmcsae, .registration=TRUE
NULL

# other namespace imports
#' @importFrom Matrix .m2sparse .updateCHMfactor bandSparse bdiag coerce
#'  diag Diagonal drop0 expand1 forceSymmetric invertPerm isDiagonal KhatriRao
#'  Matrix nnzero rsparsematrix sparseMatrix
#' @importClassesFrom Matrix CHMfactor dCHMsimpl ddiMatrix CsparseMatrix
#'  dgCMatrix dsCMatrix generalMatrix sparseMatrix
#' @importMethodsFrom Matrix %*% as.matrix as.vector Cholesky colSums crossprod
#'  determinant diag isSymmetric rowSums solve t tcrossprod unname
## do not import which() S4 generic from Matrix package as it slows down normal use of which
## @rawNamespace import(Matrix, except = which)
#' @import GIGrvg
#' @importFrom collapse allNA allv any_duplicated anyv dapply fdroplevels fmatch
#'  fmean.default fmean.matrix fquantile fsd.matrix fsum.matrix fvar.matrix qF
#'  whichNA whichv
#' @importFrom graphics abline axis legend lines matplot pairs par plot
#'  plot.new points segments
#' @importFrom methods as cbind2 new rbind2 setAs setClass setMethod signature
#'  show
#' @importFrom stats acf as.formula density fitted make.link mvfft optim
#'  pnorm predict rbeta rbinom rchisq residuals rexp rgamma rnbinom rnorm
#'  rpois runif rWishart sd setNames terms update.formula var weights
#' @importFrom utils modifyList object.size setTxtProgressBar str tail
#'  txtProgressBar
NULL
