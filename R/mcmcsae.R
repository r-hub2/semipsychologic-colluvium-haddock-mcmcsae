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
#' @importFrom Matrix .diag2sparse .m2sparse .M2m .updateCHMfactor bandSparse
#' @importFrom Matrix bdiag coerce diag Diagonal drop0 expand1 forceSymmetric
#' @importFrom Matrix invertPerm isDiagonal KhatriRao Matrix nnzero
#' @importFrom Matrix rsparsematrix sparseMatrix
#' @importClassesFrom Matrix CHMfactor dCHMsimpl ddiMatrix CsparseMatrix
#' @importClassesFrom Matrix dgCMatrix dsCMatrix generalMatrix sparseMatrix
#' @importMethodsFrom Matrix %*% as.matrix as.vector Cholesky colSums crossprod
#' @importMethodsFrom Matrix determinant diag isSymmetric rowSums solve t tcrossprod unname
#' @import GIGrvg
#' @importFrom collapse allNA allv any_duplicated anyv dapply fdroplevels fmatch
#' @importFrom collapse fmean.default fmean.matrix fquantile fsd.matrix fsum.matrix
#' @importFrom collapse fvar.matrix qF whichNA whichv
#' @importFrom graphics abline axis legend lines matplot pairs par plot
#' @importFrom graphics plot.new points segments
#' @importFrom methods as cbind2 new rbind2 setAs setClass setMethod signature show
#' @importFrom stats acf as.formula density fitted make.link mvfft nextn optim
#' @importFrom stats pnorm predict rbeta rbinom rchisq residuals rexp rgamma rnbinom rnorm
#' @importFrom stats rpois runif rWishart sd setNames terms update.formula var weights
#' @importFrom utils getFromNamespace modifyList object.size setTxtProgressBar
#' @importFrom utils str tail txtProgressBar
NULL
