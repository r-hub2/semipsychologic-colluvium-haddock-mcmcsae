## ----echo=FALSE-----------------------------------------------------------------------------------
options(width=100)  # width of output

## ----message=FALSE--------------------------------------------------------------------------------
library(survey)
data(api)
# define the regression model
model <- api00 ~ ell + meals + stype + hsg + col.grad + grad.sch
# compute corresponding population totals
XpopT <- colSums(model.matrix(model, apipop))
N <- XpopT[["(Intercept)"]]  # population size
# create the survey design object
des <- svydesign(ids=~1, data=apisrs, weights=~pw, fpc=~fpc)
# compute the calibration or GREG estimator
cal <- calibrate(des, formula=model, population=XpopT)
svymean(~ api00, des)  # equally weighted estimate
svymean(~ api00, cal)  # GREG estimate

## -------------------------------------------------------------------------------------------------
mean(apipop$api00)

## ----message=FALSE--------------------------------------------------------------------------------
library(mcmcsae)
set.seed(1)
sampler <- create_sampler(model, data=apisrs)
sim <- MCMCsim(sampler, verbose=FALSE)
(summ <- summary(sim))
compute_DIC(sim)

## -------------------------------------------------------------------------------------------------
m <- match(apisrs$cds, apipop$cds)  # population units in the sample
# use only a sample of 250 draws from each chain
predictions <- predict(
  sim, newdata=apipop[-m, ], iters=sample(1:1000, 250),
  show.progress=FALSE
)
str(predictions)
samplesum <- sum(apisrs$api00)
summary(transform_dc(
  predictions, fun = function(x) (samplesum + sum(x))/N
))

## -------------------------------------------------------------------------------------------------
summary(predict(
  sim, newdata=apipop[-m, ],
  fun=function(x, p) (samplesum + sum(x))/N,
  show.progress=FALSE
))

## -------------------------------------------------------------------------------------------------
n <- nrow(apisrs)
XsamT <- colSums(model.matrix(model, apisrs))
XpopR <- matrix(XpopT - XsamT, nrow=1) / (N - n)
predictions <- predict(
  sim, X=list(reg1=XpopR), weights = N-n,
  fun=function(x, p) (samplesum + x)/N,
  show.progress=FALSE
)
summary(predictions)

## ----fig.width=4, fig.height=4, fig.align="center"------------------------------------------------
sampler <- create_sampler(model, data=apisrs,
                          linpred=list(reg1=matrix(XpopT/N, nrow=1)),
                          compute.weights=TRUE)
sim <- MCMCsim(sampler, verbose=FALSE)
plot(weights(cal)/N, weights(sim)); abline(0, 1)
sum(weights(sim) * apisrs$api00)
print(summary(sim, "linpred_"), digits=6)

## ----fig.width=4, fig.height=4, fig.align="center", message=FALSE---------------------------------
sampler <- create_sampler(
  model,
  family = f_gaussian(var.model = ~vfac(prior=pr_invchisq(df="modeled"))),
  linpred=list(reg1=matrix(XpopR, nrow=1)),
  data=apisrs, compute.weights=TRUE
)
sim <- MCMCsim(sampler, burnin=1000, n.iter=5000, thin=2, verbose=FALSE)
(summ <- summary(sim))
plot(sim, "vfac1_df")
acceptance_rates(sim)
compute_DIC(sim)
predictions <- predict(sim, newdata=apipop[-m, ], show.progress=FALSE,
                       fun=function(x, p) (samplesum + sum(x))/N)
summary(predictions)
plot(weights(cal)/N, weights(sim)); abline(0, 1)
summary(get_means(sim, "Q_")[["Q_"]])

