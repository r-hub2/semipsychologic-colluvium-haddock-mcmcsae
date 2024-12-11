## ----echo=FALSE-----------------------------------------------------------------------------------
options(width=100)  # width of output

## -------------------------------------------------------------------------------------------------
m <- 75L  # number of areas
df <- data.frame(
  area=1:m,      # area indicator
  x=runif(m)     # covariate
)
v <- rnorm(m, sd=0.5)    # true area effects
theta <- 1 + 3*df$x + v  # quantity of interest
psi <- runif(m, 0.5, 2) / sample(1:25, m, replace=TRUE)  # given variances
df$y <- rnorm(m, theta, sqrt(psi))

## ----message=FALSE--------------------------------------------------------------------------------
library(mcmcsae)
model <- y ~ reg(~ 1 + x, name="beta") + gen(factor = ~iid(area), name="v")
sampler <- create_sampler(model, sigma.fixed=TRUE, Q0=1/psi, linpred="fitted", data=df)

## -------------------------------------------------------------------------------------------------
sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)

## -------------------------------------------------------------------------------------------------
(summ <- summary(sim))

## ----fig.show='hold'------------------------------------------------------------------------------
plot(v, summ$v[, "Mean"], xlab="true v", ylab="posterior mean"); abline(0, 1)
plot(theta, summ$linpred_[, "Mean"], xlab="true theta", ylab="estimated"); abline(0, 1)
points(theta, df$y, col=2, pch=2)

## -------------------------------------------------------------------------------------------------
compute_DIC(sim)
compute_WAIC(sim, show.progress=FALSE)

## ----fig.align="center"---------------------------------------------------------------------------
plot(df$x, residuals(sim, mean.only=TRUE), xlab="x", ylab="residual"); abline(h=0)

## -------------------------------------------------------------------------------------------------
sampler <- create_sampler(model, sigma.fixed=TRUE, Q0=1/psi,
             linpred="fitted", data=df, compute.weights=TRUE)
sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)

## ----fig.show='hold'------------------------------------------------------------------------------
plot(summ$linpred_[, "Mean"], crossprod(weights(sim), df$y),
     xlab="estimate", ylab="weighted average")
abline(0, 1)
plot(psi, diag(weights(sim)), ylab="weight")

