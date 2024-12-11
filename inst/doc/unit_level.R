## ----echo=FALSE-----------------------------------------------------------------------------------
options(width=100)  # width of output

## ----message=FALSE--------------------------------------------------------------------------------
library(survey)
data(api)
apipop$cname <- as.factor(apipop$cname)
apisrs$cname <- factor(apisrs$cname, levels=levels(apipop$cname))

## ----message=FALSE--------------------------------------------------------------------------------
library(mcmcsae)
mod <- api00 ~ 
  reg(~ ell + meals + stype + hsg + col.grad + grad.sch, name="beta") +
  gen(factor = ~ iid(cname), name="v")
sampler <- create_sampler(mod, data=apisrs)
sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)
(summary(sim))

## ----fig.width=7, fig.height=7, fig.align="center"------------------------------------------------
N <- table(apipop$cname)
samplesums <- tapply(apisrs$api00, apisrs$cname, sum)
samplesums[is.na(samplesums)] <- 0  # substitute 0 for out-of-sample areas
m <- match(apisrs$cds, apipop$cds)  # population units in the sample
res <- predict(sim, newdata=apipop, labels=names(N),
         fun=function(x, p) (samplesums + tapply(x[-m], apipop$cname[-m], sum ))/N,
         show.progress=FALSE)
(summ <- summary(res))
theta <- c(tapply(apipop$api00, apipop$cname, mean))  # true population quantities
plot_coef(summ, list(est=theta), n.se=2, est.names=c("mcmcsae", "true"), maxrows=30)

## -------------------------------------------------------------------------------------------------
apisrs$target.met <- as.numeric(apisrs$sch.wide == "Yes")
sampler <- create_sampler(update(mod, target.met ~ .), family="binomial", data=apisrs)
sim <- MCMCsim(sampler, store.all=TRUE, verbose=FALSE)
summary(sim)

## ----fig.width=7, fig.height=7, fig.align="center"------------------------------------------------
samplesums <- tapply(apisrs$target.met, apisrs$cname, sum)
samplesums[is.na(samplesums)] <- 0  # substitute 0 for out-of-sample areas
res <- predict(sim, newdata=apipop, labels=names(N),
         fun=function(x, p) (samplesums + tapply(x[-m], apipop$cname[-m], sum ))/N,
         show.progress=FALSE)
(summ <- summary(res))
theta <- c(tapply(apipop$sch.wide == "Yes", apipop$cname, mean))  # true population quantities
plot_coef(summ, list(est=theta), n.se=2, est.names=c("mcmcsae", "true"), maxrows=30)

