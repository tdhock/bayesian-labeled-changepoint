works_with_R(
  "3.5.1",
  data.table="1.11.8",
  ##icenReg="2.0.9",
  coda="0.19.2",
  MLEcens="0.1.4",
  DPpackage="1.1.7.4",
  flexsurv="1.1.1",
  penaltyLearning="2018.9.4")
library(survival)
requireGitHub::requireGitHub_package(
  "pistacliffcho",
  "icenReg_devel/Code/icenReg",
  "d43ca5c15b38d333cc720fbaeaf3ccaf89898b34",
  "icenReg")
data(neuroblastomaProcessed, package="penaltyLearning")

X.mat <- neuroblastomaProcessed$feature.mat[, c("log.n", "log.hall")]
y.mat <- neuroblastomaProcessed$target.mat
match.mat <- namedCapture::str_match_variable(
  rownames(X.mat),
  profile="[0-9]+",
  "[.]",
  chrom="[0-9]+")
table(match.mat[, "chrom"])
is.test <- match.mat[, "chrom"] == "1"
is.train <- !is.test
full.dt <- data.table(y.mat, X.mat)
train.dt <- full.dt[is.train]
test.dt <- full.dt[is.test]

## penaltyLearning squared hinge loss.
fit.penaltyLearning <- penaltyLearning::IntervalRegressionUnregularized(
  train.dt[, cbind(log.n, log.hall)],
  train.dt[, cbind(min.L, max.L)])

## survival normal AFT frequentist.
fit.survival <- survival::survreg(
  Surv(min.L, max.L, type="interval2") ~ log.n + log.hall,
  train.dt, dist="gaussian")

## icenReg normal AFT frequentist and bayesian models.
fit.par <- icenReg::ic_par(
  Surv(exp(min.L), exp(max.L), type="interval2") ~ log.n + log.hall,
  train.dt,
  model="aft",
  dist="lnorm")
fit.bayes <- icenReg::ic_bayes(
  Surv(exp(min.L), exp(max.L), type="interval2") ~ log.n + log.hall,
  train.dt,
  model="aft",
  dist="lnorm")
rbind(coef(fit.bayes), coef(fit.par))

## flexsurv frequentist.
fit.flex <- flexsurv::flexsurvreg(
  Surv(exp(min.L), exp(max.L), type="interval2") ~ log.n + log.hall,
  data=train.dt,
  inits=c(1,1,0,0),
  dist="lnorm")

## DPsurvint np-bayesian: error Lapack routine dgesv: system is
## exactly singular: U[1,1] = 0
fit.dp <- train.dt[, DPpackage::DPsurvint(
  cbind(
    ifelse(min.L==-Inf, -999, min.L),
    ifelse(max.L==Inf, -999, max.L)
  ) ~ log.n + log.hall,
  mcmc=list(
    nburn=20000,nsave=10000,nskip=10,
    ndisplay=100,tune=0.125),
  prior=list(
    alpha=1,beta0=rep(0,2),Sbeta0=diag(1000,2),
    m0=0,s0=1,tau1=0.01,tau2=0.01),
  status=TRUE)]

bayes.mat <- log(icenReg::ic_sample(fit.bayes, test.dt, samples=500))
bayes.vec <- apply(bayes.mat, 1, sd)

par.mat <- log(icenReg::ic_sample(fit.par, test.dt, samples=500))
par.vec <- apply(par.mat, 1, sd)

## predictions have constant variance, as expected.
rbind(bayes.vec, par.vec)

library(ggplot2)
ggplot()+
  geom_segment(aes(
    log.n, min.L,
    xend=log.n, yend=max.L),
    data=test.dt)

rbind(
  icenReg.bayes=log(imputeCens(fit.bayes, test.dt, imputeType="median")),
  icenReg.par=log(imputeCens(fit.par, test.dt, imputeType="median")),
  sample.mean.par=rowMeans(par.mat),
  survival=predict(fit.survival, test.dt),
  penaltyLearning=as.numeric(predict(fit.penaltyLearning, test.dt[, cbind(log.n, log.hall)])))[, 1:5]
