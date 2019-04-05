works_with_R(
  "3.5.1",
  data.table="1.11.8",
  icenReg="2.0.9",
  flexsurv="1.0",
  penaltyLearning="2018.9.4")
library(survival)

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
  dist="lnorm")

rbind(
  icenReg.bayes=log(imputeCens(fit.bayes, test.dt, imputeType="median")),
  icenReg.par=log(imputeCens(fit.par, test.dt, imputeType="median")),
  survival=predict(fit.survival, test.dt),
  penaltyLearning=as.numeric(predict(fit.penaltyLearning, test.dt[, cbind(log.n, log.hall)])))


