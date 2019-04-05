library(flexsurv)
library(penaltyLearning)
library(survival)
data(neuroblastomaProcessed, package="penaltyLearning")
X.mat <- neuroblastomaProcessed$feature.mat[, c("log.n", "log.hall")]
y.mat <- neuroblastomaProcessed$target.mat
train.df <- data.frame(X.mat, y.mat)
fit.survival <- survival::survreg(
  Surv(min.L, max.L, type="interval2") ~ log.n + log.hall,
  train.df, dist="gaussian")
fit.survival
fit.flex <- flexsurv::flexsurvreg(
  Surv(exp(min.L), exp(max.L), type="interval2") ~ log.n + log.hall,
  data=train.df,
  dist="lnorm")

