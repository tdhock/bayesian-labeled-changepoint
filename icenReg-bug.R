library(icenReg)
library(penaltyLearning)
data(neuroblastomaProcessed, package="penaltyLearning")

X.mat <- neuroblastomaProcessed$feature.mat[, c("log.n", "log.hall")]
y.mat <- neuroblastomaProcessed$target.mat
fit.par <- ic_par(
  Surv(exp(min.L), exp(max.L), type="interval2") ~ log.n + log.hall,
  data.frame(X.mat, y.mat),
  model="aft",
  dist="lnorm")
ic_sample(fit.par)
