library(DPpackage)
data(neuroblastomaProcessed, package="penaltyLearning")
x1 <- neuroblastomaProcessed$feature.mat[, "log.n"]
x2 <- neuroblastomaProcessed$feature.mat[, "log.hall"]
y <- with(neuroblastomaProcessed, cbind(
  ifelse(target.mat[,1]==-Inf, -999, target.mat[,1]),
  ifelse(target.mat[,2]==Inf, -999, target.mat[,2])))
ey <- with(neuroblastomaProcessed, cbind(
  ifelse(target.mat[,1]==-Inf, -999, exp(target.mat[,1])),
  ifelse(target.mat[,2]==Inf, -999, exp(target.mat[,2]))))

fit.dp <- DPpackage::DPsurvint(
  y ~ x1 + x2,
  mcmc=list(
    nburn=20000,nsave=10000,nskip=10,
    ndisplay=100,tune=0.125),
  prior=list(
    alpha=1,beta0=rep(0,2),Sbeta0=diag(1000,2),
    m0=0,s0=1,tau1=0.01,tau2=0.01),
  status=TRUE)
