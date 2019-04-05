source("packages.R")
data(neuroblastomaProcessed, package="penaltyLearning")

X.mat <- neuroblastomaProcessed$feature.mat[, "log2.n", drop=FALSE]
y.mat <- neuroblastomaProcessed$target.mat
match.dt <- data.table(namedCapture::str_match_variable(
  rownames(X.mat),
  profile="[0-9]+",
  "[.]",
  chromosome="[0-9]+"))
is.test <- match.dt$chromosome == "1"
is.train <- !is.test
full.dt <- match.dt[, data.table(
    profile.id=profile,
    chromosome,
    y.mat, X.mat, set=ifelse(is.train, "train", "test"))]
train.dt <- full.dt[is.train]

## survival normal AFT frequentist.
fit.survival <- survival::survreg(
  Surv(min.L, max.L, type="interval2") ~ log2.n,
  train.dt, dist="gaussian")

targets.tall <- melt(
  full.dt,
  measure.vars=c("min.L", "max.L"),
  value.name="log.penalty")[is.finite(log.penalty)]
targets.tall[, limit := sub(".L", "", variable)]
full.dt[, pred.mean := predict(fit.survival, full.dt)]
q.val <- 0.025 # 95% CI.
full.dt[, pred.lo := qnorm(q.val, pred.mean, fit.survival$scale, lower.tail=TRUE)]
full.dt[, pred.hi := qnorm(q.val, pred.mean, fit.survival$scale, lower.tail=FALSE)]
model.color <- "deepskyblue"
model.alpha <- 0.25
gg <- ggplot()+
  ggtitle("Linear model with 95% confidence interval")+
  geom_line(aes(
    log2.n, pred.mean),
    color=model.color,
    size=1,
    data=full.dt)+
  geom_ribbon(aes(
    log2.n, ymin=pred.lo, ymax=pred.hi),
    alpha=model.alpha,
    fill=model.color,
    data=full.dt)+
  geom_point(aes(
    log2.n, log.penalty, fill=limit, color=set),
    shape=21,
    data=targets.tall)+
  scale_color_manual(values=c(train="black", test="red"))+
  scale_fill_manual(values=c(min="white", max="black"))+
  ylab("log(penalty)")
png("figure-linear-error-bands-regression.png", 6, 6, units="in", res=100)
print(gg)
dev.off()


test.dt <- full.dt[set=="test"]
test.dt[, residual := penaltyLearning::targetIntervalResidual(
  cbind(min.L, max.L), pred.mean)]
test.dt[, min.margin := pred.mean-min.L]
test.dt[, max.margin := max.L-pred.mean]
test.dt[, name.margin := ifelse(
  residual==0, ifelse(
    min.margin < max.margin,
    "min", "max"),
  NA)]
test.dt[, margin := ifelse(name.margin=="min", min.margin, max.margin)]
test.dt[, sign.res := sign(residual)]
interesting <- rbind(test.dt[residual!=0, {
  data.table(.SD[order(-abs(residual))][c(1, .N)], distance=c("far", "close"))
}, by=list(sign.res)], test.dt[residual==0, {
  data.table(.SD[order(-margin)][c(1, .N)], distance=c("far", "close"))
}, by=list(name.margin)])
interesting[, Errors := ifelse(is.na(name.margin), 1, 0)]
interesting[, Limit := ifelse(
  is.na(name.margin), ifelse(
    sign.res<0, "lower", "upper"),
  ifelse(name.margin=="min", "lower", "upper"))]
int.err <- neuroblastomaProcessed$errors[interesting, on=list(profile.id, chromosome)]
int.tall <- melt(
  int.err,
  measure.vars=c("errors", "n.segments"))
dot.dt <- int.tall[min.log.lambda < pred.mean & pred.mean < max.log.lambda]
gg <- ggplot()+
  ggtitle("Linear model predictions with 95% confidence intervals for 8 test data sets")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(
    variable ~ Errors + Limit + distance + profile.id,
    labeller=label_both,
    scales="free_y")+
  geom_vline(aes(
    xintercept=pred.mean),
    color=model.color,
    data=interesting)+
  geom_tallrect(aes(
    xmin=pred.lo, xmax=pred.hi),
    fill=model.color,
    color=NA,
    alpha=model.alpha,
    data=interesting)+
  geom_segment(aes(
    min.log.lambda, value,
    xend=max.log.lambda, yend=value),
    data=int.tall)+
  geom_point(aes(
    pred.mean, value),
    color=model.color,
    shape=21,
    fill="white",
    data=dot.dt)+
  xlab("log(penalty)")
png("figure-linear-error-bands.png", 10, 6, units="in", res=100)
print(gg)
dev.off()
