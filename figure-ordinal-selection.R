source("packages.R")

data(neuroblastomaProcessed, package="penaltyLearning")

err.dt <- neuroblastomaProcessed$errors[profile.id==2]
range.vec <- err.dt[, {
  l <- c(min.log.lambda, max.log.lambda)
  range(l[is.finite(l)])
}]
pred <- seq(range.vec[1]-1, range.vec[2]+1, l=1001)

mynorm <- function(x)1/(1+exp(-10*x))

pred.dt <- err.dt[, {
  data.table(
    pred,
    prob=mynorm(max.log.lambda-pred)-mynorm(min.log.lambda-pred))
}, by=list(profile.id, chromosome, n.segments)]

pred.dt[, list(
  sum.prob=sum(prob)
), by=list(pred)][order(sum.prob)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ chromosome)+
  geom_line(aes(
    pred, prob, color=n.segments, group=n.segments),
    data=pred.dt)

max.dt <- pred.dt[, .SD[which.max(prob)], by=list(pred, profile.id, chromosome)]
max.stats <- max.dt[, list(
  mean.pred=mean(pred),
  max.pred=max(pred)
  ), by=list(profile.id, chromosome, n.segments)]
ffac <- function(x)factor(x, c("segments", "prob"))
gg <- ggplot()+
  theme_bw()+
  theme(panel.grid=element_blank(), panel.margin=grid::unit(0, "lines"))+
  facet_grid(facet ~ chromosome, scales="free")+
  scale_color_gradient(low="grey", high="red")+
  geom_segment(aes(
    min.log.lambda, n.segments,
    xend=max.log.lambda, yend=n.segments),
    size=1,
    data=data.table(facet=ffac("segments"), err.dt))+
  geom_vline(aes(
    xintercept=min.log.lambda),
    color="grey",
    data=err.dt)+
  geom_text(aes(
    max.log.lambda, 0, label=n.segments),
    hjust=1,
    color="black",
    data=data.table(facet=ffac("segments"), err.dt))+
  geom_line(aes(
    pred, prob, color=n.segments, group=n.segments),
    data=data.table(facet=ffac("prob"), pred.dt))+
  geom_point(aes(
    pred, 1.1, color=n.segments),
    data=data.table(facet=ffac("prob"), max.dt))+
  geom_text(aes(
    max.pred, 1.2, label=n.segments, color=n.segments),
    hjust=1,
    data=data.table(facet=ffac("prob"), max.stats))+
  ylab("")
png("figure-ordinal-selection.png", 40, 6, units="in", res=100)
print(gg)
dev.off()
