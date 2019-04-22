source("packages.R")

data(neuroblastoma, package="neuroblastoma")
pid <- 2
labels.dt <- data.table(neuroblastoma$annotations)[profile.id==pid]
profiles.dt <- data.table(
  neuroblastoma$profiles)[labels.dt, on=list(profile.id, chromosome)]
all.changes <- profiles.dt[, list(
  change.pos=diff(position)/2+position[-.N]
), by=list(profile.id, chromosome, min, max, annotation)]
labels.possible <- all.changes[min < change.pos & change.pos < max, list(
  possible.changes=.N), by=list(profile.id, chromosome, min, max, annotation)]

## now compute total log likelihood by summing over labels for each
## problem.
loss.dt <- profiles.dt[, {
  max.segments <- .N
  fit <- jointseg::Fpsn(logratio, max.segments)
  change.vec <- diff(position)/2+position[-.N]
  start.pos <- c(position[1], change.vec)
  end.pos <- c(change.vec, position[.N])
  seg.list <- list()
  for(n.segments in 1:max.segments){
    end.i <- fit$t.est[n.segments, 1:n.segments]
    if(n.segments==max.segments){
      exp.dt <- data.table(expected=1:.N, end.i)[expected != end.i]
      if(nrow(exp.dt)){
        print(exp.dt)
        stop("max seg model does not change after every data point")
      }
    }
    start.i <- c(1, end.i[-length(end.i)]+1)
    seg.list[[n.segments]] <- data.table(
      n.segments,
      loss=fit$cost[[n.segments]],
      segStart=list(start.pos[start.i]),
      segEnd=list(end.pos[end.i]))
  }
  do.call(rbind, seg.list)
}, by=list(profile.id, chromosome)]

segments.dt <- loss.dt[, data.table(
  segStart=segStart[[1]],
  segEnd=segEnd[[1]]
), by=list(profile.id, chromosome, n.segments)]
selection.dt <- loss.dt[, {
  penaltyLearning::modelSelection(
    .SD[, .(loss, n.segments)], "loss", "n.segments")
}, by=list(profile.id, chromosome)]
model.changes <- segments.dt[, list(
  change.pos=segStart[-1]
), by=list(profile.id, chromosome, n.segments)]

mynorm <- function(x, s)1/(1+exp(-s*x))
range.vec <- selection.dt[, {
  l <- c(min.log.lambda, max.log.lambda)
  range(l[is.finite(l)])
}]
pred <- seq(range.vec[1]-10, range.vec[2]+10, l=1001)
pred.dt <- data.table(scale=c(0.5, 1, 2, 5, 50))[, {
  data.table(pred)[, {
    data.table(selection.dt)[, {
      prob <- mynorm(max.log.lambda-pred, scale)-mynorm(min.log.lambda-pred, scale)
      prob.dt <- data.table(.SD, prob)[order(-prob)]
      prob.dt[, cum.prob := cumsum(prob)]
      result <- prob.dt[cum.prob < 0.99]
      if(nrow(result)==0){
        prob.dt[1]
      }else{
        result
      }
    }]
  }, by=list(pred)]
}, by=list(scale)]

sel.changes <- model.changes[pred.dt, on=list(profile.id, chromosome, n.segments), allow.cartesian=TRUE, nomatch=0L]

## Now we have the predicted probability of every changepoint variable.
mean.changes <- sel.changes[, {
  list(
    prob=mean(prob))
}, by=list(profile.id, chromosome, pred, scale, change.pos)]

## TODO compute mean/max prob over all changepoint variables in each labeled region.

some.error <- penaltyLearning::labelError(
  selection.dt, labels.possible, model.changes,
  change.var="change.pos",
  label.vars=c("min", "max"),
  problem.vars=c("profile.id", "chromosome"))
pred.selection <- data.table(some.error$label.errors)
ggplot()+
  ggtitle("Predicted changes in each label")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome + min + max + annotation ~ ., scales="free")+
  geom_segment(aes(
    min.log.lambda, pred.changes,
    xend=max.log.lambda, yend=pred.changes),
    data=pred.selection)

## TODO figure out why max prop is not 1 for these two:
pred.selection[, list(
  max.prop=max(pred.changes/possible.changes)
), by=list(profile.id, chromosome, min, max)]


pred.dt[, round.prob := round(mean.prob, 6)]
pred.dt[, prob0 := dbinom(0, possible.changes, round.prob) ]
pred.dt[!is.finite(prob0)]
pred.dt[, lik := ifelse(annotation=="normal", prob0, 1-prob0)]
pred.dt[, logTRUE := ifelse(
  annotation=="breakpoint",
  pbinom(0, possible.changes, round.prob, lower.tail=FALSE, log.p=TRUE),
  dbinom(0, possible.changes, round.prob, log=TRUE))]
pred.dt[, plot(log(lik), logTRUE)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(scale ~ chromosome + annotation)+
  geom_line(aes(
    pred, mean.prob),
    data=pred.dt)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(scale ~ chromosome + annotation)+
  geom_line(aes(
    pred, lik),
    data=pred.dt)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(scale ~ chromosome + annotation)+
  geom_line(aes(
    pred, logTRUE),
    data=pred.dt)+
  ylim(-10, 0)

gprob <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(scale ~ chromosome + annotation, labeller=label_both)+
  scale_color_discrete(breaks=c("prob", "logLik"))+
  geom_line(aes(
    pred, logTRUE, color=curve),
    data=data.table(curve="logLik", pred.dt))+
  geom_line(aes(
    pred, mean.prob, color=curve),
    data=data.table(curve="prob", pred.dt))+
  ylim(-3, 1)
png("figure-ordinal-loss-prob.png", 10, 4, units="in", res=100)
print(gprob)
dev.off()

##plot 01-bounded neg log like with label error.
pred.selection[, errors := fp+fn]
norm <- function(x)(x-min(x))/(max(x)-min(x))
pred.dt[, loss := norm(-log(lik+1))]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(scale ~ chromosome + annotation, labeller=label_both)+
  geom_line(aes(
    pred, loss, color=metric, size=metric),
    data=data.table(metric="-log.lik", pred.dt)[is.finite(loss)])+
  geom_segment(aes(
    min.log.lambda, errors,
    xend=max.log.lambda, yend=errors,
    color=metric, size=metric),
    data=data.table(metric="incorrect labels", pred.selection))+
  scale_size_manual(values=c(
    "incorrect labels"=0.75, "-log.lik"=1.5))
png("figure-ordinal-loss-bounded.png", 10, 4, units="in", res=100)
print(gg)
dev.off()

##plot usual neg log lik with label error.
pred.selection[, errors := fp+fn]
pred.dt[, loss := -logTRUE]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(scale ~ chromosome + annotation, labeller=label_both)+
  geom_line(aes(
    pred, loss, color=metric, size=metric),
    data=data.table(metric="-log.lik", pred.dt)[is.finite(loss)])+
  geom_segment(aes(
    min.log.lambda, errors,
    xend=max.log.lambda, yend=errors,
    color=metric, size=metric),
    data=data.table(metric="incorrect labels", pred.selection))+
  scale_size_manual(values=c(
    "incorrect labels"=0.75, "-log.lik"=1.5))+
  scale_y_continuous(limits=c(0, 10))
png("figure-ordinal-loss.png", 10, 4, units="in", res=100)
print(gg)
dev.off()
