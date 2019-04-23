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

some.error <- penaltyLearning::labelError(
  selection.dt, labels.possible, model.changes,
  change.var="change.pos",
  label.vars=c("min", "max"),
  problem.vars=c("profile.id", "chromosome"))
pred.selection <- data.table(some.error$label.errors)

sel.changes <- model.changes[selection.dt, on=list(
  profile.id, chromosome, n.segments), nomatch=0L]
sel.changes[, change0 := change.pos]
setkey(sel.changes, profile.id, chromosome, change.pos, change0)
setkey(labels.possible, profile.id, chromosome, min, max)
over.dt <- foverlaps(sel.changes, labels.possible, nomatch=0L)

mynorm <- function(x, s)1/(1+exp(-s*x))
range.vec <- selection.dt[, {
  l <- c(min.log.lambda, max.log.lambda)
  range(l[is.finite(l)])
}]
pred <- seq(range.vec[1]-10, range.vec[2]+10, l=201)
pred.dt <- data.table(scale=c(0.5, 1, 2, 5, 50))[, {
  data.table(pred)[, {
    prob.dt <- data.table(over.dt)
    prob.dt[, prob := mynorm(
      max.log.lambda-pred, scale)-mynorm(min.log.lambda-pred, scale)]
    my.dt <- prob.dt[, {
      ##print(.SD)
      ##browser()
      list(
        total.prob=sum(prob)
    )}, by=list(
      profile.id, chromosome, change.pos,
      min, max, annotation, possible.changes)]
    my.dt
  }, by=list(pred)]
}, by=list(scale)]

## Now we have the predicted probability of every changepoint variable.
pred.dt[, range(total.prob)]
pred.dt[, hist(total.prob)]

## TODO compute mean/max prob over all changepoint variables in each
## labeled region.
stats.wide <- pred.dt[, list(
  mean.prob=mean(total.prob),
  max.prob=max(total.prob)
), by=list(
  scale, pred, profile.id, chromosome,
  min, max, annotation, possible.changes)]
stats.tall <- melt(
  stats.wide,
  measure.vars=c("mean.prob", "max.prob"))
ggplot()+
  ggtitle("Predicted changes in each label")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome + min + max + annotation ~ scale, scales="free")+
  geom_point(aes(
    pred, value, color=variable),
    shape=1,
    data=stats.tall)

stats.tall[, prob0 := dbinom(0, possible.changes, value) ]
stats.tall[, lik := ifelse(annotation=="normal", prob0, 1-prob0)]
stats.tall[, logTRUE := ifelse(
  annotation=="breakpoint",
  pbinom(0, possible.changes, value, lower.tail=FALSE, log.p=TRUE),
  dbinom(0, possible.changes, value, log=TRUE))]
stats.tall[, plot(log(lik), logTRUE)]

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome + min + max + annotation ~ scale, scales="free")+
  geom_point(aes(
    pred, logTRUE, color=variable),
    shape=1,
    data=stats.tall)

ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome + min + max + annotation ~ scale, scales="free")+
  geom_line(aes(
    pred, logTRUE, color=variable),
    shape=1,
    data=stats.tall)

##plot 01-bounded neg log like with label error.
pred.selection[, errors := fp+fn]
norm <- function(x)(x-min(x))/(max(x)-min(x))
stats.tall[, loss := norm(-logTRUE)]
gg <- ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(scale ~ chromosome + annotation, labeller=label_both)+
  coord_cartesian(ylim=c(0, 5))+
  geom_line(aes(
    pred, -logTRUE, color=variable, size=variable),
    data=stats.tall)+
  geom_segment(aes(
    min.log.lambda, errors,
    xend=max.log.lambda, yend=errors,
    color=variable, size=variable),
    data=data.table(variable="incorrect labels", pred.selection))+
  scale_size_manual(values=c(
    "incorrect labels"=0.75, "mean.prob"=1.5, max.prob=1.2))+
  scale_color_manual(values=c(
    "incorrect labels"="black", "mean.prob"="red", max.prob="deepskyblue"))
png("figure-ordinal-loss-max.png", 10, 4, units="in", res=100)
print(gg)
dev.off()
