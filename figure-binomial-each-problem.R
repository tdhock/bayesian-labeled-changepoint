source("packages.R")

data(neuroblastomaDetailed, package="bams")
data(neuroblastoma, package="neuroblastoma")

all.labels <- data.table(neuroblastomaDetailed)
(label.counts <- dcast(all.labels, profile.id + chromosome ~ annotation, length))
some.probs <- label.counts[`>0breakpoints`>0 & `1breakpoint`>0 & normal>0][5:6]

some.profiles <- data.table(neuroblastoma$profiles)[some.probs, on=list(profile.id, chromosome)]
some.profiles[, pos0 := position]
possible.changes <- some.profiles[, list(
  change.pos=diff(position)/2+position[-.N]
), by=list(profile.id, chromosome)]
possible.changes[, pos0 := change.pos]
setkey(possible.changes, profile.id, chromosome, change.pos, pos0)

(some.ann <- all.labels[some.probs, on=list(profile.id, chromosome)])
setkey(some.ann, profile.id, chromosome, min, max)
over2 <- foverlaps(some.ann, possible.changes, nomatch=0L)

label2 <- over2[, list(
  possible.changes=.N
), by=list(profile.id, chromosome, min, max, annotation)]
setkey(label2, profile.id, chromosome, min, max)

ggplot()+
  ggtitle("possible changes in each label")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome ~ ., scales="free")+
  geom_tallrect(aes(
    xmin=min, xmax=max, fill=annotation),
    alpha=0.5,
    color="grey",
    data=label2)+
  geom_text(aes(
  (min+max)/2, Inf, label=possible.changes),
  vjust=1.5,
    data=label2)+
  scale_fill_manual(values=change.colors)+
  geom_point(aes(
    position, logratio),
    data=some.profiles)

prob2 <- label2[, {
  prob <- seq(0, 1, l=1001)
  prob0 <- dbinom(0, possible.changes, prob)
  lik <- if(annotation=="normal"){
    prob0
  }else if(annotation==">0breakpoints"){
    1-prob0
  }else if(annotation=="1breakpoint"){
    dbinom(1, possible.changes, prob)
  }else{
    stop("likelihood undefined for annotation ", annotation)
  }
  data.table(prob, lik)
}, by=list(profile.id, chromosome, min, max, annotation)]
prob.max <- prob2[, {
  .SD[lik==max(lik)]
}, by=list(profile.id, chromosome, min, max, annotation)]
gg <- ggplot()+
  ggtitle("Likelihood function to maximize for each label")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome + min + max ~ ., scales="free")+
  geom_line(aes(
    prob, lik),
    data=prob2)+
  geom_point(aes(
    prob, lik),
    data=prob.max)+
  geom_text(aes(
    0.1, 0.5, label=annotation),
    data=label2)

prob.dots <- label2[, {
  pred.changes <- 0:possible.changes
  prob <- pred.changes/possible.changes
  prob0 <- dbinom(0, possible.changes, prob)
  lik <- if(annotation=="normal"){
    prob0
  }else if(annotation==">0breakpoints"){
    1-prob0
  }else if(annotation=="1breakpoint"){
    dbinom(1, possible.changes, prob)
  }else{
    stop("likelihood undefined for annotation ", annotation)
  }
  data.table(prob, lik, pred.changes)
}, by=list(profile.id, chromosome, min, max, annotation)]
prob.dots[, prediction := ifelse(lik==max(lik), "optimal", "sub-optimal"), by=list(
  profile.id, chromosome, min, max, annotation)]
gg <- ggplot()+
  ggtitle("Likelihood function to maximize for each label")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(. ~ profile.id + chromosome + annotation + min + max, scales="free")+
  geom_line(aes(
    prob, lik),
    data=prob2)+
  geom_point(aes(
    prob, lik, fill=prediction),
    shape=21,
    data=prob.dots)+
  scale_x_continuous(breaks=c(0.5, 1))+
  scale_fill_manual(values=c(optimal="black", "sub-optimal"="white"))
png("figure-binomial-each-problem-prob-dots.png", 13, 4, units="in", res=100)
print(gg)
dev.off()

## now compute total log likelihood by summing over labels for each
## problem.
some.loss <- some.profiles[, {
  max.segments <- 150
  fit <- Segmentor3IsBack::Segmentor(logratio, model=2, Kmax=max.segments)
  change.vec <- diff(position)/2+position[-.N]
  start.pos <- c(position[1], change.vec)
  end.pos <- c(change.vec, position[.N])
  seg.list <- list()
  for(n.segments in 1:max.segments){
    end.i <- fit@breaks[n.segments, 1:n.segments]
    start.i <- c(1, end.i[-length(end.i)]+1)
    seg.list[[n.segments]] <- data.table(
      n.segments,
      loss=fit@likelihood[[n.segments]],
      segStart=list(start.pos[start.i]),
      segEnd=list(end.pos[end.i]))
  }
  do.call(rbind, seg.list)
}, by=list(profile.id, chromosome)]
some.segments <- some.loss[, data.table(
  segStart=segStart[[1]],
  segEnd=segEnd[[1]]
), by=list(profile.id, chromosome, n.segments)]
some.selection <- some.loss[, {
  penaltyLearning::modelSelection(.SD[, .(loss, n.segments)], "loss", "n.segments")
}, by=list(profile.id, chromosome)]
model.changes <- some.segments[, list(
  change.pos=segStart[-1]
), by=list(profile.id, chromosome, n.segments)]
some.error <- penaltyLearning::labelError(
  some.selection, label2, model.changes,
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

## now plot probability as a function of penalty.
prob.pen <- prob.dots[pred.selection, on=list(
  profile.id, chromosome, min, max, annotation, pred.changes)]
pred.selection[, errors := fp+fn]
pred.selection[, res.max := pred.changes-max.changes]
pred.selection[, res.min := min.changes-pred.changes]
pred.selection[, res.changes := ifelse(
  0<res.max, res.max, ifelse(
    0<res.min, res.min, 0))]
pred.selection[, res.thresh := ifelse(res.changes<5, res.changes, 5)]
seg.size <- 1
prob.pen[, neg.lik := -lik+1]
gg <- ggplot()+
  ggtitle("Probability of each label per penalty")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(metric ~ profile.id + chromosome + annotation + min + max, scales="free")+
  coord_cartesian(ylim=c(-0.5, 3.5))+
  geom_segment(aes(
    min.log.lambda, neg.lik,
    xend=max.log.lambda, yend=neg.lik),
    size=seg.size,
    data=data.table(metric="-likelihood", prob.pen))+
  geom_segment(aes(
    min.log.lambda, errors,
    xend=max.log.lambda, yend=errors),
    size=seg.size,
    data=data.table(metric="label errors", pred.selection))+
  geom_segment(aes(
    min.log.lambda, res.changes,
    xend=max.log.lambda, yend=res.changes),
    size=seg.size,
    data=data.table(metric="label residual", pred.selection))+
  ylab("")
png("figure-binomial-each-problem-label.png", 12, 4.5, units="in", res=100)
print(gg)
dev.off()

## plot total log likelihood.
label.lik <- pred.selection[, lik := {
  prob <- pred.changes/possible.changes
  prob0 <- dbinom(0, possible.changes, prob)
  ifelse(
    annotation=="normal", prob0, ifelse(
      annotation==">0breakpoints", 1-prob0,
      dbinom(1, possible.changes, prob)))
}]
nan <- label.lik[!is.finite(lik)]
if(nrow(nan)){
  print(nan)
  stop("some NaN in likelihood, this usually indicates counts>possible")
}
label.lik[, log.lik := log(lik+1)]
mean.dt <- label.lik[, list(
  mean.log.lik=mean(log.lik)
), by=list(profile.id, chromosome, n.segments, min.log.lambda, max.log.lambda)]
res.dt <- pred.selection[, list(
  total.res=sum(res.changes)
), by=list(profile.id, chromosome, n.segments, min.log.lambda, max.log.lambda)]
res.dt[, res.thresh := ifelse(total.res<10, total.res, 10)]
mean.dt[, errors := -mean.log.lik]
interval.dt <- penaltyLearning::targetIntervals(
  mean.dt, c("profile.id", "chromosome"))
gg <- ggplot()+
  ggtitle("Label error and log likelihood per penalty")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(metric ~ profile.id + chromosome, scales="free", labeller=label_both)+
  xlab("")+
  ylab("")+
  penaltyLearning::geom_tallrect(aes(
    xmin=min.log.lambda, xmax=max.log.lambda),
    color=NA,
    fill="grey",
    data=interval.dt)+
  geom_segment(aes(
    min.log.lambda, -mean.log.lik,
    xend=max.log.lambda, yend=-mean.log.lik),
    size=seg.size,
    data=data.table(metric="
-log(likelihood)", mean.dt))+
  geom_segment(aes(
    min.log.lambda, res.thresh,
    xend=max.log.lambda, yend=res.thresh),
    size=seg.size,
    data=data.table(metric="
label residual", res.dt))+
  geom_segment(aes(
    min.log.lambda, errors,
    xend=max.log.lambda, yend=errors),
    size=seg.size,
    data=data.table(metric="
label errors", some.error$model.errors))
png("figure-binomial-each-problem.png", 6, 4, units="in", res=100)
print(gg)
dev.off()
