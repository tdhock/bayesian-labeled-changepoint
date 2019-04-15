source("packages.R")

pred <- seq(-2, 4, by=0.01)
thresh <- 1
mynorm <- function(x)pnorm(x, sd=0.1)
pred.dt <- rbind(
  data.table(k=0, pred, prob=mynorm(-pred)),
  data.table(k=1, pred, prob=mynorm(thresh-pred)-mynorm(-pred)),
  data.table(k=2, pred, prob=1-mynorm(thresh-pred)))
pred.dt[, list(
  sum.prob=sum(prob)
), by=list(pred)][order(sum.prob)]
pred.dt[, K := factor(k)]
ggplot()+
  geom_line(aes(
    pred, prob, color=K, group=K),
    data=pred.dt)

max.dt <- pred.dt[, .SD[max(prob)==prob], by=list(pred)]
max.dt[, list(n.max=.N), by=list(pred)][order(n.max)]
max.dt[k==1]

ggplot()+
  geom_point(aes(
    pred, k),
    data=max.dt)

pred <- seq(-2, 4, by=0.01)
thresh <- 1
mynorm <- function(x)1/(1+exp(-10*x))
pred.dt <- rbind(
  data.table(k=0, pred, prob=mynorm(-pred)),
  data.table(k=1, pred, prob=mynorm(thresh-pred)-mynorm(-pred)),
  data.table(k=2, pred, prob=1-mynorm(thresh-pred)))
pred.dt[, list(
  sum.prob=sum(prob)
), by=list(pred)][order(sum.prob)]
pred.dt[, K := factor(k)]
ggplot()+
  geom_line(aes(
    pred, prob, color=K, group=K),
    data=pred.dt)

max.dt <- pred.dt[, .SD[max(prob)==prob], by=list(pred)]
max.dt[, list(n.max=.N), by=list(pred)][order(n.max)]
max.dt[k==1]

ggplot()+
  geom_point(aes(
    pred, k),
    data=max.dt)

## plot both.
pred <- seq(-2, 4, by=0.01)
thresh <- 1
fun.list <- list(
  sigmoid=function(x, s)1/(1+exp(-x/s)),
  pnorm=function(x, s)pnorm(x, sd=s))
pred.dt <- data.table(link.fun=names(fun.list))[, {
  mynorm <- fun.list[[link.fun]]
  data.table(scale=c(1, 0.1))[, rbind(
    data.table(k=0, pred, prob=mynorm(-pred, scale)),
    data.table(k=1, pred, prob=mynorm(thresh-pred, scale)-mynorm(-pred, scale)),
    data.table(k=2, pred, prob=1-mynorm(thresh-pred, scale))
  ), by=list(scale)]
}, by=list(link.fun)]
pred.dt[, K := factor(k)]
max.dt <- pred.dt[, .SD[which.max(prob)], by=list(pred, link.fun, scale)]
gg <- ggplot()+
  ggtitle("ordinal regression does not always predict the right class")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(link.fun ~ scale, labeller=label_both)+
  geom_line(aes(
    pred, prob, color=K, group=K),
    data=pred.dt)+
  geom_point(aes(
    pred, -0.1, color=K),
    data=max.dt)+
  scale_color_manual(values=c("blue", "red", "black"))
png("figure-ordinal.png", 6, 6, units="in", res=100)
print(gg)
dev.off()
