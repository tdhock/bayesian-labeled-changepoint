source("packages.R")
data(neuroblastoma, package="neuroblastoma")

pid <- 2
one.pro <- data.table(neuroblastoma$profiles)[profile.id==pid]
one.pro[, pos0 := position]
(one.ann <- data.table(neuroblastoma$annotations)[profile.id==pid])
setkey(one.pro, profile.id, chromosome, position, pos0)
setkey(one.ann, profile.id, chromosome, min, max)
over.dt <- foverlaps(one.ann, one.pro, nomatch=0L)
label.dt <- over.dt[, list(
  possible.changes=.N-1
), by=list(profile.id, chromosome, min, max, annotation)]
labeled.pro <- one.pro[label.dt, on=list(profile.id, chromosome)]
change.dt <- one.pro[, list(change.pos=diff(position)/2+position[-.N]), by=list(profile.id, chromosome)]
change.dt[, pos0 := change.pos]
setkey(change.dt, profile.id, chromosome, change.pos, pos0)
over2 <- foverlaps(one.ann, change.dt, nomatch=0L)
label2 <- over2[, list(
  possible.changes=.N
), by=list(profile.id, chromosome, min, max, annotation)]

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
    data=labeled.pro)

prob.dt <- label2[, {
  changes <- 0:possible.changes
  prob0 <- dbinom(changes, possible.changes, 0)
  data.table(changes, prob=if(annotation=="normal")prob0 else 1-prob0)
}, by=list(profile.id, chromosome, min, max, annotation)]

ggplot()+
  ggtitle("possible changes in each label")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome ~ ., scales="free")+
  geom_point(aes(
    changes, prob),
    data=prob.dt)

prob2 <- label2[, {
  prob <- seq(0, 1, l=1001)
  prob0 <- dbinom(0, possible.changes, prob)
  lik <- if(annotation=="normal")prob0 else 1-prob0
  data.table(prob, lik)
}, by=list(profile.id, chromosome, min, max, annotation)]

gg <- ggplot()+
  ggtitle("Likelihood function to maximize")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome ~ ., scales="free")+
  geom_line(aes(
    prob, lik),
    data=prob2)+
  geom_text(aes(
    0.1, 0.5, label=annotation),
    data=label2)
png("figure-binomial-each-label-prob.png", 6, 6, units="in", res=100)
print(gg)
dev.off()

sigmoid <- function(x)1/(1+exp(-x))
inv.sigmoid <- function(p)-log(1/p-1)
gg <- ggplot()+
  ggtitle("Likelihood function to maximize")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(chromosome ~ ., scales="free")+
  geom_line(aes(
    inv.sigmoid(prob), lik),
    data=prob2[is.finite(inv.sigmoid(prob))])+
  geom_text(aes(
    inv.sigmoid(0.1), 0.5, label=annotation),
    data=label2)
png("figure-binomial-each-label.png", 6, 6, units="in", res=100)
print(gg)
dev.off()
