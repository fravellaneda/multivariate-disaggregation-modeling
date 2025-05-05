#Simulations

library(sf)
library(INLA)

source("aux_fun.R")
    
# Stack for estimation
# Create the stack for l1
stk.l1 <- inla.stack(
  data = list(l = l.1),
  A = list(A.m),
  effects = list(list(
    spatial.field.sh = 1:nv,
    Intercept.l1 = rep(1, nv)
  )),
  tag = "l1")

# Create the stack for l2
stk.l2 <- inla.stack(
  data = list(l = l.2),
  A = list(A.m),
  effects = list(list(
    sh.copy.l2 = 1:nv,
    spatial.field.l2 = 1:nv,
    Intercept.l2 = rep(1, nv)
  )),
  tag = "l2")

# Create the stack for l3
stk.l3 <- inla.stack(
  data = list(l = l.3),
  A = list(A.m),
  effects = list(list(
    sh.copy.l3 = 1:nv,
    l2.copy.l3 = 1:nv,
    spatial.field.l3 = 1:nv,
    Intercept.l3 = rep(1, nv)
  )),
  tag = "l3")

A.pr <- inla.spde.make.A(mesh = mesh, loc = dp)

# Stack for predicting l1
stk.l1.pre <- inla.stack(tag = "l1.pre",
                         data = list(l = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                         A = list(A.pr),
                         effects = list(list(
                           Intercept.l1 = rep(1, nv),
                           spatial.field.sh = 1:nv
                         )))

# Stack for predicting l2
stk.l2.pre <- inla.stack(tag = "l2.pre",
                         data = list(l = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                         A = list(A.pr),
                         effects = list(list(
                           Intercept.l2 = rep(1, nv),
                           sh.copy.l2 = 1:nv,
                           spatial.field.l2 = 1:nv
                         )))

# Stack for predicting l3
stk.l3.pre <- inla.stack(tag = "l3.pre",
                         data = list(l = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                         A = list(A.pr),
                         effects = list(list(
                           Intercept.l3 = rep(1, nv),
                           sh.copy.l3 = 1:nv,
                           l2.copy.l3 = 1:nv,
                           spatial.field.l3 = 1:nv
                         )))

# Stack with shared effect
stk.shared <- inla.stack(tag = "shared",
                         data = list(l = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                         A = list(A.pr),
                         effects = list(spatial.field.sh = 1:nv))

stk.l2.spec <- inla.stack(tag = "l2.spec",
                          data = list(l = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                          A = list(A.pr),
                          effects = list(spatial.field.l2 = 1:nv))

stk.l3.spec <- inla.stack(tag = "l3.spec",
                          data = list(l = matrix(NA, nrow=nrow(dp),ncol=n.pp)),
                          A = list(A.pr),
                          effects = list(spatial.field.l3 = 1:nv))


join.stack <- inla.stack(
  stk.l1, stk.l2, stk.l3,
  stk.l1.pre, stk.l2.pre, stk.l3.pre,
  stk.shared,stk.l2.spec, stk.l3.spec)

hyper.eps <- list(hyper = list(prec = list(prior = 'pc.prec',
                                           param = c(1, 0.1))))

hyper <- list(beta = list(prior = 'normal', param = c(0, 10)))

form <- l ~ 0 + Intercept.l1 + Intercept.l2 + Intercept.l3 +
  f(spatial.field.sh, model = spde) +
  f(sh.copy.l2, copy = "spatial.field.sh", fixed = FALSE, hyper = hyper) +
  f(spatial.field.l2, model = spde) +
  f(sh.copy.l3, copy = "spatial.field.sh", fixed = FALSE, hyper = hyper) +
  f(l2.copy.l3, copy = "spatial.field.l2", fixed = FALSE, hyper = hyper) +
  f(spatial.field.l3, model = spde)

print("modeling")

t1 <-Sys.time()
res <- inla(
  formula=form, verbose = FALSE,
  data = inla.stack.data(join.stack, spde = spde),
  family = rep("gaussian", 3),
  control.family = list(hyper.eps, hyper.eps, hyper.eps),
  control.predictor = list(A = inla.stack.A(join.stack),
                           compute = TRUE, link = 1),
  control.compute = list(dic = TRUE,config=TRUE),
  #control.mode = list(theta = theta.ini.e, restart = TRUE),
  #control.inla = list(int.strategy = 'eb')
)
t2 <- Sys.time()
print("Time modeling")
print(t2-t1)

#res$mode$theta

idx.l1 <- inla.stack.index(join.stack, 'l1.pre')$data
idx.l2 <- inla.stack.index(join.stack, 'l2.pre')$data
idx.l3 <- inla.stack.index(join.stack, 'l3.pre')$data

idx.sh <- inla.stack.index(join.stack, 'shared')$data
idx.l2.spec <- inla.stack.index(join.stack, 'l2.spec')$data
idx.l3.spec <- inla.stack.index(join.stack, 'l3.spec')$data

# Prediction

dpdf$mean.l1 <- res$summary.linear.predictor[idx.l1, "mean"]
dpdf$mean.l2 <- res$summary.linear.predictor[idx.l2, "mean"]
dpdf$mean.l3 <- res$summary.linear.predictor[idx.l3, "mean"]

dpdf$l1.lower <- res$summary.linear.predictor[idx.l1, "0.025quant"]
dpdf$l2.lower <- res$summary.linear.predictor[idx.l2, "0.025quant"]
dpdf$l3.lower <- res$summary.linear.predictor[idx.l3, "0.025quant"]

dpdf$l1.upper <- res$summary.linear.predictor[idx.l1, "0.975quant"]
dpdf$l2.upper <- res$summary.linear.predictor[idx.l2, "0.975quant"]
dpdf$l3.upper <- res$summary.linear.predictor[idx.l3, "0.975quant"]

dpdf$shared <- res$summary.linear.predictor[idx.sh, "mean"]
dpdf$l2.spec <- res$summary.linear.predictor[idx.l2.spec, "mean"]
dpdf$l3.spec <- res$summary.linear.predictor[idx.l3.spec, "mean"]

p.sd <- lapply(res$internal.marginals.hyperpar[1:3],varF)

tabcrp1 <- cbind(true = alpha, res$summary.fixed[, c(1:3, 5)])
# Precision of the errors
tabcrp2 <- cbind(
  true = c(e = e.sd),
  t(sapply(p.sd, function(m)
    unlist(inla.zmarginal(m, silent = TRUE))[c(1:3, 7)])))
colnames(tabcrp2) <- colnames(tabcrp1)
# Copy parameters
tabcrp3 <- cbind(
  true = beta, res$summary.hyperpar[10:12, c(1:3, 5)])
tabcrp <- rbind(tabcrp1, tabcrp2, tabcrp3)

range.est <- c()
for (name in c("sh","l2","l3")){
  spde.est <- inla.spde2.result(inla = res,
                                name = paste0("spatial.field.",name),
                                spde = spde, do.transf = TRUE)
  temp <- try(inla.zmarginal(spde.est$marginals.range.nominal[[1]],
                             silent = TRUE))
  range.est <- rbind(range.est,unlist(temp)[c(1:3, 7)])
}
tabran1 <- cbind(true = c(ran = range),range.est)
colnames(tabran1) <- colnames(tabcrp1)

m.var.est <- c()
for (name in c("sh","l2","l3")){
  spde.est <- inla.spde2.result(inla = res,
                                name = paste0("spatial.field.",name),
                                spde = spde, do.transf = TRUE)
  temp <- try(inla.zmarginal(spde.est$marginals.variance.nominal[[1]],
                             silent = TRUE))
  m.var.est <- rbind(m.var.est,unlist(temp)[c(1:3, 7)])
}
tabran2 <- cbind(true = c(var = m.var),m.var.est)
colnames(tabran2) <- colnames(tabcrp1)

tabran <- rbind(tabran1,tabran2)

# Create the directory if it doesn't exist
if (!dir.exists(paste0("./results/sim.results.",j,"/"))) {
  dir.create(paste0("./results/sim.results.",j,"/"), recursive = TRUE)
}

# Create the directory if it doesn't exist
if (!dir.exists(paste0("./results/sim.results.",j,"/maps/"))) {
  dir.create(paste0("./results/sim.results.",j,"/maps/"), recursive = TRUE)
}

allEffects <- as.data.frame(rbind(tabcrp,tabran[c(1,4,2,5,3,6),]))
write.csv(dpdf,paste0("./results/sim.results.",j,"/pred_",seed,".csv"),
          row.names = FALSE)
write.csv(allEffects,paste0("./results/sim.results.",j,"/effects_",seed,".csv"),
          row.names = TRUE)
st_write(map,paste0("./results/sim.results.",j,"/maps/map_",seed,".shp"),
         append = FALSE)

sink(paste0("./results/sim.results.",j,"/texto_",seed,".txt"))
print(summary(res))
cat("Fixed Effects","\n","\n")
print(tabcrp)
cat("\n","Random Effects","\n","\n")
print(tabran)
cat("\n","Formula","\n","\n")
print(form)
cat("\n","time","\n","\n")
print(t2-t1)
sink()

