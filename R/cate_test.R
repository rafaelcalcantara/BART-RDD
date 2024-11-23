set.seed(7)
setwd("~/Git/BART-RDD/")
library(XBART)
source("R/tree_depth.R")
### Setup--------------------------------------------------
par(bty="L")
n <- 500
c <- 0
Owidth        <- 0.2
Omin          <- 5
Opct          <- 0.7
ntrees_con    <- 10
ntrees_mod    <- 10
ntrees        <- 5
Nmin          <- 5
max_depth     <- 50
num_sweeps    <- 150
burnin        <- 50
p_categorical <- 0
###
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w) w*sd(mu0.x(x))/sd(w)
tau0.x <- function(x,c) (-exp(x)+exp(c))*0.1/(1+exp(2*x))
tau0.w <- function(w,tau.bar,level) {
  w <- as.integer(w)
  if (level == 1) het <- w < mean(w)
  if (level == 2) het <- w %% 2 == 0
  if (level == 3) het <- w %% 2 == 0 | w %% 3 == 0
  if (level == 4) het <- (w %% 2 == 0 | w %% 3 == 0) & w %% 5 == 0
  if (level == 5) het <- (w %% 2 == 0) + (w %% 3 == 0) + (w %% 5 == 0)
  # het <- w
  het <- het-mean(het)
  # k <- 2*tau.bar/classes
  k <- 0.7*tau.bar/abs(min(het))
  return(k*het)
  # het <- cut(w,quantile(w,probs=seq(0,1,length.out=classes+1)),include.lowest=T)
  # het <- as.numeric(het)
  # het <- het-mean(het)
  # return(het/classes)
}
mu <- function(x,w,kappa,delta,tau) (mu0.x(x) + kappa*mu0.w(w))/sd(mu0.x(x) + kappa*mu0.w(w))*delta*sd(tau)
tau <- function(x,w,kappa,tau.bar,c,classes,level) tau.bar + tau0.x(x,c) + kappa*tau0.w(w,tau.bar,level)
#### Parameters
ate <- 0.2
delta_mu <- 0.5
delta_sig <- 0.5
kappa <- 1
p <- 0.4
#####
# pdf("Figures/sim_example_simple.pdf")
# par(mfrow=c(2,3))
### 1)-----------------------------------------------------
#### Parameters
classes <- 4
level <- 2
#### Data
x <- 2*rbeta(n,2,4)-0.75
z <- as.numeric(x>=c)
w <- runif(n,0,classes+1)
t <- tau(x,w,kappa,ate,c,classes,level)
y <- mu(x,w,kappa,delta_mu,t) + t*z
Ey <- y
y <- y + rnorm(n,0,delta_sig*sd(t))
# table(t)
# plot(x,mu(x,w,kappa,delta_mu,t))
# plot(x,t)
# plot(x,y)
## Fit models
cate <- tau(c,w,kappa,ate,c,classes,level)
test <- -Owidth+c <= x & x <= c+Owidth
cate.test <- cate[test]
mean(cate)
### BART-RDD
Opct <- 1
fit <- XBART::XBCF.rd(y, w, x, c,
                      Owidth = Owidth, Omin = Omin, Opct = Opct,
                      num_trees_mod = ntrees_mod,
                      num_trees_con = ntrees_con,
                      num_cutpoints = n,
                      num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical, max_depth = max_depth,
                      tau_con = 2*var(y)/ntrees,
                      tau_mod = 0.5*var(y)/ntrees,parallel=T,nthread=10)
###
# trees <- jsonlite::fromJSON(fit$tree_json_con)
print("BART-RDD")
# br <- rowMeans(sapply(trees$trees, function(i) sapply(i,nleaf)))
#### Plot
pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
par(mfrow=c(2,2))
# pdf("Figures/cate_bart_rdd_1.pdf")
matplot(w[test],cbind(cate.test,rowMeans(pred.bart.rdd)),
        col=c("black","blue"),ylab=bquote(tau(X==c,W)),xlab="W",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L",main="BART-RDD")
# legend("topleft",bty="n",col=c("black","blue"),pch=19,
#        legend=c(bquote(tau(X==c,W)),"BART-RDD fit"),cex=1)
# dev.off()
### S-BART
fit <- XBART::XBART(y, cbind(x,w,z), num_trees = ntrees,
                    num_cutpoints = n, num_sweeps = num_sweeps,
                    burnin = burnin, Nmin = Nmin,
                    p_categorical = p_categorical+1, max_depth = max_depth,
                    tau = var(y)/ntrees,parallel=T,nthread=10)
test.sample <- cbind(c,w,z)
###
trees <- jsonlite::fromJSON(fit$tree_json)
print("S-BART")
sb <- rowMeans(sapply(trees$trees, function(i) sapply(i,nleaf)))
#### Plot
test.sample.1 <- cbind(c,w,1)[test,]
test.sample.0 <- cbind(c,w,0)[test,]
pred1 <- XBART::predict.XBART(fit,test.sample.1)[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit,test.sample.0)[,(burnin+1):num_sweeps]
pred.sbart <- pred1-pred0
# pdf("Figures/cate_sbart_1.pdf")
matplot(w[test],cbind(cate.test,rowMeans(pred.sbart)),
        col=c("black","blue"),ylab=bquote(tau(X==c,w)),xlab="W",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L",main="S-BART")
# legend("topleft",bty="n",col=c("black","blue"),pch=19,
#        legend=c(bquote(tau(X==c,W)),"S-BART fit"),cex=1)
# dev.off()
### T-BART
fit1 <- XBART::XBART(y[z==1], cbind(x,w)[z==1,], num_trees = ntrees,
                     num_cutpoints = sum(z==1), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical, max_depth = max_depth,
                     tau = var(y[z==1])/ntrees,parallel=T,nthread=10)
fit0 <- XBART::XBART(y[z==0], cbind(x,w)[z==0,], num_trees = ntrees,
                     num_cutpoints = sum(z==0), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical, max_depth = max_depth,
                     tau = var(y[z==0])/ntrees,parallel=T,nthread=10)
test.sample <- cbind(c,w)[test,]
###
trees0 <- jsonlite::fromJSON(fit$tree_json)
trees1 <- jsonlite::fromJSON(fit$tree_json)
print("T-BART")
t0 <- rowMeans(sapply(trees0$trees, function(i) sapply(i,nleaf)))
t1 <- rowMeans(sapply(trees1$trees, function(i) sapply(i,nleaf)))
#### Plot
pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
pred.tbart <- pred1-pred0
# pdf("Figures/cate_tbart_1.pdf")
matplot(w[test],cbind(cate.test,rowMeans(pred.tbart)),
        col=c("black","blue"),ylab=bquote(tau(X==c,W)),xlab="W",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L",main="T-BART")
# legend("topleft",bty="n",col=c("black","blue"),pch=19,
#        legend=c(bquote(tau(X==c,W)),"T-BART fit"),cex=1)
# dev.off()
# mean((rowMeans(pred.bart.rdd)-cate.test)^2)
# mean((rowMeans(pred.sbart)-cate.test)^2)
# mean((rowMeans(pred.tbart)-cate.test)^2)
# plot(cate.test,rowMeans(pred.bart.rdd))
# plot(cate.test,rowMeans(pred.sbart))
# plot(cate.test,rowMeans(pred.tbart))
# ### 2)-----------------------------------------------------
# #### Parameters
# classes <- 4
# level <- 2
# #### Data
# x <- 2*rbeta(n,2,4)-0.75
# z <- as.numeric(x>=c)
# w <- runif(n,0,classes+1)
# t <- tau(x,w,kappa,ate,c,classes,level)
# y <- mu(x,w,kappa,delta_mu,t) + t*z
# Ey <- y
# y <- y + rnorm(n,0,delta_sig*sd(t))
# # table(t)
# # plot(x,mu(x,w,kappa,delta_mu,t))
# # plot(x,t)
# # plot(x,y)
# ## Fit models
# cate <- tau(c,w,kappa,ate,c,classes,level)
# test <- -Owidth+c <= x & x <= c+Owidth
# cate.test <- cate[test]
# mean(cate)
# ### BART-RDD
# fit <- XBART::XBCF.rd(y, w, x, c,
#                       Owidth = Owidth, Omin = Omin, Opct = Opct,
#                       num_trees_mod = ntrees,
#                       num_trees_con = ntrees,
#                       num_cutpoints = n,
#                       num_sweeps = num_sweeps,
#                       burnin = burnin, Nmin = Nmin,
#                       p_categorical_con = p_categorical,
#                       p_categorical_mod = p_categorical,
#                       tau_con = 2*var(y)/ntrees,
#                       tau_mod = 0.5*var(y)/ntrees,parallel=T,nthread=10)
# #### Plot
# pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
# pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
# # pdf("Figures/cate_bart_rdd_1.pdf")
# matplot(w[test],cbind(cate.test,rowMeans(pred.bart.rdd)),
#         col=c("black","blue"),ylab=bquote(tau(X==c,W)),xlab="W",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L",main="BART-RDD")
# # legend("topleft",bty="n",col=c("black","blue"),pch=19,
# #        legend=c(bquote(tau(X==c,W)),"BART-RDD fit"),cex=1)
# # dev.off()
# ### S-BART
# fit <- XBART::XBART(y, cbind(x,w,z), num_trees = ntrees,
#                     num_cutpoints = n, num_sweeps = num_sweeps,
#                     burnin = burnin, Nmin = Nmin,
#                     p_categorical = p_categorical+1,
#                     tau = var(y)/ntrees,parallel=T,nthread=10)
# test.sample <- cbind(x,w,z)
# #### Plot
# test.sample.1 <- cbind(x,w,1)[test,]
# test.sample.0 <- cbind(x,w,0)[test,]
# pred1 <- XBART::predict.XBART(fit,test.sample.1)[,(burnin+1):num_sweeps]
# pred0 <- XBART::predict.XBART(fit,test.sample.0)[,(burnin+1):num_sweeps]
# pred.sbart <- pred1-pred0
# # pdf("Figures/cate_sbart_1.pdf")
# matplot(w[test],cbind(cate.test,rowMeans(pred.sbart)),
#         col=c("black","blue"),ylab=bquote(tau(X==c,w)),xlab="W",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L",main="S-BART")
# # legend("topleft",bty="n",col=c("black","blue"),pch=19,
# #        legend=c(bquote(tau(X==c,W)),"S-BART fit"),cex=1)
# # dev.off()
# ### T-BART
# fit1 <- XBART::XBART(y[z==1], cbind(x,w)[z==1,], num_trees = ntrees,
#                      num_cutpoints = sum(z==1), num_sweeps = num_sweeps,
#                      burnin = burnin, Nmin = Nmin,
#                      p_categorical = p_categorical,
#                      tau = var(y[z==1])/ntrees,parallel=T,nthread=10)
# fit0 <- XBART::XBART(y[z==0], cbind(x,w)[z==0,], num_trees = ntrees,
#                      num_cutpoints = sum(z==0), num_sweeps = num_sweeps,
#                      burnin = burnin, Nmin = Nmin,
#                      p_categorical = p_categorical,
#                      tau = var(y[z==0])/ntrees,parallel=T,nthread=10)
# test.sample <- cbind(x,w)[test,]
# #### Plot
# pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
# pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
# pred.tbart <- pred1-pred0
# # pdf("Figures/cate_tbart_1.pdf")
# matplot(w[test],cbind(cate.test,rowMeans(pred.tbart)),
#         col=c("black","blue"),ylab=bquote(tau(X==c,W)),xlab="W",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L",main="T-BART")
# # legend("topleft",bty="n",col=c("black","blue"),pch=19,
# #        legend=c(bquote(tau(X==c,W)),"T-BART fit"),cex=1)
# # dev.off()
# ######
# dev.off()
# ###
boxplot(cbind(`BART-RDD`=br,`S-BART`=sb,`T-BART0`=t0,`T-BART1`=t1))