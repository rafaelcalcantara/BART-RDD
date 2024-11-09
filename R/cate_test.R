set.seed(7)
setwd("~/Git/BART-RDD/")
library(XBART)
### Setup--------------------------------------------------
par(bty="L")
n <- 1000
c <- 0
Owidth        <- 0.2
Omin          <- 1
Opct          <- 0.9
ntrees        <- 5
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 2
###
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w) w*sd(mu0.x(x))/sd(w)
tau0.x <- function(x) -0.1*x
tau0.w <- function(w,tau.bar) {
  het <- w < 5
  tau.bar*het*sd(tau0.x(x))/sd(het)
}
mu <- function(x,w,kappa,delta) (mu0.x(x) + kappa*mu0.w(w))/sd(mu0.x(x) + kappa*mu0.w(w))*delta
tau <- function(x,w,kappa,delta,tau.bar) tau.bar + (tau0.x(x) + kappa*tau0.w(w,tau.bar))/sd(tau0.x(x) + kappa*tau0.w(w,tau.bar))*delta
### 1)-----------------------------------------------------
#### Parameters
ate <- 0.2
delta_mu <- 0.5
delta_tau <- 0.3
kappa <- 1
classes <- 10
p <- 0.4
#### Data
x <- 2*rbeta(n,2,4)-0.75
z <- as.numeric(x>=c)
test <- -0.07 <= x & x <= 0.07
w <- rbinom(n,classes,p) + 1
y <- mu(x,w,kappa,delta_mu) + tau(x,w,kappa,delta_tau,ate)*z
Ey <- y
y <- y + rnorm(n)
## Fit models
cate <- tau(c,w,kappa,delta_tau,ate)
test <- -Owidth+c <= x & x <= c+Owidth
cate.test <- cate[test]
### BART-RDD
fit <- XBART::XBCF.rd(y, w, x, c,
                      Owidth = Owidth, Omin = Omin, Opct = Opct,
                      num_trees_mod = ntrees,
                      num_trees_con = ntrees,
                      num_cutpoints = n,
                      num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical,
                      tau_con = 2*var(y)/ntrees,
                      tau_mod = 0.5*var(y)/ntrees,parallel=T,nthread=10)
#### Plot
pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
# pdf("Figures/cate_bart_rdd_1.pdf")
matplot(cbind(sort(cate.test),rowMeans(pred.bart.rdd)[order(cate.test)]),
        col=c("black","blue"),ylab=bquote(tau),xlab="i",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L")
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X==c,W)),"BART-RDD fit"),cex=1)
# dev.off()
### S-BART
fit <- XBART::XBART(y, cbind(x,w,z), num_trees = ntrees,
                    num_cutpoints = n, num_sweeps = num_sweeps,
                    burnin = burnin, Nmin = Nmin,
                    p_categorical = p_categorical,
                    tau = var(y)/ntrees,parallel=T,nthread=10)
test.sample <- cbind(x,w,z)
#### Plot
test.sample.1 <- cbind(x,w,1)[test,]
test.sample.0 <- cbind(x,w,0)[test,]
pred1 <- XBART::predict.XBART(fit,test.sample.1)[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit,test.sample.0)[,(burnin+1):num_sweeps]
pred.sbart <- pred1-pred0
# pdf("Figures/cate_sbart_1.pdf")
matplot(cbind(sort(cate.test),rowMeans(pred.sbart)[order(cate.test)]),
        col=c("black","blue"),ylab=bquote(tau),xlab="i",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L")
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X==c,W)),"S-BART fit"),cex=1)
# dev.off()
### T-BART
fit1 <- XBART::XBART(y[z==1], cbind(x,w)[z==1,], num_trees = ntrees,
                     num_cutpoints = sum(z==1), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==1])/ntrees,parallel=T,nthread=10)
fit0 <- XBART::XBART(y[z==0], cbind(x,w)[z==0,], num_trees = ntrees,
                     num_cutpoints = sum(z==0), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==0])/ntrees,parallel=T,nthread=10)
test.sample <- cbind(x,w)[test,]
#### Plot
pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
pred.tbart <- pred1-pred0
# pdf("Figures/cate_tbart_1.pdf")
matplot(cbind(sort(cate.test),rowMeans(pred.tbart)[order(cate.test)]),
        col=c("black","blue"),ylab=bquote(tau),xlab="i",pch=19,
        cex.axis=1,cex.lab=1,cex=1,bty="L")
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X==c,W)),"T-BART fit"),cex=1)
# dev.off()
### 2)-----------------------------------------------------
#### Parameters
# ate <- 0.5
# delta_mu <- 1.25
# delta_tau <- 0.3
# kappa <- 1
# #### Data
# x <- 2*rbeta(n,2,4)-0.75
# z <- as.numeric(x>=c)
# w <- cbind(runif(n,-0.1,0.1),rnorm(n,0,0.2),rbinom(n,1,0.4)-0.4,rbinom(n,1,dnorm(x,c,0.5)))
# w[,4] <- w[,4]-mean(w[,4])
# y <- mu(x,w,kappa,delta_mu) + tau(x,w,kappa,delta_tau,ate)*z
# Ey <- y
# y <- y + rnorm(n)
# ## Fit models
# cate <- tau(c,w,kappa,delta_tau,ate)
# test <- -Owidth+c <= x & x <= c+Owidth
# cate.test <- cate[test]
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
# pred <- XBART::predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
# pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
# pdf("Figures/cate_bart_rdd_2.pdf")
# matplot(cbind(sort(cate.test),rowMeans(pred.bart.rdd)[order(cate.test)]),
#         col=c("black","blue"),ylab=bquote(tau),xlab="i",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L")
# legend("topleft",bty="n",col=c("black","blue"),pch=19,
#        legend=c(bquote(tau(X==c,W)),"BART-RDD fit"),cex=1)
# dev.off()
# ### S-BART
# fit <- XBART::XBART(y, cbind(x,w,z), num_trees = ntrees,
#                     num_cutpoints = n, num_sweeps = num_sweeps,
#                     burnin = burnin, Nmin = Nmin,
#                     p_categorical = p_categorical,
#                     tau = var(y)/ntrees,parallel=T,nthread=10)
# test.sample <- cbind(x,w,z)
# #### Plot
# test.sample.1 <- cbind(x,w,1)[test,]
# test.sample.0 <- cbind(x,w,0)[test,]
# pred1 <- XBART::predict.XBART(fit,test.sample.1)[,(burnin+1):num_sweeps]
# pred0 <- XBART::predict.XBART(fit,test.sample.0)[,(burnin+1):num_sweeps]
# pred.sbart <- pred1-pred0
# pdf("Figures/cate_sbart_2.pdf")
# matplot(cbind(sort(cate.test),rowMeans(pred.sbart)[order(cate.test)]),
#         col=c("black","blue"),ylab=bquote(tau),xlab="i",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L")
# legend("topleft",bty="n",col=c("black","blue"),pch=19,
#        legend=c(bquote(tau(X==c,W)),"S-BART fit"),cex=1)
# dev.off()
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
# pdf("Figures/cate_tbart_2.pdf")
# matplot(cbind(sort(cate.test),rowMeans(pred.tbart)[order(cate.test)]),
#         col=c("black","blue"),ylab=bquote(tau),xlab="i",pch=19,
#         cex.axis=1,cex.lab=1,cex=1,bty="L")
# legend("topleft",bty="n",col=c("black","blue"),pch=19,
#        legend=c(bquote(tau(X==c,W)),"T-BART fit"),cex=1)
# dev.off()