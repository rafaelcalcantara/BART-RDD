set.seed(0)
setwd("~/../Git/BART-RDD/")
library(XBART)
###
n <- 1000
c <- 0
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w) {
  if(is.list(w)) sapply(w, function(i) rowMeans(i))
  else rowMeans(w)
}
tau0.x <- function(x) -0.1*x
tau0.w <- function(w) {
  if(is.list(w)) sapply(w, function(i) rowMeans(i[,3:4]))
  else rowMeans(w[,3:4])
}
mu <- function(x,w,kappa,delta) (mu0.x(x) + kappa*mu0.w(w))/sd(mu0.x(x) + kappa*mu0.w(w))*delta
tau <- function(x,w,kappa,delta,tau.bar) tau.bar + (tau0.x(x) + kappa*tau0.w(w))/sd(tau0.x(x) + kappa*tau0.w(w))*delta
#### Parameters
ate <- 0.5
delta_mu <- 1
delta_tau <- 5
kappa <- 1
#### Data
x <- 2*rbeta(n,2,4)-0.75
z <- as.numeric(x>=c)
w <- cbind(runif(n,-0.1,0.1),rnorm(n,0,0.2),rbinom(n,1,0.4)-0.4,rbinom(n,1,dnorm(x,c,0.5)))
w[,4] <- w[,4]-mean(w[,4])
y <- mu(x,w,kappa,delta_mu) + tau(x,w,kappa,delta_tau,ate)*z
Ey <- y
y <- y + rnorm(n)
#### Checking data
cate.x <- tau(x,w,kappa,delta_tau,ate)
cate <- tau(0,w,kappa,delta_tau,ate)
plot(x,cate.x,bty="n",pch=19)
## Fit models
Owidth        <- 0.1
Omin          <- 10
Opct          <- 0.6
ntrees        <- 5
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 2
test          <- -Owidth+c <= x & x <= c+Owidth
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
                      tau_mod = 0.5*var(y)/ntrees)
#### Plots
##### Full data
pred <- XBART::predict.XBCFrd(fit,w,x)
pred.bart.rdd <- pred$yhats.adj[,(burnin+1):num_sweeps]
matplot(x,cbind(Ey,rowMeans(pred.bart.rdd)),
        col=c("black","blue"),pch=19,bty="n",ylab="E[Y|X,W,Z]",xlab="X",
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
abline(v=Owidth,lty=2)
abline(v=-Owidth,lty=2)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c("E[Y|X,W,Z]","BART-RDD fit"),cex=0.75)
##### Mu
pred.bart.rdd <- pred$mu.adj[,(burnin+1):num_sweeps]
matplot(x,cbind(mu(x,w,kappa,delta_mu),rowMeans(pred.bart.rdd)),bty="n",
        col=c("black","blue"),ylab=bquote(mu),xlab="X",pch=19,
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
abline(v=Owidth,lty=2)
abline(v=-Owidth,lty=2)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(mu(X,W)),"BART-RDD fit"),cex=0.75)
##### Tau
pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
matplot(x,cbind(tau(x,w,kappa,delta_tau,ate),rowMeans(pred.bart.rdd)),
        col=c("black","blue"),bty="n",ylab=bquote(tau),xlab="X",pch=19,
        cex.axis=0.75,cex.lab=0.75,
        ylim=c(min(cbind(tau(x,w,kappa,delta_tau,ate),rowMeans(pred.bart.rdd))),
               max(cbind(tau(x,w,kappa,delta_tau,ate),rowMeans(pred.bart.rdd)))+3))
abline(v=c,lty=3)
abline(v=Owidth,lty=2)
abline(v=-Owidth,lty=2)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X,W)),"BART-RDD fit"),cex=0.75)
####
pred <- XBART::predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
pred.bart.rdd <- pred$tau.adj[,(burnin+1):num_sweeps]
####
matplot(cbind(sort(cate.test),rowMeans(pred.bart.rdd)[order(cate.test)]),
        col=c("black","blue"),bty="n",ylab=bquote(tau),xlab="i",pch=19,
        cex.axis=0.75,cex.lab=0.75)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X==c,W)),"BART-RDD fit"),cex=0.75)
### S-BART
fit <- XBART::XBART(y, cbind(x,w,z), num_trees = ntrees,
                    num_cutpoints = n, num_sweeps = num_sweeps,
                    burnin = burnin, Nmin = Nmin,
                    p_categorical = p_categorical,
                    tau = var(y)/ntrees, parallel=F)
test.sample <- cbind(x,w,z)
#### Plots
##### Full data
pred.sbart <- XBART::predict.XBART(fit,test.sample)[,(burnin+1):num_sweeps]
matplot(x,cbind(Ey,rowMeans(pred.sbart)),
        col=c("black","blue","blue"),pch=19,bty="n",ylab="E[Y|X,W,Z]",xlab="X",
        ylim=c(min(cbind(Ey,pred.sbart)),max(cbind(Ey,pred.sbart))+5),
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
legend("topleft",bty="n",col=c("black","blue"),pch=19,legend=c("E[Y|X,W,Z]","S-BART fit"),cex=0.75)
#### mu
test.sample <- cbind(x,w,0)
pred.sbart <- XBART::predict.XBART(fit,test.sample)[,(burnin+1):num_sweeps]
matplot(x,cbind(mu(x,w,kappa,delta_mu),rowMeans(pred.sbart)),
        col=c("black","blue","blue"),pch=19,bty="n",ylab=bquote(mu),xlab="X",
        ylim=c(min(cbind(Ey,pred.sbart)),max(cbind(Ey,pred.sbart))+5),
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(mu(X,W)),"S-BART fit"),cex=0.75)
#### tau
test.sample.1 <- cbind(x,w,1)
test.sample.0 <- cbind(x,w,0)
pred1 <- XBART::predict.XBART(fit,test.sample.1)[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit,test.sample.0)[,(burnin+1):num_sweeps]
pred.sbart <- pred1-pred0
matplot(x,cbind(tau(x,w,kappa,delta_tau,ate),rowMeans(pred.sbart)),
        col=c("black","blue","blue"),pch=19,bty="n",ylab=bquote(tau),xlab="X",
        ylim=c(min(cbind(Ey,pred.sbart)),max(cbind(Ey,pred.sbart))+5),
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X,W)),"S-BART fit"),cex=0.75)
####
test.sample.1 <- cbind(0,w,1)
test.sample.0 <- cbind(0,w,0)
pred1 <- XBART::predict.XBART(fit,test.sample.1[test,])[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit,test.sample.0[test,])[,(burnin+1):num_sweeps]
pred.sbart <- pred1-pred0
####
matplot(cbind(sort(cate.test),rowMeans(pred.sbart)[order(cate.test)]),
        col=c("black","blue"),bty="n",ylab=bquote(tau),xlab="i",pch=19,
        cex.axis=0.75,cex.lab=0.75)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X==c,W)),"S-BART fit"),cex=0.75)
### T-BART
fit1 <- XBART::XBART(y[z==1], cbind(x,w)[z==1,], num_trees = ntrees,
                     num_cutpoints = sum(z==1), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==1])/ntrees, parallel=F)
fit0 <- XBART::XBART(y[z==0], cbind(x,w)[z==0,], num_trees = ntrees,
                     num_cutpoints = sum(z==0), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==0])/ntrees, parallel=F)
test.sample <- cbind(x,w)
#### Plots
##### Full data
pred1 <- XBART::predict.XBART(fit1,test.sample[z==1,])[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,test.sample[z==0,])[,(burnin+1):num_sweeps]
pred.tbart <- c(rowMeans(pred0),rowMeans(pred1))
matplot(x,cbind(Ey,pred.tbart),
        col=c("black","blue","blue"),pch=19,bty="n",ylab="E[Y|X,W,Z]",xlab="X",
        ylim=c(min(cbind(Ey,pred.tbart)),max(cbind(Ey,pred.tbart))+5),
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
legend("topleft",bty="n",col=c("black","blue"),pch=19,legend=c("E[Y|X,W,Z]","T-BART fit"),cex=0.75)
#### mu
pred1 <- XBART::predict.XBART(fit0,test.sample[z==1,])[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,test.sample[z==0,])[,(burnin+1):num_sweeps]
pred.tbart <- c(rowMeans(pred0),rowMeans(pred1))
matplot(x,cbind(mu(x,w,kappa,delta_mu),pred.tbart),
        col=c("black","blue","blue"),pch=19,bty="n",ylab=bquote(mu),xlab="X",
        ylim=c(min(cbind(Ey,pred.tbart)),max(cbind(Ey,pred.tbart))+5),
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(mu(X,W)),"T-BART fit"),cex=0.75)
#### tau
pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
pred.tbart <- pred1-pred0
matplot(x,cbind(tau(x,w,kappa,delta_tau,ate),rowMeans(pred.tbart)),
        col=c("black","blue","blue"),pch=19,bty="n",ylab=bquote(tau),xlab="X",
        ylim=c(min(cbind(Ey,pred.tbart)),max(cbind(Ey,pred.tbart))+5),
        cex.axis=0.75,cex.lab=0.75)
abline(v=c,lty=3)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X,W)),"T-BART fit"),cex=0.75)
####
pred1 <- XBART::predict.XBART(fit1,test.sample[test,])[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,test.sample[test,])[,(burnin+1):num_sweeps]
pred.tbart <- pred1-pred0
####
matplot(cbind(sort(cate.test),rowMeans(pred.tbart)[order(cate.test)]),
        col=c("black","blue"),bty="n",ylab=bquote(tau),xlab="i",pch=19,
        cex.axis=0.75,cex.lab=0.75)
legend("topleft",bty="n",col=c("black","blue"),pch=19,
       legend=c(bquote(tau(X==c,W)),"T-BART fit"),cex=0.75)
### Comparing
#### CATE predictions
boxplot(cbind(`BART-RDD`=rowMeans(pred.bart.rdd)-cate.test,
              `S-BART`=rowMeans(pred.sbart)-cate.test,
              `T-BART`=rowMeans(pred.tbart)-cate.test),
        cex.axis=0.75,cex.lab=0.75,main="CATE bias")
abline(h=0,lty=2)
boxplot(cbind(`BART-RDD`=rowMeans(pred.bart.rdd)-cate.test,
              `S-BART`=rowMeans(pred.sbart)-cate.test,
              `T-BART`=rowMeans(pred.tbart)-cate.test),
        cex.axis=0.75,cex.lab=0.75,main="CATE variance")
####
mean(sqrt((rowMeans(pred.bart.rdd)-cate.test)^2))
mean(sqrt((rowMeans(pred.tbart)-cate.test)^2))
mean(sqrt((rowMeans(pred.sbart)-cate.test)^2))
abs(mean(colMeans(pred.bart.rdd))-mean(cate.test))
abs(mean(colMeans(pred.tbart))-mean(cate.test))
abs(mean(colMeans(pred.sbart))-mean(cate.test))
mean(cate.test)
# summary(colMeans(pred.bart.rdd))
# summary(colMeans(pred.tbart))
# aggregate(rowMeans(pred.bart.rdd)-cate.test,list(CATE=cate.test),mean)
# aggregate(rowMeans(pred.tbart)-cate.test,list(CATE=cate.test),mean)
# mean(rowMeans(pred.bart.rdd)-cate.test)
# mean(rowMeans(pred.tbart)-cate.test)