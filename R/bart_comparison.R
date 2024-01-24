set.seed(7)
library(XBART)
setwd("~/../Git/BART-RDD/")
### Parameters
Omin          <- 5
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 2
## Functions
mu <- function(x=NULL,w=NULL,het=NULL,nonlin=NULL)
{
  if (het==1) out <- 0.1875*sin((w[,1]+w[,2])*pi) + 1 + 1.875*x - 1.25*x^2 + 1.75*x^3 + nonlin*dnorm(x-c,0,1)*sin(2*pi*dnorm(x-c,0,0.05))
  if (het==2) out <- 0.1*(w[,4]-mean(w[,4]))*sin((w[,1]+w[,2])*pi) + 1 + 0.9*(w[,4]-mean(w[,4])) + w[,4]*x - 0.9*w[,4]*x^2 + w[,4]*x^3 + nonlin*(w[,4]-mean(w[,4]))*dnorm(x-c,0,1)*sin(2*pi*dnorm(x-c,0,0.05))
  if (het==3)
  {
    f <- rep(0,length(x))
    for (i in 1:length(f))
    {
      if (w[i,3]==1 & w[i,4]==1) f[i] <- 0.43
      if (w[i,3]==2 & w[i,4]==1) f[i] <- 0.27
      if (w[i,3]==3 & w[i,4]==1) f[i] <- 0.1
      if (w[i,3]==1 & w[i,4]==2) f[i] <- 0.77
      if (w[i,3]==2 & w[i,4]==2) f[i] <- 0.93
      if (w[i,3]==3 & w[i,4]==2) f[i] <- 1.1
    }
    out <- 0.2*(f-mean(f))*sin((w[,1]+w[,2])*pi) + 1 + 2*(f-mean(f)) + 2.27*f*x - 1.13*f*x^2 + 2*f*x^3 + nonlin*(f-mean(f))*dnorm(x-c,0,1)*sin(2*pi*dnorm(x-c,0,0.05))
  }
  return(out)
}
tau <- function(x=NULL,w=NULL,het=NULL,ate=NULL)
{
  if (het==1) out <- (2*w[,1] + w[,2]^2 + 10/(5-exp(-10*w[,2])) - 4.35) + 0.025*cos((w[,1]+w[,2])*pi) + ate - 2.8*x + 1.4*x^2 - 0.14*x^3
  if (het==2) out <- 2*(w[,4]-mean(w[,4]))*w[,1] + (w[,4]-mean(w[,4]))*w[,2]^2 + (w[,4]-mean(w[,4]))*10/(5-exp(-10*w[,2])) + 0.0125*(w[,4]-mean(w[,4]))*cos((w[,1]+w[,2])*pi) + ate + ate*(w[,4]-mean(w[,4])) - 1.8*w[,4]*x + 0.9*w[,4]*x^2 - 0.09*w[,4]*x^3
  if (het==3)
  {
    f <- rep(0,length(x))
    for (i in 1:length(f))
    {
      if (w[i,3]==1 & w[i,4]==1) f[i] <- 0.43
      if (w[i,3]==2 & w[i,4]==1) f[i] <- 0.27
      if (w[i,3]==3 & w[i,4]==1) f[i] <- 0.1
      if (w[i,3]==1 & w[i,4]==2) f[i] <- 0.77
      if (w[i,3]==2 & w[i,4]==2) f[i] <- 0.93
      if (w[i,3]==3 & w[i,4]==2) f[i] <- 1.1
    }
    out <- 2*(f-mean(f))*w[,1] + (f-mean(f))*w[,2]^2 + (f-mean(f))*10/(5-exp(-10*w[,2])) + 0.1*(f-mean(f))*cos((w[,1]+w[,2])*pi) + ate + 1.5*ate*(f-mean(f)) - 3.4*f*x + 1.7*f*x^2 - 0.17*f*x^3
  }
  return(out)
}
##
n <- 1000
s <- 1
c <- 0
t <- 0.1
het <- 1
nonlin <- 1
u <- matrix(runif(n*s),n,s)
## x <- 2*rbeta(n,2,4)-u+0.2
## x <- rnorm(n,c-0.5+u,0.5)
## x <- runif(n,c-2,c+2)
## x <- 2*rbeta(n,2,4)-1
## w1 <- matrix(runif(n*s,as.numeric(u),as.numeric(u+1)),n,s)
## w1 <- 0.5 - x
w1 <- matrix(runif(n*s,0,1),n,s)
w2 <- matrix(runif(n*s,0,2),n,s)
## w3 <- matrix(rbinom(n*s,2,as.numeric(u))+1,n,s)
w3 <- matrix(rbinom(n*s,2,0.7),n,s)
w4 <- matrix(rbinom(n*s,1,0.6)+1,n,s)
w <- cbind(w1,w2,w3,w4)
p <- function(w) pnorm(rowSums(w),1.5,0.5)
x <- qnorm(p(cbind(w1,w2)),0,0.2)
z <- as.numeric(x>c)
y <- mu(x,w,het,nonlin) + tau(x,w,het,t)*z + rnorm(n,0,0.5*sd(mu(x,w,het,nonlin)+tau(x,w,het,t)*z))
h <- 0.5*sd(x)
xtest <- x[c-2*sd(x)<x & x<c+2*sd(x)]
xtest <- c(aggregate(xtest[xtest<c],list(cut(xtest[xtest<c],2)),mean)[,2],
           c,
           aggregate(xtest[c<xtest],list(cut(xtest[c<xtest],2)),mean)[,2])
ate <- sapply(xtest,function(i) mean(tau(rep(i,n),w,het,t)))
### Plot data
summary(tau(rep(c,n),w,het,t))
cor(x,tau(x,w,het,t))
cor(w1,tau(x,w,het,t))
cor(w2,tau(x,w,het,t))
cor(w3,tau(x,w,het,t))
cor(w4,tau(x,w,het,t))
cor(x,w1)
par(mfrow=c(2,3))
plot(x,tau(x,w,het,t))
plot(w1,tau(x,w,het,t))
plot(w2,tau(x,w,het,t))
plot(w3,tau(x,w,het,t))
plot(w4,tau(x,w,het,t))
plot(x,w1)
## BART-RDD
fit <- XBART::XBCF.rd(y, w, x, c,
                      Owidth = h, Omin = Omin, Opct = Opct,
                      num_trees_mod = ntrees,
                      num_trees_con = ntrees,
                      num_cutpoints = n,
                      num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical,
                      tau_con = 2*var(y)/ntrees,
                      tau_mod = 0.5*var(y)/ntrees,
                      parallel = T, nthread = 11)
bart.rdd <- sapply(xtest,function(i) colMeans(XBART::predict.XBCFrd(fit,w,rep(i,n))$tau.adj[,(burnin+1):num_sweeps]))
bart.rdd <- t(apply(bart.rdd,2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
bart.rdd.cate <- rowMeans(XBART::predict.XBCFrd(fit,w,rep(c,n))$tau.adj[,(burnin+1):num_sweeps])
## BCF
fit <- XBART::XBCF.discrete(y=y, Z=z, X_con = cbind(x,w), X_mod = cbind(x,w),
                            num_trees_mod = ntrees,
                            num_trees_con = ntrees,
                            num_cutpoints = n,
                            num_sweeps = num_sweeps,
                            burnin = burnin, Nmin = Nmin,
                            p_categorical_con = p_categorical,
                            p_categorical_mod = p_categorical,
                            tau_con = 2*var(y)/ntrees,
                            tau_mod = 0.5*var(y)/ntrees,
                            parallel = T, nthread = 11)
bcf <- sapply(xtest,function(i) colMeans(XBART::predict.XBCFdiscrete(fit,X_con = cbind(i,w), X_mod = cbind(i,w),Z=z,pihat=z,burnin=burnin)$tau.adj[,(burnin+1):num_sweeps]))
bcf <- t(apply(bcf,2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
bcf.cate <- rowMeans(XBART::predict.XBCFdiscrete(fit,X_con = cbind(c,w), X_mod = cbind(c,w),Z=z,pihat=z,burnin=burnin)$tau.adj[,(burnin+1):num_sweeps])
## SBART
fit <- XBART::XBART(y, cbind(x,w,z), num_trees = ntrees,
                     num_cutpoints = n, num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical+1,
                     tau = var(y)/ntrees, parallel=T,nthread = 11)
sbart.pred <- function(i)
{
  pred1 <- XBART::predict.XBART(fit,cbind(i,w,rep(1,n)))[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit,cbind(i,w,rep(0,n)))[,(burnin+1):num_sweeps]
  colMeans(pred1-pred0)
}
sbart <- sapply(xtest,sbart.pred)
sbart <- t(apply(sbart,2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
pred1 <- XBART::predict.XBART(fit,cbind(c,w,rep(1,n)))[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit,cbind(c,w,rep(0,n)))[,(burnin+1):num_sweeps]
sbart.cate <- rowMeans(pred1-pred0)
## TBART
fit1 <- XBART::XBART(y[z==1], cbind(x,w)[z==1,], num_trees = ntrees,
                     num_cutpoints = sum(z==1), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==1])/ntrees, parallel=T, nthread = 11)
fit0 <- XBART::XBART(y[z==0], cbind(x,w)[z==0,], num_trees = ntrees,
                     num_cutpoints = sum(z==0), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==0])/ntrees, parallel=T, nthread = 11)
tbart.pred <- function(i)
{
  pred1 <- XBART::predict.XBART(fit1,cbind(i,w))[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit0,cbind(i,w))[,(burnin+1):num_sweeps]
  colMeans(pred1-pred0)
}
tbart <- sapply(xtest,tbart.pred)
tbart <- t(apply(tbart,2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
pred1 <- XBART::predict.XBART(fit1,cbind(c,w))[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,cbind(c,w))[,(burnin+1):num_sweeps]
tbart.cate <- rowMeans(pred1-pred0)
## Plots
par(mfrow=c(2,3))
## BART-RDD
matplot(xtest,cbind(ate,bart.rdd),type=c("l","b","l","l"),lty=c(1,1,2,2),bty="n",pch=20,col=c(1,2,2,2),main="BART-RDD",xlab="X",ylab="")
abline(v=c,lty=2)
abline(v=c-h,lty=2)
abline(v=c+h,lty=2)
## BCF
matplot(xtest,cbind(ate,bcf),type=c("l","b","l","l"),lty=c(1,1,2,2),bty="n",pch=20,col=c(1,2,2,2),main="BCF",xlab="X",ylab="")
abline(v=c,lty=2)
##
plot.new()
legend("center",legend=c("True","Fit"),col=1:2,lwd=2,lty=1)
## SBART
matplot(xtest,cbind(ate,sbart),type=c("l","b","l","l"),lty=c(1,1,2,2),bty="n",pch=20,col=c(1,2,2,2),main="S-BART",xlab="X",ylab="")
abline(v=c,lty=2)
## TBART
matplot(xtest,cbind(ate,tbart),type=c("l","b","l","l"),lty=c(1,1,2,2),bty="n",pch=20,col=c(1,2,2,2),main="T-BART",xlab="X",ylab="")
abline(v=c,lty=2)
## CATE
cate.error <- (cbind(bart.rdd.cate,bcf.cate,sbart.cate,tbart.cate)-tau(rep(c,n),w,het,t))^2
colnames(cate.error) <- c("BART-RDD","BCF","S-BART","T-BART")
boxplot(cate.error,main="CATE Error")
##
### ATE Interval Length
sapply(list(bart.rdd,bcf,sbart,tbart),function(i) i[3,3]-i[3,2])
### CATE RMSE
sqrt(colMeans(cate.error))