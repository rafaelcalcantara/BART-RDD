# Setup------------------------------------------------------------------------
library(stochtree)
## DGP parameters
n <- 2000
ate <- 0.25
k1 <- 8 ## variability in mu0.x
k2 <- 1 ## amplitude of tau0.w relative to ATE (set from 0 to 1)
# k3 <- 1 ## sd of mu0.w relative to sd of tau
sig_error <- 0.25 ## relative to sd of tau
p <- 2 # Dim of w
rho <- 0
c <- 0
pts_in_window <- 150
x.center <- 0 ## In case we want to center x at a different location
print.params <- T ## sanity check
## DGP functions
mu0.x <- function(x,c) 0.1*(x-c+1)^3
mu0.w <- function(w) (rowMeans(w)+2)^2
# tau0.x <- function(x,c) 0 # cheating to make this constant?
# tau0.w <- function(w,k) k*(2*pnorm(rowMeans(w),0,0.3)-1)
mu <- function(x,w,k1,k2,c) k1*mu0.x(x,c)+k2*mu0.w(w)
# mu0 <- function(x,w,k1,k3,c) k1*mu0.x(x,c) + k3*mu0.w(w)
# tau0 <- function(x,w,k2,c) tau0.w(w,k2) + tau0.x(x,c)
tau0 <- function(w,ate) ate*(pnorm(rowMeans(w)+2,0,1) + pnorm(rowMeans(w)+1,0,1)/4 + dnorm(2*rowMeans(w),0,1) + dnorm(rowMeans(w)-1,0,1))
# mu <- function(x,w,k1,k3,mu.bar,c) mu0(x,w,k1,k3,c)
# tau <- function(x,w,ate,k2,tau.bar,c) tau0(x,w,k2,c) - tau.bar + ate
tau <- function(w,ate,tau.bar) tau0(w,ate) - tau.bar + ate
h.grid <- function(x,c,grid)
{
  abs.x <- sort(abs(x-c))
  out <- rep(0,length(grid))
  names(out) <- grid
  x.right <- sum(c < x)
  x.left <- sum(x < c)
  x.tot <- length(x)
  for(total in grid)
  {
    i <- 1
    sum.right <- sum.left <- 0
    while(sum.right < total | sum.left < total) 
    {
      sum.left <- sum(c-abs.x[i] <= x & x < c)
      sum.right <- sum(c < x & x <= c+abs.x[i])
      if (sum.left == sum(x<c) & sum.right == sum(c<x)) break
      i <- i+1
    }
    out[as.character(total)] <- abs.x[i]
  }
  return(out)
}
## Set parameters and demean data
x0 <- rnorm(n,x.center,1)
# x0 <- 2*rbeta(n,2,4)-1
w0 <- matrix(rnorm(n*p,rho*x.center,sqrt(1-rho^2)),n,p)
k2 <- k2*sd(tau0(w0,ate))/sd(mu(c,w0,k1,1,c))
mu.bar <- mean(mu(c,w0,k1,k2,c))
tau.bar <- mean(tau0(w0,ate))
sig_error <- sig_error*max(abs(mean(tau(w0,ate,tau.bar))),2*sd(tau(w0,ate,tau.bar)))
## Samples
x <- rnorm(n,x.center,1)
# x <- 2*rbeta(n,2,4)-1
z <- as.numeric(x>=c)
w <- matrix(rnorm(n*p,rep(x,p)*rho,sqrt(1-rho^2)),n,p)
prog <- mu(x,w,k1,k2,c)
cate <- tau(w,ate,tau.bar)
y <- prog + cate*z + rnorm(n,0,sig_error)
c <- c/sd(x)
x <- x/sd(x)
# c <- 0
# z <- as.numeric(x>=c)
prog <- mu(c,w,k1,k2,c)
cate <- tau(w,ate,tau.bar)
if (print.params)
{
  print(paste(c("K1: ","K2: ","sigma: ","ATE: ","sd(mu): ","sd(tau): "),
              c(k1,k2,sig_error,mean(cate),sd(prog),sd(cate))))
}
## Plotting the data
# Owidth <- h.grid(x,c,pts_in_window)
Owidth <- 0.05
test <- -Owidth+c<=x & x<=Owidth+c
par(mfrow=c(2,2))
plot(x,y,col=z+1,pch=19,bty="n")
abline(v=c,lty=2)
plot(x[test],y[test],col=z[test]+1,pch=19,bty="n")
abline(v=c,lty=2)
plot(rowMeans(w)[test],prog[test])
plot(rowMeans(w)[test],cate[test])
# Estimation-------------------------------------------------------------------
B <- cbind(z*x, x, z,1)
B1 <- cbind(rep(c,n), rep(c,n), rep(1,n), rep(1,n))
B0 <- cbind(rep(0,n), rep(c,n), rep(0,n), rep(1,n))
### If using stochtree prior to separate param lists
parmlist <- list(num_trees_mean =50,sigma2_init = 0.01,
                 sample_sigma_global=TRUE,min_samples_leaf_mean=20, alpha_mean = 0.95, beta_mean=2,
                 standardize = T, max_depth_mean = 10)
### If using stochtree with separate param lists
global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.01)
mean.parmlist <- list(num_trees=50, min_samples_leaf=20, alpha=0.95, beta=10, max_depth=10, sample_sigma2_leaf=FALSE)
## BARDDT----------------------------------------------------------------------
### If using stochtree prior to separate param lists
# time.barddt.fit <- system.time({
#   barddt.fit = bart(X_train= as.matrix(cbind(x*(!test),w)), y_train=y,
#                     W_train = B, params=parmlist,
#                     num_mcmc=1000,num_gfr=30)
# })
### If using stochtree with separate param lists
time.barddt.fit <- system.time({
  barddt.fit = bart(X_train= as.matrix(cbind(x*(!test),w)), y_train=y,
                    W_train = B, mean_forest_params=mean.parmlist,
                    general_params=global.parmlist,
                    num_mcmc=1000,num_gfr=30)
})
B1 <- B1[test,]
B0 <- B0[test,]
xmat_test <- cbind(rep(c,n),w)[test,]
time.barddt.pred <- system.time({
  pred1 <- predict(barddt.fit,xmat_test,B1)$y_hat
  pred0 <- predict(barddt.fit,xmat_test,B0)$y_hat
})
barddt.pred <- pred1 - pred0
## T-BART----------------------------------------------------------------------
### If using stochtree prior to separate param lists
# time.bart.fit <- system.time({
#   bart.fit = bart(X_train= as.matrix(cbind(x,w)), y_train=y,
#                   params=parmlist,
#                   num_mcmc=1000,num_gfr=30)
# })
### If using stochtree with separate param lists
time.tbart.fit <- system.time({
  tbart.fit.0 = bart(X_train= as.matrix(cbind(x,w)[z==0,]), y_train=y[z==0],
                     mean_forest_params=mean.parmlist,
                     general_params=global.parmlist,num_mcmc=1000,num_gfr=30)
  tbart.fit.1 = bart(X_train= as.matrix(cbind(x,w)[z==1,]), y_train=y[z==1],
                     mean_forest_params=mean.parmlist,
                     general_params=global.parmlist,num_mcmc=1000,num_gfr=30)
})
xmat_test <- cbind(c,w)[test,]
time.tbart.pred <- system.time({
  pred1 <- predict(tbart.fit.1,xmat_test)$y_hat
  pred0 <- predict(tbart.fit.0,xmat_test)$y_hat
})
tbart.pred <- pred1 - pred0
## S-BART----------------------------------------------------------------------
### If using stochtree prior to separate param lists
# time.bart.fit <- system.time({
#   bart.fit = bart(X_train= as.matrix(cbind(x,w)), y_train=y,
#                   params=parmlist,
#                   num_mcmc=1000,num_gfr=30)
# })
### If using stochtree with separate param lists
time.sbart.fit <- system.time({
  sbart.fit = bart(X_train= as.matrix(cbind(x,z,w)), y_train=y,
                   mean_forest_params=mean.parmlist,
                   general_params=global.parmlist,num_mcmc=1000,num_gfr=30)
})
xmat_test.1 <- cbind(c,1,w)[test,]
xmat_test.0 <- cbind(c,0,w)[test,]
time.sbart.pred <- system.time({
  pred1 <- predict(sbart.fit,xmat_test.1)$y_hat
  pred0 <- predict(sbart.fit,xmat_test.0)$y_hat
})
sbart.pred <- pred1 - pred0
## Polynomial------------------------------------------------------------------
poly.data <- data.frame(y=y,x=x,w1=w[,1],w2=w[,2],z=z)
bw <- rdrobust::rdbwselect(y,x,c,covs=w)$bws[2]
poly.data <- subset(poly.data, c-bw <= x & x <= c+bw)
time.polynomial.fit <- system.time({
  poly.fit <- lm(y~(poly(x,degree=2,raw=T)*(poly(w1,degree=2,raw=T)*poly(w2,degree=2,raw=T)))+(poly(w1,degree=3,raw=T)*poly(w2,degree=3,raw=T))*z,data=poly.data)
})
xmat_test.1 <- data.frame(x=c,w1=w[test,1],w2=w[test,2],z=1)
xmat_test.0 <- data.frame(x=c,w1=w[test,1],w2=w[test,2],z=0)
time.polynomial.pred <- system.time({
  pred1 <- predict(poly.fit,xmat_test.1)
  pred0 <- predict(poly.fit,xmat_test.0)
})
poly.pred <- pred1-pred0
## Results
print(mget(ls(pattern="time")))
cate <- cate[test]
### Plots
matplot(rowMeans(w)[test],cbind(cate,rowMeans(barddt.pred),rowMeans(tbart.pred),rowMeans(sbart.pred),poly.pred),type="n")
ylim <- par("usr")[3:4]
par(mfrow=c(2,2))
matplot(rowMeans(w)[test],
        cbind(cate,rowMeans(barddt.pred)),
        pch=20,col=c("black","orange"),
        xlab="W",ylab="CATE",bty="n",main="Leaf-level RDD",ylim = ylim)
matplot(rowMeans(w)[test],
        cbind(cate,rowMeans(tbart.pred)),
        pch=20,col=c("black","maroon"),
        xlab="W",ylab="CATE",bty="n",main="T-BART",ylim = ylim)
matplot(rowMeans(w)[test],
        cbind(cate,rowMeans(sbart.pred)),
        pch=20,col=c("black","red"),
        xlab="W",ylab="CATE",bty="n",main="S-BART",ylim = ylim)
matplot(rowMeans(w)[test],
        cbind(cate,poly.pred),
        pch=20,col=c("black","salmon"),
        xlab="W",ylab="CATE",bty="n",main="Polynomial",ylim = ylim)
####
par(mfrow=c(2,2))
matplot(rowMeans(w)[test],
        cbind(cate,rowMeans(barddt.pred)),
        pch=20,col=c("black","orange"),
        xlab="W",ylab="CATE",bty="n",main="Leaf-level RDD")
matplot(rowMeans(w)[test],
        cbind(cate,rowMeans(tbart.pred)),
        pch=20,col=c("black","maroon"),
        xlab="W",ylab="CATE",bty="n",main="T-BART")
matplot(rowMeans(w)[test],
        cbind(cate,rowMeans(sbart.pred)),
        pch=20,col=c("black","red"),
        xlab="W",ylab="CATE",bty="n",main="S-BART")
matplot(rowMeans(w)[test],
        cbind(cate,poly.pred),
        pch=20,col=c("black","salmon"),
        xlab="W",ylab="CATE",bty="n",main="Polynomial")
#####
par(mfrow=c(1,1))
matplot(rowMeans(w)[test],
        cbind(cate,rowMeans(barddt.pred),rowMeans(tbart.pred)),
        pch=20,col=c("black","orange","maroon"),
        xlab="W",ylab="CATE",bty="n")