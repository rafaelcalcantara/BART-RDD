## Setup
set.seed(0)
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
## Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
##
n             <- 1000
c             <- 0
Omin          <- 2
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 150
burnin        <- 50
p_categorical <- 0
##
fit <- function(s)
{
    tau <- rep(0,length(s))
    foreach(i=1:length(s),.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            h <- abs(s[i])
            fit <- XBCF.rd(y, w, x, c,
                           Owidth = h, Omin = Omin, Opct = Opct,
                           num_trees_mod = ntrees,
                           num_trees_con = ntrees,
                           num_cutpoints = n,
                           num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical,
                           p_categorical_mod = p_categorical,
                           tau_con = 2*var(y)/ntrees,
                           tau_mod = 0.5*var(y)/ntrees)
            test <- -h<=x & x<=h
            print(paste0("i: ",i,"; test: ",sum(test)))
            if (sum(test)>1)
            {
                pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
                pred <- pred$tau.adj[,(burnin+1):num_sweeps]
                pred <- colMeans(pred)
                c(mean(pred),var(pred),2,h,
                  mean(pred)-sd(pred),mean(pred)+sd(pred),sum(test))
            } else
            {
                rep(NA,6)
            }
        }
}
opt.fit <- function(h)
{
    fit <- XBCF.rd(y, w, x, c,
                   Owidth = h, Omin = Omin, Opct = Opct,
                   num_trees_mod = ntrees,
                   num_trees_con = ntrees,
                   num_cutpoints = n,
                   num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_categorical,
                   p_categorical_mod = p_categorical,
                   tau_con = 2*var(y)/ntrees,
                   tau_mod = 0.5*var(y)/ntrees)
    pred <- predict.XBCFrd(fit,w,x)
    return(pred$tau.adj.mean)
}
## Data 1
x <- 2*rbeta(n,2,4)-1
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*W + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + 2*z + rnorm(n)
s <- quantile(x,seq(0.1,0.9,0.025))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau1 <- fit(s)
tau1 <- do.call("rbind",tau1)
h1 <- tau1[which(tau1[,2]==max(tau1[,2],na.rm=T)),4]
opt.pred1 <- opt.fit(h1)
true.tau1 <- 2+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
### Plot tau hat
matplot(x=sort(x),y=cbind(opt.pred1,true.tau1)[order(x),],type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau))
abline(v=0,lty=2)
p1 <- recordPlot()
## Data 2
x <- 2*rbeta(n,2,4)-0.5
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*W + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + 2*z + rnorm(n)
s <- quantile(x,seq(0.1,0.9,0.025))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau2 <- fit(s)
tau2 <- do.call("rbind",tau2)
h2 <- tau2[which(tau2[,2]==max(tau2[,2],na.rm=T)),4]
opt.pred2 <- opt.fit(h2)
true.tau2 <- 2+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
### Plot tau hat
matplot(x=sort(x),y=cbind(opt.pred2,true.tau2)[order(x),],type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau))
abline(v=0,lty=2)
p2 <- recordPlot()
## Data 3
x <- runif(n,-1,1)
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*W + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + 2*z + rnorm(n)
s <- quantile(x,seq(0.1,0.9,0.025))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau3 <- fit(s)
tau3 <- do.call("rbind",tau3)
h3 <- tau3[which(tau3[,2]==max(tau3[,2],na.rm=T)),4]
opt.pred3 <- opt.fit(h3)
true.tau3 <- 2+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
### Plot tau hat
matplot(x=sort(x),y=cbind(opt.pred3,true.tau3)[order(x),],type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau))
abline(v=0,lty=2)
p3 <- recordPlot()
## Data 4
x <- rnorm(n,0,0.25)
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*W + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + 2*z + rnorm(n)
s <- quantile(x,seq(0.1,0.9,0.025))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau4 <- fit(s)
tau4 <- do.call("rbind",tau4)
h4 <- tau4[which(tau4[,2]==max(tau4[,2],na.rm=T)),4]
opt.pred4 <- opt.fit(h4)
true.tau4 <- 2+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
### Plot tau hat
matplot(x=sort(x),y=cbind(opt.pred4,true.tau4)[order(x),],type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau))
abline(v=0,lty=2)
p4 <- recordPlot()
##
xaxis <- as.numeric(gsub("%","",names(s)))
opt.h <- c(h1,h2,h3,h4)
##
par(mfrow=c(2,2))
matplot(x=xaxis,y=tau1[,c(1,2,3)],col=c("blue","red","black"),type="l",lty=c(1,1,2,2,2),xlab="Quantile",ylab=expression(tau))
matplot(x=xaxis,y=tau2[,c(1,2,3)],col=c("blue","red","black"),type="l",lty=c(1,1,2,2,2),xlab="Quantile",ylab=expression(tau))
matplot(x=xaxis,y=tau3[,c(1,2,3)],col=c("blue","red","black"),type="l",lty=c(1,1,2,2,2),xlab="Quantile",ylab=expression(tau))
matplot(x=xaxis,y=tau4[,c(1,2,3)],col=c("blue","red","black"),type="l",lty=c(1,1,2,2,2),xlab="Quantile",ylab=expression(tau))
##
par(mfrow=c(2,2))
p1
p2
p3
p4
