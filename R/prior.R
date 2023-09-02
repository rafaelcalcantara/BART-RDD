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
k             <- 0.25
n             <- 500
c             <- 0
Omin          <- 1
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 100
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
                           tau_mod = 0.5*var(y)/ntrees,
                           random_seed=0)
            test <- -h<=x & x<=h
            print(paste0("i: ",i,"; test: ",sum(test)))
            if (sum(test)>1)
            {
                pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
                pred <- pred$tau.adj[,(burnin+1):num_sweeps]
                pred <- colMeans(pred)
                c(mean(pred),var(pred),k,h,
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
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + k*z + rnorm(n)
s <- quantile(x,seq(0.3,0.9,0.05))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau1 <- fit(s)
tau1 <- do.call("rbind",tau1)
h1 <- tau1[which(tau1[,2]==max(tau1[,2],na.rm=T)),4]
opt.pred1 <- opt.fit(h1)
opt.pred1 <- predict(loess(y~x,data=data.frame(y=opt.pred1,x=x)))
true.tau1 <- k+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
true.tau1 <- predict(loess(y~x,data=data.frame(y=true.tau1,x=x)))
### Plot tau hat
xaxis <- quantile(x,seq(0.3,0.9,0.05))
x1 <- sort(x)
y1 <- cbind(opt.pred1,true.tau1)[order(x),]
## Data 2
x <- 2*rbeta(n,2,4)-0.5
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*W + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + k*z + rnorm(n)
s <- quantile(x,seq(0.3,0.9,0.05))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau2 <- fit(s)
tau2 <- do.call("rbind",tau2)
h2 <- tau2[which(tau2[,2]==max(tau2[,2],na.rm=T)),4]
opt.pred2 <- opt.fit(h2)
opt.pred2 <- predict(loess(y~x,data=data.frame(y=opt.pred2,x=x)))
true.tau2 <- k+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
true.tau2 <- predict(loess(y~x,data=data.frame(y=true.tau2,x=x)))
### Plot tau hat
xaxis <- cbind(xaxis,quantile(x,seq(0.3,0.9,0.05)))
x2 <- sort(x)
y2 <- cbind(opt.pred2,true.tau2)[order(x),]
## Data 3
x <- runif(n,-1,1)
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*W + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + k*z + rnorm(n)
s <- quantile(x,seq(0.3,0.9,0.05))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau3 <- fit(s)
tau3 <- do.call("rbind",tau3)
h3 <- tau3[which(tau3[,2]==max(tau3[,2],na.rm=T)),4]
opt.pred3 <- opt.fit(h3)
opt.pred3 <- predict(loess(y~x,data=data.frame(y=opt.pred3,x=x)))
true.tau3 <- k+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
true.tau3 <- predict(loess(y~x,data=data.frame(y=true.tau3,x=x)))
### Plot tau hat
xaxis <-cbind(xaxis,quantile(x,seq(0.3,0.9,0.05)))
x3 <- sort(x)
y3 <- cbind(opt.pred3,true.tau3)[order(x),]
## Data 4
x <- rnorm(n,0,0.25)
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*W + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
y <- mu.fun(w, x)/var(mu.fun(w, x)) + (tau.fun(w, x) - mean(tau.fun(w,0)))*z/var(tau.fun(w, x)) + k*z + rnorm(n)
s <- quantile(x,seq(0.3,0.9,0.05))
## index <- which(s[s<0]==max(s[s<0]))
## s <- s[(index-5):(index+5)]
##
tau4 <- fit(s)
tau4 <- do.call("rbind",tau4)
h4 <- tau4[which(tau4[,2]==max(tau4[,2],na.rm=T)),4]
opt.pred4 <- opt.fit(h4)
opt.pred4 <- predict(loess(y~x,data=data.frame(y=opt.pred4,x=x)))
true.tau4 <- k+(tau.fun(w,x)-mean(tau.fun(w,0)))/var(tau.fun(w,x))
true.tau4 <- predict(loess(y~x,data=data.frame(y=true.tau4,x=x)))
### Plot tau hat
x4 <- sort(x)
y4 <- cbind(opt.pred4,true.tau4)[order(x),]
##
## xaxis <- as.numeric(gsub("%","",names(s)))
xaxis <- cbind(xaxis,quantile(x,seq(0.3,0.9,0.05)))
opt.h <- c(h1,h2,h3,h4)
##
par(mfrow=c(2,2))
plot(x=xaxis[,1],y=tau1[,1],col=ifelse(tau1[,2]==max(tau1[,2]),"blue","red"),type="p",pch=19,ylim=c(k-0.5,max(tau1[,1])+0.5),bty="n",xlab="X",ylab=expression(tau),main="Skewed to the left")
abline(h=k,lty=2)
plot(x=xaxis[,2],y=tau2[,1],col=ifelse(tau2[,2]==max(tau2[,2]),"blue","red"),type="p",pch=19,ylim=c(k-0.5,max(tau2[,1])+0.5),bty="n",xlab="X",ylab=expression(tau),main="Skewed to the right")
abline(h=k,lty=2)
plot(x=xaxis[,3],y=tau3[,1],col=ifelse(tau3[,2]==max(tau3[,2]),"blue","red"),type="p",pch=19,ylim=c(k-0.5,max(tau3[,1]+1)),bty="n",xlab="X",ylab=expression(tau),main="Uniform")
abline(h=k,lty=2)
plot(x=xaxis[,4],y=tau4[,1],col=ifelse(tau4[,2]==max(tau4[,2]),"blue","red"),type="p",pch=19,ylim=c(k-0.5,max(tau4[,1]+0.5)),bty="n",xlab="X",ylab=expression(tau),main="Gaussian")
abline(h=k,lty=2)
## ##
## par(mfrow=c(2,2))
## matplot(x=x1,y=y1,type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau),main="Skewed to the left")
## abline(v=0,lty=2)
## matplot(x=x2,y=y2,type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau),main="Skewed to the right")
## abline(v=0,lty=2)
## matplot(x=x3,y=y3,type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau),main="Uniform")
## abline(v=0,lty=2)
## matplot(x=x4,y=y4,type="l",lty=1,col=c("blue","black"),xlab="X",ylab=expression(tau),main="Gaussian")
## abline(v=0,lty=2)
