## Setup
set.seed(0)
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
library(rdrobust)
## Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
##
s1            <- 11
s2            <- 50
n             <- 500
c             <- 0
Omin          <- 5
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 2
##
### DGP
## x <- matrix(rnorm(n*s2,-0.2,1),n,s2)
## x <- matrix(rnorm(n*s2,0,1),n,s2)
## x <- matrix(rnorm(n*s2,0.2,1),n,s2)
## x <- matrix(2*rbeta(n*s2,2,4)-1,n,s2)
## p <- apply(x,2,function(i) (i-min(i))/(max(i)-min(i)))
sig <- 0.1
u <- runif(n,0,1)
w <- vector("list",s2)
for (i in 1:s2)
{
    w1 <- rnorm(n,u,0.15)
    w2 <- runif(n,0,0.5)
    w3 <- rbinom(n,4,u)+1
    w4 <- rbinom(n,1,0.7)+0.5
    w[[i]] <- cbind(w1,w2,w3,w4)
}
x <- matrix(rnorm(n*s2,u-0.5,0.5),n,s2)
z <- apply(x,2,function(i) as.numeric(i>=c))
e <- matrix(rnorm(n*s2,0,sig),n,s2)
mu <- function(x,w) 0.1*w[,3]*w[,4]*sin(w[,1]+w[,2]) + 0.2*w[,3]*w[,4] + 0.1*w[,3]*w[,4]*x + 0.1*w[,3]*w[,4]*(x+0.1)^2 + 0.2*w[,3]*w[,4]*x^3
tau <- function(x,w) 0.01*w[,3]*w[,4]*cos(w[,1]+w[,2]) + 0.05*w[,3]*w[,4] - 0.2*w[,3]*w[,4]*x + 0.05*w[,3]*w[,4]*(x+0.1)^2 - 0.1*w[,3]*w[,4]*x^3
y <- matrix(0,n,s2)
for (i in 1:s2) y[,i] <- mu(x[,i],w[[i]]) + tau(x[,i],w[[i]])*z[,i]
## Plotting data and functions
## d1 <- density(x[w[[1]][,4]==1.5,1])
## d0 <- density(x[w[[1]][,4]==0.5,1])
## plot(x[,1],y[,1],xlab="X",ylab="Y",bty="n",col=w[[1]][,4]+0.5,pch=20,ylim=c(min(y[,1],0),max(y[,1],d1$y,d0$y)),xlim=c(min(x[,1],d1$x,d0$x),max(x[,1],d1$x,d0$x)))
## lines(density(x[w[[1]][,4]==1.5,1]),col=1)
## lines(density(x[w[[1]][,4]==0.5,1]),col=2)
## legend("topright",col=1:2,lty=1,legend=c("W=1.5","W=0.5"))
## abline(v=0,lty=2)
plot(x[,1],y[,1],xlab="X",ylab="Y",bty="n",col=2*w[[1]][,3]*w[[1]][,4],pch=20)
###
y <- y + e
plot(x[,1],y[,1],bty="n",pch=20,col=z[,1]+1,xlab="X",ylab="Y")
## ###
## x1 <- sort(x[,1])
## y1 <- y[order(x[,1]),1]
## e1 <- e[order(x[,1]),1]
## w1 <- w[[1]][order(x[,1]),]
## t1 <- tau(x1,w1)
## mp <- cbind(y1-e1,y1-e1-sig,y1-e1+sig)
## ###
## plot(x1[w1[,2]==0.5],t1[w1[,2]==0.5],bty="n",type="l",lty=1,col=1,xlab="X",ylab="Y",ylim=c(min(t1),max(t1)))
## lines(x1[w1[,2]==1.5],t1[w1[,2]==1.5],lty=1,col=2)
## legend("topright",col=1:2,legend=c("W=0.5","W=1.5"),lty=1,cex=0.75)
## abline(v=0,lty=2)
## ###
## par(mfrow=c(1,2))
## plot(x1,y1,bty="n",col=w1[,2]+0.5,pch=20,xlab="X",ylab="Y")
## abline(v=0,lty=2)
## legend("topleft",col=1:2,legend=c("W=0.5","W=1.5"),lty=1,cex=0.75)
## matplot(x1[w1[,2]==0.5],mp[w1[,2]==0.5,],bty="n",type="l",lty=1,lwd=c(2,1,1),col=1,xlab="X",ylab="Y",ylim=c(min(mp),max(mp)))
## matlines(x1[w1[,2]==1.5],mp[w1[,2]==1.5,],lty=1,lwd=c(2,1,1),col=2)
## legend("topleft",col=1:2,legend=c("W=0.5","W=1.5"),lty=1,cex=0.75)
## abline(v=0,lty=2)
## Prior predictive
pars <- runif(6,0,0.3)
mu.prior <- function(x,w) as.vector(cbind(1,x,x^2+0.1)%*%pars[1:3])
tau.prior <- function(x,w) as.vector(cbind(1,x,x^2+0.1)%*%pars[4:6])
ys <- vector("list",s2)
for (i in 1:s2) ys[[i]] <- mu.prior(x[,i],w[[i]]) + tau.prior(x[,i],w[[i]])*z[,i] + matrix(rnorm(n*s1,0,0.1),n,s1)
plot(x[,1],ys[[1]][,1],bty="n",pch=20,col=z[,1]+1,xlab="X",ylab=expression(Y[prior]))
## Functions
opt.h.int <- function(s,h,ys,x,w,z)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            y <- ys[,i]
            print(paste0("Sample: ",i))
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
            test <- -h+c<=x & x<=h+c
            pred <- predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
opt.h <- function(s,y,x,w,z,h)
{
    out <- vector("list",length(h))
    for (i in 1:length(h))
    {
        print(paste0("h: ",i))
        temp <- opt.h.int(s,h[i],y,x,w,z)
        out[[i]] <- temp
    }
    return(out)
}
h <- round(rbind(0.5*apply(x,2,sd),apply(x,2,sd),2*apply(x,2,sd)),2)
h1 <- vector("list",s2)
bart.prior.time <- system.time({
    for (j in 1:s2) h1[[j]] <- opt.h(s1,ys[[j]],x[,j],w[[j]],z[,j],h[,j])
})
saveRDS(h1,"Results/prior_h.rds")
h1 <- readRDS("Results/prior_h.rds")
out <- matrix(0,s2,nrow(h))
for (i in 1:s2)
{
    for (j in 1:nrow(h))
    {
        out[i,j] <- sqrt(mean(sapply(h1[[i]][[j]],function(x) (mean(colMeans(x))-mean(tau.prior(0,w[[i]])))^2)))
    }
}
h.opt <- mapply(function(i,j) h[j,i],i=1:s2,j=apply(out,1,function(x) which(x==min(x))))
##
fit <- function(s,h,y,x,w,z)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Sample: ",i))
            ys <- y[,i]
            hs <- h[i]
            ws <- w[[i]]
            xs <- x[,i]
            fit <- XBCF.rd(ys, ws, xs, c,
                           Owidth = hs, Omin = Omin, Opct = Opct,
                           num_trees_mod = ntrees,
                           num_trees_con = ntrees,
                           num_cutpoints = n,
                           num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical,
                           p_categorical_mod = p_categorical,
                           tau_con = 2*var(ys)/ntrees,
                           tau_mod = 0.5*var(ys)/ntrees)
            test <- -hs+c<=xs & xs<=hs+c
            pred <- predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
##
ate <- sapply(1:s2,function(i) mean(tau(0,w[[i]])))
bart.time <- system.time({
    bart.fit <- fit(s2,h.opt,y,x,w,z)
})
saveRDS(bart.fit,"Results/prior_test_bart_rdd.rds")
bart.fit <- readRDS("Results/prior_test_bart_rdd.rds")
ate.bart <- colMeans(sapply(bart.fit,colMeans))
rmse.bart <- sqrt(mean((ate.bart-ate)^2))
bias.bart <- mean(ate.bart-ate)
var.bart <- var(ate.bart)
ci.bart <- apply(sapply(bart.fit,colMeans),2,quantile,c(0.025,0.975))
ci.bart <- t(ci.bart)
### CATE
ind <- sapply(1:s2, function(i) -h.opt[i] <= x[,i] & x[,i] <= h.opt[i])
rmse.cate <- mean(sapply(1:s2,function(i) sqrt(mean((rowMeans(bart.fit[[i]])-tau(0,w[[i]])[ind[,i]])^2))))
ci.cate <- lapply(1:s2,function(i) t(apply(bart.fit[[i]],1,quantile,c(0.025,0.975))))
cov.cate <- mean(sapply(1:s2, function(i) mean(ci.cate[[i]][,1] <= tau(0,w[[i]])[ind[,i]] & tau(0,w[[i]])[ind[,i]] <= ci.cate[[i]][,2])))
###
ckt.fit <- function(s,y,x,w)
{
    foreach(i=1:s,.multicombine=T,.export="c") %dopar%
        {
            rdrobust(y[,i],x[,i],c=c,covs=w[[i]])
        }
}
ckt.time <- system.time({
    ckt <- ckt.fit(s2,y,x,w)
})
saveRDS(ckt,"Results/prior_test_polynomial.rds")
ckt <- readRDS("Results/prior_test_polynomial.rds")
ate.ckt <- sapply(ckt, function(i) i$coef[3])
rmse.ckt <- sqrt(mean((ate.ckt-ate)^2))
bias.ckt <- mean(ate.ckt-ate)
var.ckt <- var(ate.ckt)
ci.ckt <- sapply(ckt, function(i) i$ci[3,])
ci.ckt <- t(ci.ckt)
##
print(paste0("RMSE BART: ",rmse.bart))
print(paste0("RMSE Polynomial: ",rmse.ckt))
print(paste0("Bias BART: ",bias.bart))
print(paste0("Bias Polynomial: ",bias.ckt))
print(paste0("Variance BART: ",var.bart))
print(paste0("Variance Polynomial: ",var.ckt))
print(paste0("Coverage BART: ", mean(ci.bart[,1] <= ate & ate <= ci.bart[,2])))
print(paste0("Coverage Polynomial: ", mean(ci.ckt[,1] <= ate & ate <= ci.ckt[,2])))
print(paste0("Interval Size BART: ", mean(- ci.bart[,1] + ci.bart[,2])))
print(paste0("Interval Size Polynomial: ", mean(- ci.ckt[,1] + ci.ckt[,2])))
print("Distance to tau")
print(summary(abs(rowSums(cbind(ifelse(ate-ci.bart[,1]>0,0,ate-ci.bart[,1]),ifelse(ci.bart[,2]-ate>0,0,ci.bart[,2]-ate))))))
print("95-th percentile of distance to tau")
print(quantile(abs(rowSums(cbind(ifelse(ate-ci.bart[,1]>0,0,ate-ci.bart[,1]),ifelse(ci.bart[,2]-ate>0,0,ci.bart[,2]-ate)))),0.95))
print("Optimal h")
print(table(h.opt))
print(paste0("sd(X) = ",sd(x)))
print("Summary ATE:")
print(summary(ate))
print("Summary ATE BART:")
print(summary(ate.bart))
print("Summary ATE Polynomial:")
print(summary(ate.ckt))
print("CATE summary")
print("RMSE:")
print(rmse.cate)
print("Coverage:")
print(cov.cate)
print(paste0("BART prior tuning time: ",bart.prior.time))
print(paste0("BART time: ",bart.time))
print(paste0("Polynomial time: ",ckt.time))
print("Done!")