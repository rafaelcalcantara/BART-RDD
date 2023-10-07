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
s1            <- 20
s2            <- 500
n             <- 500
c             <- 0
Omin          <- 5
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 1
##
### DGP
sig <- 0.15
x <- matrix(0.6*rnorm(n*s2,-1,0.5) + 0.4*rnorm(n*s2,1,0.5),n,s2)
w <- matrix(rbinom(n*s2,3,0.7)+1,n,s2)
z <- apply(x,2,function(i) as.numeric(i>=c))
e <- matrix(rnorm(n*s2,0,sig),n,s2)
mu <- function(x,w) 0.1*w + 0.05*w*x - 0.2*w*x^2 + 0.1*(w-1)*x^3
tau <- function(x,w) 0.01*(w+1) - 0.02*w*x - 0.05*w*x^2 + 0.5*w*x^3
y <- mu(x,w) + tau(x,w)*z + e
plot(x[,1],y[,1],bty="n",pch=20,col=z+1,xlab="X",ylab="Y")
## Prior predictive
mu.prior <- function(x,w) (1+exp(1))/(1+exp(-x)) - 2 + 0.2*w
tau.prior <- function(x) 0.25*x - 0.08*x^3 + 0.05
ys <- vector("list",s2)
for (i in 1:s2) ys[[i]] <- mu.prior(x[,i],w[,i]) + tau.prior(x[,i])*z[,i] + matrix(rnorm(n*s1,0,0.1),n,s1)
plot(x[,1],ys[[1]][,1],bty="n",pch=20,col=z+1,xlab="X",ylab=expression(Y[prior]))
## w <- vector("list",s2)
## for (i in 1:s2) w[[i]] <- cbind(rbinom(n,1,0.8),runif(n,-0.5,0.5))
## z <- apply(x,2,function(i) as.numeric(i>=c))
## e <- matrix(rnorm(n*s2,0,sig),n,s2)
## mu <- function(x,w,i)
## {
##     if (i==1) 0.4 + 0.5*x - 0.2*x^2 + 0.3*x^3 + 0.02*w[,2]
##     else pnorm(x,0,0.5) + 0.5*sin(x^2) + 0.1*w[,2]
## }
## tau <- function(x,w,i)
## {
##     if (i==1) 0.05
##     else if (i==2) 0.05 - 0.1*x + 0.05*x^2 + 0.05*w[,1] - 0.02*(1-w[,1])
##     else if (i==3) 0.05 - 0.1*x + 0.05*x^2 + 0.02*exp(-w[,2])
##     else 0.05 - 0.1*x + 0.5*x^2 + 0.05*w[,1] - 0.02*(1-w[,1]) + (0.05*w[,1] + 0.02*(1-w[,1]))*exp(-w[,2])
## }
## y <- vector("list",s2)
## for (i in 1:s2)
## {
##     index <- 0
##     for (j in 1:2)
##     {
##         for (k in 1:4)
##         {
##             index <- index+1
##             temp <- mu(x[,i],w[[i]],j) + tau(x[,i],w[[i]],k)*z[,i] + e[,i]
##             y[[i]] <- cbind(y[[i]],temp)
##         }
##     }
## }
## ## Sample y for prior predictive
## ate.prior <- function(x,w) 0.05 + 0.05*w[,1] - 0.1*w[,2]
## ys <- vector("list",s2)
## for (j in 1:s2)
## {
##     temp <- matrix(0,n,s1)
##     for (i in 1:s1) temp[,i] <- ate.prior(x[,j],w[[j]])*z[,j] + 0.1*w[[j]][,2] + 1/(1+exp(-x[,j])) + rnorm(n,0,0.1)
##     ys[[j]] <- temp
## }
## ##
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
            pred <- predict.XBCFrd(fit,w[test],rep(c,sum(test)))
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
h <- apply(x,2,function(i) unique(sort(abs(round(quantile(i,seq(0,1,0.025)),2))))[3:7])
h1 <- vector("list",s2)
for (j in 1:s2) h1[[j]] <- opt.h(s1,ys[[j]],x[,j],w[,j],z[,j],h[,j])
saveRDS(h1,"Results/prior_h.rds")
out <- matrix(0,s2,nrow(h))
for (i in 1:s2)
{
    for (j in 1:nrow(h))
    {
        out[i,j] <- sqrt(mean(sapply(h1[[i]][[j]],function(x) (mean(colMeans(x))-mean(tau.prior(x[,i],w[,i])))^2)))
    }
}
h.opt <- mapply(function(i,j) h[j,i],i=1:s1,j=apply(out,1,function(x) which(x==min(x))))
##
fit <- function(s,h,y,x,w,z)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("Sample: ",i))
            ys <- y[,i]
            hs <- h.opt[i]
            ws <- w[,i]
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
            pred <- predict.XBCFrd(fit,ws[test],rep(c,sum(test)))
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
##
bart.fit <- fit(s2,h.opt,y,x,w,z)
saveRDS(bart.fit,"Results/prior_test_bart_rdd.rds")
## ate.hat <- sapply(lapply(bart.fit, function(x) sapply(x,colMeans)),colMeans)
## rmse <- ate.hat
## for (i in 1:4) rmse[,c(i,i+4)] <- (rmse[,c(i,i+4)]-sapply(1:s2,function(x) mean(tau(0,w[[x]],i))))^2
## rmse <- colMeans(rmse)
## rmse <- sqrt(rmse)
###
ckt.fit <- function(s,y,x,w)
{
    foreach(i=1:s,.multicombine=T,.export="c") %dopar%
        {
            rdrobust(y[,i],x[,i],c=c,covs=w[,i])
        }
}
ckt <- ckt.fit(s2,y,x,w)
saveRDS(ckt,"Results/prior_test_polynomial.rds")
## ate.ckt <- sapply(ckt, function(i) sapply(i,function(j) j$coef[3]))
## rmse.ckt <- ate.ckt
## for (i in 1:4) rmse.ckt[,c(i,i+4)] <- (rmse.ckt[,c(i,i+4)]-sapply(1:s2,function(x) mean(tau(0,w[[x]],i))))^2
## rmse.ckt <- colMeans(rmse.ckt)
## rmse.ckt <- sqrt(rmse.ckt)
## ###
## prior.point.estimates <- t(sapply(lapply(h1,function(x) sapply(x,colMeans)),colMeans))
## matplot(x=h,(prior.point.estimates-0.05-0.05*mean(w[,1])+0.1*mean(w[,2]))^2,type="b",lty=2,pch=20,col="blue",ylab=expression(tau),bty="n")
## h2 <- opt.h(ncol(y[[1]]),y[[1]],x,w,z,h)
## post.point.estimates <- t(sapply(lapply(h2,function(x) sapply(x,colMeans)),colMeans))
## for (i in 1:4) post.point.estimates[,c(i,i+4)] <- (post.point.estimates[,c(i,i+4)] - mean(tau(0,w,i)))^2
## matplot(x=h,post.point.estimates,type="b",lty=2,pch=20,col="blue",ylab=expression(tau),bty="n")
## ## h2 <- opt.h(s,ate2,y2,x,w,z,h)
## ## h3 <- opt.h(s,ate3,y3,x,w,z,h)
## ## saveRDS(list(h1,h2,h3),"Results/prior.rds")
## ## ###
## ## prior <- readRDS("Results/prior.rds")
## ## h1 <- prior[[1]]
## ## h2 <- prior[[2]]
## ## h3 <- prior[[3]]
## ## ###
## rmse1 <- sapply(h1,function(x) sapply(x,function(y) abs(mean(colMeans(y))-ate1)/ate1))
## ## rmse2 <- sapply(h2,function(x) sapply(x,function(y) abs(mean(colMeans(y))-ate2)/ate2))
## ## rmse3 <- sapply(h3,function(x) sapply(x,function(y) abs(mean(colMeans(y))-ate3)/ate3))
## rmse1 <- as.data.frame(rmse1)
## ## rmse2 <- as.data.frame(rmse2)
## ## rmse3 <- as.data.frame(rmse3)
## ## names(rmse1) <- names(rmse2) <- names(rmse3) <- h
## names(rmse1) <- h
## ate.hat <- sapply(h1,function(x) sapply(x, function(y) mean(colMeans(y))))
## ate.hat <- data.frame(ate.hat)
## names(ate.hat) <- h
## ###
## cate.hat1 <- matrix(0,s,length(h))
## for (i in 1:length(h))
## {
##     cate.hat1[,i] <- sapply(h1[[i]],function(a) mean(colMeans(a[w[x<=h[i]+c & x>=-h[i]+c]==1,])))
## }
## cate.hat1 <- data.frame(cate.hat1)
## names(cate.hat1) <- h
## cate.hat0 <- matrix(0,s,length(h))
## for (i in 1:length(h))
## {
##     cate.hat0[,i] <- sapply(h1[[i]],function(a) mean(colMeans(a[w[x<=h[i]+c & x>=-h[i]+c]==0,])))
## }
## cate.hat0 <- data.frame(cate.hat0)
## names(cate.hat0) <- h
## ## par(mfrow=c(1,3))
## barplot(table(h[apply(rmse1,1,function(x) which(x==min(x)))]))
## ## barplot(table(h[apply(rmse2,1,function(x) which(x==min(x)))]))
## ## barplot(table(h[apply(rmse3,1,function(x) which(x==min(x)))]))
## ###
## ## par(mfrow=c(1,3))
## matplot(x=h,y=t(rmse1),pch=19,col="blue",type="b",ylab="Rel. error",main=paste(expression(tau),"=3"),lty=2)
## ## matplot(x=h,y=t(rmse2),pch=19,col="blue",type="b",ylab="Rel. error",main=paste(expression(tau),"=0.2"),lty=2)
## ## matplot(x=h,y=t(rmse3),pch=19,col="blue",type="b",ylab="Rel. error",main=paste(expression(tau),"=0.5"),lty=2)
## ## ###
## boxplot(ate.hat)
## abline(h=ate1,lty=2)
## ## ###
## par(mfrow=c(1,2))
## boxplot(cate.hat1,main="w=1")
## abline(h=6,lty=2)
## boxplot(cate.hat0,main="w=0")
## abline(h=3,lty=2)
## ###
## par(mfrow=c(3,1))
## barplot(table(h1$h[apply(h1$rmse,1,function(x) which(x==min(x)))]))
## barplot(table(h1$h[apply(h2$rmse,1,function(x) which(x==min(x)))]))
## barplot(table(h1$h[apply(h3$rmse,1,function(x) which(x==min(x)))]))
## ###
## mse <- data.frame(h1a$rmse)
## names(mse) <- h1a$h
## h.plot <- t(sapply(h1b,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
## par(mfrow=c(1,3))
## boxplot(sqrt(mse),main="RMSE")
## plot(density(rowMeans(h1a$ate)),main="ATE")
## abline(v=mean(tau.fun(w,0)),lty=2)
## matplot(x=h1a$h,h.plot,col="blue",lty=c(1,2,2),type=c("b","l","l"),pch=21,main="Fit for true y")
## abline(h=mean(tau.fun(w,0)),lty=2)
## ###
## pdf("Figures/prior1.pdf")
## par(mfrow=c(2,2))
## plot(h1a$h,sqrt(colMeans(h1a$rmse)),bg=ifelse(colMeans(h1a$rmse)==min(colMeans(h1a$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Negative")
## plot(h1b$h,sqrt(colMeans(h1b$rmse)),bg=ifelse(colMeans(h1b$rmse)==min(colMeans(h1b$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Positive")
## plot(h1c$h,sqrt(colMeans(h1c$rmse)),bg=ifelse(colMeans(h1c$rmse)==min(colMeans(h1c$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Zero")
## dev.off()
## ## Data 2
## x <- 2*rbeta(n,2,4)-0.5
## w <- matrix(rnorm(2*n,0,0.25),n,2)
## z <- x>=0
## mu.fun <- function(W, X){return(0.1*rowSums(W) + 1/(1+exp(-5*X)))} 
## tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
## y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n)
## h2a <- opt.h(s,x,w,z,1)
## h2b <- opt.h(s,x,w,z,2)
## h2c <- opt.h(s,x,w,z,3)
## h2d <- unlist(fit(h2a$h))
## h2d <- abs(h2d-mean(tau.fun(w,0)))
## ###
## ## pdf("Figures/prior2.pdf")
## par(mfrow=c(2,2))
## plot(h2a$h,sqrt(colMeans(h2a$rmse)),bg=ifelse(colMeans(h2a$rmse)==min(colMeans(h2a$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Negative")
## plot(h2b$h,sqrt(colMeans(h2b$rmse)),bg=ifelse(colMeans(h2b$rmse)==min(colMeans(h2b$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Positive")
## plot(h2c$h,sqrt(colMeans(h2c$rmse)),bg=ifelse(colMeans(h2c$rmse)==min(colMeans(h2c$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Zero")
## plot(h2c$h,h2d,bg=ifelse(h2d==min(h2d),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Real data")
## ## dev.off()
## ## Data 3
## x <- runif(n,-1,1)
## w <- matrix(rnorm(2*n,0,0.25),n,2)
## z <- x>=0
## mu.fun <- function(W, X){return(0.1*rowSums(W) + 1/(1+exp(-5*X)))} 
## tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
## y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n)
## h3a <- opt.h(s,x,w,z,1)
## h3b <- opt.h(s,x,w,z,2)
## h3c <- opt.h(s,x,w,z,3)
## h3d <- unlist(fit(h3a$h))
## h3d <- abs(h3d-mean(tau.fun(w,0)))
## ###
## ## pdf("Figures/prior3.pdf")
## par(mfrow=c(2,2))
## plot(h3a$h,sqrt(colMeans(h3a$rmse)),bg=ifelse(colMeans(h3a$rmse)==min(colMeans(h3a$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Negative")
## plot(h3b$h,sqrt(colMeans(h3b$rmse)),bg=ifelse(colMeans(h3b$rmse)==min(colMeans(h3b$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Positive")
## plot(h3c$h,sqrt(colMeans(h3c$rmse)),bg=ifelse(colMeans(h3c$rmse)==min(colMeans(h3c$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Zero")
## plot(h3c$h,h3d,bg=ifelse(h3d==min(h3d),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Real data")
## ## dev.off()
## ## Data 4
## x <- rnorm(n,0,0.25)
## w <- matrix(rnorm(2*n,0,0.25),n,2)
## z <- x>=0
## mu.fun <- function(W, X){return(0.1*rowSums(W) + 1/(1+exp(-5*X)))} 
## tau.fun <- function(W, X) return( sin(mu.fun(W, X)) ) # make sure the treatment effect is non-zero
## y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n)
## h4a <- opt.h(s,x,w,z,1)
## h4b <- opt.h(s,x,w,z,2)
## h4c <- opt.h(s,x,w,z,3)
## h4d <- unlist(fit(h4a$h))
## h4d <- abs(h4d-mean(tau.fun(w,0)))
## ###
## ## pdf("Figures/prior4.pdf")
## par(mfrow=c(2,2))
## plot(h4a$h,sqrt(colMeans(h1a$rmse)),bg=ifelse(colMeans(h4a$rmse)==min(colMeans(h4a$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Negative")
## plot(h4b$h,sqrt(colMeans(h4b$rmse)),bg=ifelse(colMeans(h4b$rmse)==min(colMeans(h4b$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Positive")
## plot(h4c$h,sqrt(colMeans(h4c$rmse)),bg=ifelse(colMeans(h4c$rmse)==min(colMeans(h4c$rmse)),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Zero")
## plot(h4c$h,h4d,bg=ifelse(h4d==min(h4d),"blue","red"),bty="n",xlab="h",ylab="RMSE",pch=21,main="Real data")
## ## dev.off()
## ####
print("Done!")
