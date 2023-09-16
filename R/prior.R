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
s             <- 10
n             <- 750
c             <- 0
Omin          <- 5
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 100
burnin        <- 20
p_categorical <- 0
##
### Prior model
ysamp <- function(x,w,U,tau)
{
    theta1 <- rep(0.5,5)
    theta2 <- c(tau,tau,tau)
    muX <- cbind(1,scale(x),scale(x^2),scale(sin(w)),U)
    tauX <- cbind(1,scale(sqrt(U)),scale(w))
    ys <- vector("list",s)
    for (i in 1:s) ys[[i]] <- muX%*%theta1 + z*tauX%*%theta2 + rnorm(n,0,0.5)
    return(ys)
}
ate <- function(tau,w,U) tau*(1 + mean(scale(sqrt(U))) + mean(scale(w)))
opt.h.int <- function(s,h,ate,ys,x,w,z)
{
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            y <- ys[[i]]
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
                           tau_mod = 0.5*var(y)/ntrees,
                           random_seed=0)
            test <- -h<=x & x<=h
            pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
            pred$tau.adj[,(burnin+1):num_sweeps]
        }
}
opt.h <- function(s,ate,y,x,w,z,h)
{
    out <- vector("list",length(h))
    for (i in 1:length(h))
    {
        print(paste0("h: ",i))
        temp <- opt.h.int(s,h[i],ate,y,x,w,z)
        out[[i]] <- temp
    }
    return(out)
}
## Data 1
U <- runif(n,0,1)
x <- U + 2*rbeta(n,2,4)-1.5
w <- rnorm(n,0,0.25)
z <- x>=0
y1 <- ysamp(x,w,U,0.05)
y2 <- ysamp(x,w,U,0.2)
y3 <- ysamp(x,w,U,0.5)
ate1 <- ate(0.05,w,U)
ate2 <- ate(0.2,w,U)
ate3 <- ate(0.5,w,U)
###
h <- sort(abs(quantile(x,seq(0,1,0.025))))
h <- round(h,2)
h <- unique(h)
h <- h[1:10]
h <- seq(0.05,0.2,0.025)
h1 <- opt.h(s,ate1,y1,x,w,z,h)
h2 <- opt.h(s,ate2,y2,x,w,z,h)
h3 <- opt.h(s,ate3,y3,x,w,z,h)
saveRDS(list(h1,h2,h3),"Results/prior.rds")
## ###
## prior <- readRDS("Results/prior.rds")
## h1 <- prior[[1]]
## h2 <- prior[[2]]
## h3 <- prior[[3]]
## ###
rmse1 <- sapply(h1,function(x) sapply(x,function(y) abs(mean(colMeans(y))-ate1)/ate1))
rmse2 <- sapply(h2,function(x) sapply(x,function(y) abs(mean(colMeans(y))-ate2)/ate2))
rmse3 <- sapply(h3,function(x) sapply(x,function(y) abs(mean(colMeans(y))-ate3)/ate3))
rmse1 <- as.data.frame(rmse1)
rmse2 <- as.data.frame(rmse2)
rmse3 <- as.data.frame(rmse3)
names(rmse1) <- names(rmse2) <- names(rmse3) <- h
###
par(mfrow=c(1,3))
barplot(table(h[apply(rmse1,1,function(x) which(x==min(x)))]))
barplot(table(h[apply(rmse2,1,function(x) which(x==min(x)))]))
barplot(table(h[apply(rmse3,1,function(x) which(x==min(x)))]))
###
par(mfrow=c(1,3))
matplot(x=h,y=t(rmse1),pch=19,col="blue",type="b",ylab="Rel. error",main=paste(expression(tau),"=0.05"),lty=2)
matplot(x=h,y=t(rmse2),pch=19,col="blue",type="b",ylab="Rel. error",main=paste(expression(tau),"=0.2"),lty=2)
matplot(x=h,y=t(rmse3),pch=19,col="blue",type="b",ylab="Rel. error",main=paste(expression(tau),"=0.5"),lty=2)
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
## print("Done!")
