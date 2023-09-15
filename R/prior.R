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
n             <- 500
c             <- 0
Omin          <- 5
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 5
num_sweeps    <- 100
burnin        <- 20
p_categorical <- 0
##
opt.h.int <- function(s,h,x,w,z)
{
    ## Prior model
    theta1 <- matrix(rnorm(s*6,0.5,0.1),s,6)
    theta2 <- matrix(rnorm(s*6,0.6,0.1),s,6)
### New data matrix
    Xx <- cbind(1,x,x^2,x^3,x^4,sin(w))
### Generate sample data
    ys <- vector("list",s)
    ate <- theta2[,1]-theta1[,1]+(theta2[,6]-theta1[,6])*mean(sin(w))
    for (i in 1:s)
    {
        ys[[i]] <- (!z)*Xx%*%theta1[i,1:6] + z*Xx%*%theta2[i,1:6] + rnorm(length(x),0,0.5)
    }
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            y <- ys[[i]]
            true.ate <- ate[i]
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
            pred <- pred$tau.adj
            pred <- mean(colMeans(pred))
            error <- (true.ate-pred)^2
            c(Error=error,ATE=ate[i])
        }
}
opt.h <- function(s,x,w,z)
{
    h <- sort(abs(quantile(x,seq(0,1,0.025))))
    h <- round(h,2)
    h <- unique(h)
    h <- h[1:10]
    ## h <- seq(min(h),max(h),length.out=10)
    rmse <- matrix(0,s,length(h))
    ate <- matrix(0,s,length(h))
    for (i in 1:length(h))
    {
        print(paste0("h: ",i))
        temp <- opt.h.int(s,h[i],x,w,z)
        rmse[,i] <- do.call("rbind",temp)[,"Error"]
        ate[,i] <- do.call("rbind",temp)[,"ATE"]
    }
    return(list(rmse=rmse,h=h,ate=ate))
}
fit <- function(h,y,x)
{
    foreach(i=1:length(h),.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            print(paste0("h: ",i))
            fit <- XBCF.rd(y, w, x, c,
                           Owidth = h[i], Omin = Omin, Opct = Opct,
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
            test <- -h[i]<=x & x<=h[i]
            pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
            pred <- pred$tau.adj
            pred <- colMeans(pred)
            pred
        }
}
## Data 1
x <- 2*rbeta(n,2,4)-1
w <- rnorm(n,0,0.25)
z <- x>=0
tau <- function(W, X) return(0.2 + sin(W) + as.numeric(X>=0))
mu1 <- function(W, X) return(X - 0.1*X^2 + 0.5*W)
mu2 <- function(W, X) return(0.1*W + exp(W) + 1/(1+exp(-5*X)))
mu3 <- function(W, X) return((X<0.1)*cos(W) + abs(W) + X^2 + exp(X))
y1 <- mu1(w, x) + tau(w, x)*z1 + rnorm(n,0,0.5)
y2 <- mu2(w, x) + tau(w, x)*z2 + rnorm(n,0,0.5)
y3 <- mu3(w, x) + tau(w, x)*z3 + rnorm(n,0,0.5)
###
h1 <- opt.h(s,x,w,z)
fit1 <- fit(h1$h,y1,x)
h2 <- opt.h(s,x,w,z)
fit2 <- fit(h2$h,y2,x)
h3 <- opt.h(s,x,w,z)
fit3 <- fit(h3$h,y3,x)
###
h.plot1 <- t(sapply(fit1,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
h.plot2 <- t(sapply(fit2,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
h.plot3 <- t(sapply(fit3,function(x) c(mean(x),quantile(x,c(0.025,0.975),na.rm=T))))
par(mfrow=c(3,2))
hist(h1$h[apply(h1$rmse,1,function(x) which(x==min(x)))],main="h",xlab="")
matplot(x=h1$h,h.plot1,col="blue",bg="blue",lty=c(1,2,2),type=c("b","l","l"),pch=21,main="")
abline(h=mean(tau(w,0)),lty=2)
hist(h2$h[apply(h2$rmse,1,function(x) which(x==min(x)))],main="h",xlab="")
matplot(x=h2$h,h.plot2,col="blue",bg="blue",lty=c(1,2,2),type=c("b","l","l"),pch=21,main="")
abline(h=mean(tau(w,0)),lty=2)
hist(h3$h[apply(h3$rmse,1,function(x) which(x==min(x,na.rm=T)))],main="h",xlab="")
matplot(x=h3$h,h.plot3,col="blue",bg="blue",lty=c(1,2,2),type=c("b","l","l"),pch=21,main="")
abline(h=mean(tau(w,0)),lty=2)
###
mse <- data.frame(h1a$rmse)
names(mse) <- h1a$h
h.plot <- t(sapply(h1b,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
par(mfrow=c(1,3))
boxplot(sqrt(mse),main="RMSE")
plot(density(rowMeans(h1a$ate)),main="ATE")
abline(v=mean(tau.fun(w,0)),lty=2)
matplot(x=h1a$h,h.plot,col="blue",lty=c(1,2,2),type=c("b","l","l"),pch=21,main="Fit for true y")
abline(h=mean(tau.fun(w,0)),lty=2)
###
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
