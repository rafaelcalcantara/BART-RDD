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
p_categorical <- 1
##
### DGP
sig <- 0.15
## x <- matrix(0.6*rnorm(n*s2,-1,0.5) + 0.4*rnorm(n*s2,1,0.5),n,s2)
## x <- matrix(rnorm(n*s2),n,s2)
x <- matrix(0.4*rnorm(n*s2,-1,0.5) + 0.6*rnorm(n*s2,1,0.5),n,s2)
w <- matrix(rbinom(n*s2,1,0.7)+0.5,n,s2)
z <- apply(x,2,function(i) as.numeric(i>=c))
e <- matrix(rnorm(n*s2,0,sig),n,s2)
mu <- function(x,w) 0.2*w + 0.1*w*x - 0.2*w*x^2 + 0.1*w*x^3
tau <- function(x,w) 0.06*w - 0.2*w*x - 0.5*w*x^2 + 0.5*w*x^3
y <- mu(x,w) + tau(x,w)*z + e
plot(x[,1],y[,1],bty="n",pch=20,col=z[,1]+1,xlab="X",ylab="Y")
## Prior predictive
pars <- rnorm(8,0,0.1)
mu.prior <- function(x,w) as.vector(cbind(1,x,x^2,x^3)%*%pars[1:4])
tau.prior <- function(x,w) as.vector(cbind(1,x,x^2,x^3)%*%pars[5:8])
ys <- vector("list",s2)
for (i in 1:s2) ys[[i]] <- mu.prior(x[,i],w[,i]) + tau.prior(x[,i],w[,i])*z[,i] + matrix(rnorm(n*s1,0,0.1),n,s1)
plot(x[,1],ys[[1]][,1],bty="n",pch=20,col=z[,1]+1,xlab="X",ylab=expression(Y[prior]))
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
h <- round(rbind(0.5*apply(x,2,sd),apply(x,2,sd),2*apply(x,2,sd)),2)
h1 <- vector("list",s2)
bart.prior.time <- system.time({
    for (j in 1:s2) h1[[j]] <- opt.h(s1,ys[[j]],x[,j],w[,j],z[,j],h[,j])
})
## saveRDS(h1,"Results/prior_h.rds")
## h1 <- readRDS("Results/prior_h.rds")
out <- matrix(0,s2,nrow(h))
for (i in 1:s2)
{
    for (j in 1:nrow(h))
    {
        out[i,j] <- sqrt(mean(sapply(h1[[i]][[j]],function(x) (mean(colMeans(x))-mean(tau.prior(0,w[,i])))^2)))
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
ate <- sapply(1:ncol(w),function(i) mean(tau(0,w[,i])))
bart.time <- system.time({
    bart.fit <- fit(s2,h.opt,y,x,w,z)
})
## saveRDS(bart.fit,"Results/prior_test_bart_rdd.rds")
## bart.fit <- readRDS("Results/prior_test_bart_rdd.rds")
ate.bart <- colMeans(sapply(bart.fit,colMeans))
rmse.bart <- sqrt(mean((ate.bart-ate)^2))
bias.bart <- mean(ate.bart-ate)
var.bart <- var(ate.bart)
ci.bart <- apply(sapply(bart.fit,colMeans),2,quantile,c(0.025,0.975))
ci.bart <- t(ci.bart)
###
ckt.fit <- function(s,y,x,w)
{
    foreach(i=1:s,.multicombine=T,.export="c") %dopar%
        {
            rdrobust(y[,i],x[,i],c=c,covs=cbind(w[,i],x[,i]*w[,i]))
        }
}
ckt.time <- system.time({
    ckt <- ckt.fit(s2,y,x,w)
})
## saveRDS(ckt,"Results/prior_test_polynomial.rds")
## ckt <- readRDS("Results/prior_test_polynomial.rds")
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
print(paste0("BART prior tuning time: ",bart.prior.time))
print(paste0("BART time: ",bart.time))
print(paste0("Polynomial time: ",ckt.time))
print("Done!")
