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
c             <- 0
Omin          <- 1
h             <- 0.1
Opct          <- 0.9
ntrees        <- 10
Nmin          <- 10
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 0
## h
s <- 11
## Data 1
n <- 500
x <- 2*rbeta(n,2,4)-1
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*w + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) +1) # make sure the treatment effect is non-zero
y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n)
fit <- function(s)
{
    tau <- rep(0,s)
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            h <- 0.05 + (i-1)/100
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
                           parallel=T,num_cores = no_cores)
            test <- -h<=x & x<=h
            pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
            pred <- pred$tau.adj[,(burnin+1):num_sweeps]
            pred <- colMeans(pred)
            mean(pred)
        }
}
##
tau1 <- fit(s)
tau1 <- unlist(tau1)
## Data 2
n <- 500
x <- 2*rbeta(n,2,4)-0.5
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*w + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) +1) # make sure the treatment effect is non-zero
y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n)
fit <- function(s)
{
    tau <- rep(0,s)
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            h <- 0.05 + (i-1)/100
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
                           parallel=T,num_cores = no_cores)
            test <- -h<=x & x<=h
            pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
            pred <- pred$tau.adj[,(burnin+1):num_sweeps]
            pred <- colMeans(pred)
            mean(pred)
        }
}
##
tau2 <- fit(s)
tau2 <- unlist(tau2)
## Data 3
n <- 500
x <- runif(n,-1,1)
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*w + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) +1) # make sure the treatment effect is non-zero
y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n)
fit <- function(s)
{
    tau <- rep(0,s)
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            h <- 0.1 + (i-1)/100
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
                           parallel=T,num_cores = no_cores)
            test <- -h<=x & x<=h
            pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
            pred <- pred$tau.adj[,(burnin+1):num_sweeps]
            pred <- colMeans(pred)
            mean(pred)
        }
}
##
tau3 <- fit(s)
tau3 <- unlist(tau3)
## Data 4
n <- 500
x <- rnorm(n,0,0.25)
w <- rnorm(n,0,0.25)
z <- x>=0
mu.fun <- function(W, X){return(0.1*w + 1/(1+exp(-5*X)))} 
tau.fun <- function(W, X) return( sin(mu.fun(W, X)) +1) # make sure the treatment effect is non-zero
y <- mu.fun(w, x) + tau.fun(w, x)*z + rnorm(n)
fit <- function(s)
{
    tau <- rep(0,s)
    foreach(i=1:s,.multicombine=T,.export=c("n","c","Omin","h","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            h <- 0.01 + (i-1)/100
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
                           parallel=T,num_cores = no_cores)
            test <- -h<=x & x<=h
            pred <- predict.XBCFrd(fit,w[test],rep(0,sum(test)))
            pred <- pred$tau.adj[,(burnin+1):num_sweeps]
            pred <- colMeans(pred)
            mean(pred)
        }
}
##
tau4 <- fit(s)
tau4 <- unlist(tau4)
##
par(mfrow=c(2,2))
plot(x=seq(0.05,0.15,0.01),y=tau1,type="b",col="blue",bty="n")
abline(h=mean(tau.fun(w,0)),lty=2)
plot(x=seq(0.05,0.15,0.01),y=tau2,type="b",col="blue",bty="n")
abline(h=mean(tau.fun(w,0)),lty=2)
plot(x=seq(0.01,0.11,0.01),y=tau3,type="b",col="blue",bty="n")
abline(h=mean(tau.fun(w,0)),lty=2)
plot(x=seq(0.1,0.2,0.01),y=tau4,type="b",col="blue",bty="n")
abline(h=mean(tau.fun(w,0)),lty=2)
