## Install latest version of XBART package if available
library(devtools)
install_github("JingyuHe/XBART@XBCF-RDD")
## Install latest version of HighDimRD package (KR)
#library(remotes)
#remotes::install_github("kolesarm/RDHonest")
#library(RDHonest)
install_github("Rafael-C-Alcantara/HighDimRD")
library(HighDimRD)

## Helper functions
### Extract 95% CI and posterior mean from pred.XBCFrd
tau.pred <- function(x,burnin,num_sweeps)
{
    ate.posterior <- colMeans(x$tau.adj[,(burnin+1):num_sweeps])
    return(c(quantile(ate.posterior,c(.025,.975)),mean(ate.posterior)))
}
### Minimize MSE for XBCF-RDD for a given grid of Owidth
min.mse <- function(h)
{
    pred.list <- vector("list",length(h))
    for (i in 1:length(h))
    {
        ## XBCF
        fit.XBCFrd <- XBCF.rd(y, w, x, c, Owidth = h[i], Omin = Omin, Opct = Opct,
                              pcat_con = 0, pcat_mod = 0,
                              num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                              num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                              burnin = burnin, Nmin = Nmin)

### Given the original (w,x) used for fitting the model, we take (w,x=0) as the test set
        pred.list[[i]] <- predict.XBCFrd(fit.XBCFrd, W = w, X = rep(0,n))
    }

    ## ATE summary for pred.list
    return(do.call("rbind",lapply(pred.list,tau.pred,burnin,num_sweeps)))
}
## Setup
set.seed(000)
library(XBART)

## Example data
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1)
n <- 2000
p <- 2
c <- 0 # Cutoff
w <- matrix(rnorm(n*p), n, p)
x <- rnorm(n,sd=.5)
z <- x >= c
y <- mu(w, x) + tau(w, x)*z + rnorm(n, 0, 0.1)
true.ate <- mean(tau(w,0))
## Estimation settings
### Owidth: overlap bandiwdth
### Omin: minimum number of observations inside overlap region for each leaf node
### Opct: If nb is the number of obs in node b, Opct*nb of them have to be inside overlap region
Owidth        <- seq(0.01,0.1,by=0.005)
Omin          <- 5
Opct          <- 0.9
num_trees_mod <- 10
num_trees_con <- 10
num_cutpoints <- n
Nmin          <- 20
num_sweeps    <- 20
burnin        <- 10

### Test run to check reasonable range for Owidth
ate <- min.mse(Owidth)
ylim <- c(min(range(true.ate,ate)), max(range(true.ate,ate)))
##
matplot(Owidth,ate,ylim=ylim,type=c("l","l","b"),lty=c(2,2,1),pch=16,cex=0.75,col="blue")
abline(h=mean(tau(w,0)),lty=2)
### Get 3 smallest MSE for that run and set Owidth to that
ind <- (ate[,3]-true.ate)^2
ind <- which(ind %in% head(sort(ind),n=3))
Owidth <- Owidth[ind]
Owidth
## Simulations
s <- 10 ## samples
###
results <- list(mse=NA,cont.tau=NA,cont.zero=NA,int.length=NA,pe=NA)
dims <- list(NULL,c("XBCF-RDD (1)","XBCF-RDD (2)","XBCF-RDD (3)",
                    "CGS","KR","FH"))
results$mse <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$cont.tau <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$cont.zero <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$int.length <- matrix(0,nrow=s,ncol=6,dimnames=dims)
results$pe <- matrix(0,nrow=s,ncol=6,dimnames=dims)
###
for (i in 1:s)
{
    w <- matrix(rnorm(n*p), n, p)
    x <- rnorm(n,sd=.5)
    z <- x >= c
    y <- mu(w, x) + tau(w, x)*z + rnorm(n, 0, 0.1)
    true.ate <- mean(tau(w,0))
    ## Get interactions and expansions for KR
    w1 <- fourier_basis(w,4)
    w_HighDim <- cbind(w,interaction_terms(w),w1,interaction_terms(w1))
    ## Estimation
    ate.xbcf <- min.mse(Owidth)
    ate.kr  <- HighDim_rd(y,x,w_HighDim,tpc="CV" ,rd="robust")
    ate.kr <- c(ate.kr$rd$Estimate[,"tau.bc"],
                ate.kr$rd$Estimate[,"tau.bc"]-1.96*ate.kr$rd$Estimate[,"se.rb"],
                ate.kr$rd$Estimate[,"tau.bc"]+1.96*ate.kr$rd$Estimate[,"se.rb"])
    ## Store results
### mse
    results$mse[i,1] <- (ate.xbcf[1,3]-true.ate)^2
    results$mse[i,2] <- (ate.xbcf[2,3]-true.ate)^2
    results$mse[i,3] <- (ate.xbcf[3,3]-true.ate)^2
    results$mse[i,"KR"] <- (ate.kr[1]-true.ate)^2
### cont.tau
    results$cont.tau[i,1] <- true.ate >= ate.xbcf[1,1] & true.ate <= ate.xbcf[1,2]
    results$cont.tau[i,2] <- true.ate >= ate.xbcf[2,1] & true.ate <= ate.xbcf[2,2]
    results$cont.tau[i,3] <- true.ate >= ate.xbcf[3,1] & true.ate <= ate.xbcf[3,2]
    results$cont.tau[i,"KR"] <- true.ate >= ate.kr[2] & true.ate <= ate.kr[3]
### cont.zero
    results$cont.zero[i,1] <- 0 >= ate.xbcf[1,1] & 0 <= ate.xbcf[1,2]
    results$cont.zero[i,2] <- 0 >= ate.xbcf[2,1] & 0 <= ate.xbcf[2,2]
    results$cont.zero[i,3] <- 0 >= ate.xbcf[3,1] & 0 <= ate.xbcf[3,2]
    results$cont.zero[i,"KR"] <- 0 >= ate.kr[2] & 0 <= ate.kr[3]
### int.length
    results$int.length[i,1] <- ate.xbcf[1,2] - ate.xbcf[1,1]
    results$int.length[i,2] <- ate.xbcf[2,2] - ate.xbcf[2,1]
    results$int.length[i,3] <- ate.xbcf[3,2] - ate.xbcf[3,1]
    results$int.length[i,"KR"] <- ate.kr[3] - ate.kr[2]
### pe (point.estimate)
    results$pe[i,1] <- ate.xbcf[1,3]
    results$pe[i,2] <- ate.xbcf[2,3]
    results$pe[i,3] <- ate.xbcf[3,3]
    results$pe[i,"KR"] <- ate.kr[1]
}
## Load CGS results and merge
cgs <- readRDS("results_cgs.rds")
## Plot results
par(mfrow=c(1,2))
boxplot(sqrt(results$mse+cgs$mse),cex.axis=.75)
boxplot(results$pe+cgs$pe,cex.axis=.75)
par(mfrow=c(1,1))
