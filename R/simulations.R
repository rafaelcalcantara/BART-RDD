## Install latest version of XBART package if available
library(devtools)
install_github("JingyuHe/XBART@XBCF-RDD")
## Install latest version of HighDimRD package (KR)
#library(remotes)
#remotes::install_github("kolesarm/RDHonest")
#library(RDHonest)
install_github("Rafael-C-Alcantara/HighDimRD")
library(HighDimRD)
## Source code from FH
## install.packages(c("np","rdd","matrixStats","xtable","boot"))
## source("R_code_for_Including_covariates_in_the_regression_.R")
## Helper functions
### Extract 95% CI and posterior mean from pred.XBCFrd
tau.pred <- function(x,burnin,num_sweeps)
{
    ate.posterior <- colMeans(x$tau.adj[,(burnin+1):num_sweeps])
    return(c(mean(ate.posterior),quantile(ate.posterior,c(.025,.975))))
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
### Table results
resTab <- function(x)
{
    
    return(do.call("rbind",lapply(x,colMeans)))
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
matplot(Owidth,ate,ylim=ylim,type=c("l","b","l"),lty=c(2,1,2),pch=16,cex=0.75,col="blue")
abline(h=mean(tau(w,0)),lty=2)
### Get 3 smallest MSE for that run and set Owidth to that
ind <- (ate[,3]-true.ate)^2
ind <- which(ind %in% head(sort(ind),n=3))
Owidth <- Owidth[ind]
Owidth
## Simulations
s <- 100 ## samples
###
results <- list("XBCF-RDD (1)"=matrix(0,s,3),"XBCF-RDD (2)"=matrix(0,s,3),
                "XBCF-RDD (3)"=matrix(0,s,3),"CGS"=NA,"KR"=matrix(0,s,3),"FH"=NA)
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
    ### ate.fh <- rdd.x(y,x,w)
    ## Store results
    results[["XBCF-RDD (1)"]][i,] <- ate.xbcf[1,]
    results[["XBCF-RDD (2)"]][i,] <- ate.xbcf[2,]
    results[["XBCF-RDD (3)"]] <- ate.xbcf[3,]
    results[["KR"]][i,] <- ate.kr
    ## results[["FH"]][i,] <- ate.fh
}
## Load CGS results and merge
cgs <- readRDS("results_cgs.rds")
results[["CGS"]] <- cgs$CGS
results[["FH"]] <- cgs$FH
saveRDS(results,"results.rds")
## Plot MSE
# boxplot(sqrt(results$mse+cgs$mse),cex.axis=.75)
## Table results
# round(resTab(results)+resTab(cgs),digits=4)
