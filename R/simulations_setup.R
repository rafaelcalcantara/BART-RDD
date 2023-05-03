library(devtools)
## Install latest version of XBART package if available
install_github("Rafael-C-Alcantara/XBART@XBCF-RDD")
library(XBART)
## Install CGS package if on Windows
if (Sys.info()["sysname"] == "Windows")
{
    ## Install and load CGS package
    install.packages("http://apps.olin.wustl.edu/faculty/chib/rpackages/rdd/bayesrdd_1.0.zip",
                     repo=NULL,source=T)
    library(bayesrdd)
}
## Install latest version of HighDimRD package (KR)
install_github("kolesarm/RDHonest")
library(RDHonest)
install_github("akreiss/HighDimRD")
library(HighDimRD)
## Source code from FH
### install.packages(c("np","rdd","matrixStats","xtable","boot"))
source("R/FH.R")
## Helper functions
### Extract 95% CI and posterior mean from pred.XBCFrd
tau.pred <- function(x,burnin,num_sweeps)
{
    ate.posterior <- colMeans(x$tau.adj[,(burnin+1):num_sweeps])
    return(c(mean(ate.posterior),quantile(ate.posterior,c(.025,.975))))
}
### Return ATE prediction for XBCF-RDD
pred.ate.xbcf <- function()
{
    ## XBCF
    fit.XBCFrd <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                          num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                          num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                          burnin = burnin, Nmin = Nmin,
                          p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                          random_seed=0)

### Given the original (w,x) used for fitting the model, we take (w,x=0) as the test set
    pred <- predict.XBCFrd(fit.XBCFrd, W = w, X = rep(0,n))

    ## ATE summary for pred.list
    return(tau.pred(pred,burnin,num_sweeps))
}
### Return available cutpoints at root node for a given Owidth
cutpoints <- function(Owidth)
{
    fit <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                   num_trees_mod = 1, num_trees_con = 1,
                   num_cutpoints = num_cutpoints, num_sweeps = 1,
                                             p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                   burnin = 0, Nmin = Nmin)
    trees_json <- jsonlite::fromJSON(fit$tree_json_mod,simplifyVector=F)
    return(trees_json$trees[["0"]][["0"]][["avail.cutpoints"]])
}
### Find minimum Owidth that leads to more than 1 available cutppoint at root noden
findOwidth <- function(seq)
{
    i <- 0.01
    ## Do this to avoid i s.t. there are no obs in the bandwidth
    while(sum(x>=-i & x<=i)<15) i <- i+0.01
    while(cutpoints(i)==1) i <- i+seq
    return(i)
}
### Plot results
resPlot <- function(x)
{
    plotMat <- do.call("cbind",lapply(x[-1],function(y) (y[,1]-x$ATE)^2))
    boxplot(plotMat)
}
### Table results
resTab <- function(x)
{
    return(colMeans(do.call("cbind",lapply(x[-1],function(y) (y[,1]-x$ATE)^2))))
}
## DGP functions
### 1) CGS
dgp.cgs <- function(n,p)
{
    ## Additional covariates are unnecessary
    w  <- matrix(rnorm(n*p),n,p)
    x  <- 2*rbeta(n,2,4)-1
    z  <- x>=0
    g0 <- 0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5
    g1 <- 0.52 + 0.84*x - 3*x^2 + 7.99*x^3 - 9.01*x^4 + 3.56*x^5
    y  <- z*g1 + (1-z)*g0 + 0.1295*rt(n,3)
    tau <- 0.04
    return(list(y=y,x=x,z=z,w=w,tau=tau))
}
### 2) KR
dgp.kr <- function(n,p=200)
{
    sig.e <- 0.1295^2
    sig.z <- 0.1353^2
    v     <- (0.8*sqrt(6)*sig.e/pi)/1:p
    var   <- matrix(0, p+1, p+1)
    var[1,1] <- sig.e
    var[1,2:(p+1)] <- v
    var[2:(p+1),1] <- v
    var[2:(p+1),2:(p+1)] <- sig.z*diag(p)
    eZ <- mnormt::rmnorm(n,varcov=var)
    e <- eZ[,1]
    w <- eZ[,-1]
    a <- 2/(1:p)^2
    x  <- 2*rbeta(n,2,4)-1
    z  <- x>=0
    y0 <- 0.36 + 0.96*x + 5.47*x^2 + 15.28*x^3 + 15.87*x^4 + 5.14*x^5 + 0.22*w%*%a
    y1 <- 0.38 + 0.62*x - 2.84*x^2 + 8.42*x^3 - 10.24*x^4 + 4.31*x^5 + 0.28*w%*%a
    y  <- as.vector(z*y1 + (1-z)*y0 + e)
    tau <- 0.02
    return(list(y=y,x=x,z=z,w=w,tau=tau))
}
### 3) FH
dgp.fh <- function(n,p=2)
{
    ## We consider only the case where the distribution of W varies by treatment state and covariates also affect Y
    x <- rnorm(n)
    latent.vars <- matrix(rnorm(3*n),n,3)
    z <- x >= 0
    w1 <- 0.2*z + 0.5*latent.vars[,1]
    w2 <- 0.2*z + 0.5*latent.vars[,2]
    w <- cbind(w1,w2)
    y <- z + 0.5*x - 0.25*z*x + 0.25*x^2 + 0.4*(w1+w2) + 0.2*(w1^2+w2^2) + latent.vars[,3]
    tau <- 1
    return(list(y=y,x=x,z=z,w=w,tau=tau))
}
## Generate 1 dataset
dgp <- function(n,p,type)
{
    if (type=="cgs") data <- dgp.cgs(n,p)
    if (type=="kr") data <- dgp.kr(n,p)
    if (type=="fh") data <- dgp.fh(n,p)
    list2env(data,globalenv())
}
