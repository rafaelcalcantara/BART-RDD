library(devtools)
## Install CGS package if on Windows
if (Sys.info()["sysname"] == "Windows")
{
    ## Install and load CGS package
    install.packages("http://apps.olin.wustl.edu/faculty/chib/rpackages/rdd/bayesrdd_1.0.zip",
                     repo=NULL,source=T)
    library(bayesrdd)
} else
{
  ## Install latest version of XBART package if available
  install_github("rafaelcalcantara/XBART@XBCF-RDD")
  library(XBART)
  ## Install latest version of HighDimRD package (KR)
  install_github("kolesarm/RDHonest")
  library(RDHonest)
  install_github("akreiss/HighDimRD")
  library(HighDimRD)
}
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
                          tau_con = var(y)/num_trees_con,
                          tau_mod = var(y)/num_trees_mod,
                          random_seed=0)

### Given the original (w,x) used for fitting the model, we take (w,x=0) as the test set
    pred <- predict.XBCFrd(fit.XBCFrd, W = w, X = rep(0,n))

    ## ATE summary for pred.list
    return(tau.pred(pred,burnin,num_sweeps))
}

#### For CGS with no covariates
pred.ate.xbcf0 <- function()
{
    ## XBCF
    fit.XBCFrd <- XBCF.rd(y, NULL, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                          num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                          num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                          burnin = burnin, Nmin = Nmin,
                          p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                          tau_con = var(y)/num_trees_con,
                          tau_mod = var(y)/num_trees_mod,
                          random_seed=0)

### Given the original (w,x) used for fitting the model, we take (w,x=0) as the test set
    pred <- predict.XBCFrd(fit.XBCFrd, W = NULL, X = rep(0,n))

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
                   tau_con = var(y)/num_trees_con,
                   tau_mod = var(y)/num_trees_mod,
                   burnin = 0, Nmin = Nmin)
    trees_json <- jsonlite::fromJSON(fit$tree_json_mod,simplifyVector=F)
    return(trees_json$trees[["0"]][["0"]][["avail.cutpoints"]])
}
### Find minimum Owidth that leads to more than 1 available cutppoint at root node
findOwidth <- function(seq)
{
    i <- 0.01
    ## Do this to avoid i s.t. there are no obs in the bandwidth
    while(sum(x>=-i & x<0)<0.04*n & sum(x>=0 & x<=i)<0.04*n) i <- i+0.01
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
    return(do.call("cbind",
                   lapply(list(xbcf,cgs,kr,fh),
                          function(x) do.call("rbind",
                                              lapply(x, function(y) mean((y[,1]-y[,4])^2))
                                              )
                          )
                   )
           )
}
## DGP functions
### 1) CGS
dgp.cgs <- function(n,p)
{
    ## Additional covariates are unnecessary
    w  <- matrix(rnorm(n*p),n,p)
    x  <- 2*rbeta(n,2,4)-1
    z  <- x>=0
    mu <- function(x)
    {
        0.48 + 1.27*x + 7.18*x^2 + 20.21*x^3 + 21.54*x^4 + 7.33*x^5
    }
    tau <- function(x)
    {
        0.04 - 0.43*x - 10.18*x^2 - 12.22*x^3 - 30.55*x^4 - 3.77*x^5
    }
    y  <- mu(x) + tau(x)*z + 0.1295*rt(n,3)
    ate <- 0.04
    tau.x <- tau(x)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
### 2) FH
dgp.fh <- function(n,p=2)
{
    ## We consider only the case where the covariates affect Y but don't change across treatment states
    x <- 2*rbeta(n,2,4)-1
    z <- x >= 0
    w <- matrix(rnorm(n*p,sd=0.5),n,p)
    a <- 2/(1:p)^2
    mu <- function(x,w)
    {
        0.5*x + 0.25*x^2 + 0.4*w%*%a + 0.2*w^2%*%a
    }
    tau <- function(x)
    {
        1 - 0.25*x
    }
    y <- mu(x,w) + tau(x)*z + rnorm(n,sd=0.25)
    ate <- 1
    tau.x <- tau(x)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
### 3) KR
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
    mu <- function(x,w,a)
    {
        0.36 + 0.96*x + 5.47*x^2 + 15.28*x^3 + 5.14*x^5 + as.vector(w%*%a)
    }
    tau <- function(x,w,a)
    {
        0.02 - 0.34*x - 8.31*x^2 - 6.86*x^3 - 0.83*x^5 + as.vector(0.06*w%*%a)
    }
    y  <- mu(x,w,a) + tau(x,w,a)*z + e
    ate <- mean(tau(0,w,a))
    tau.x <- tau(x,w,a)
    return(list(y=y,x=x,z=z,w=w,ate=ate,tau.x=tau.x))
}
## Generate 1 dataset
dgp <- function(n,p,type)
{
    if (type=="cgs") data <- dgp.cgs(n,p)
    if (type=="kr") data <- dgp.kr(n,p)
    if (type=="fh") data <- dgp.fh(n,p)
    list2env(data,globalenv())
}
