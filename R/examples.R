## Install latest version of XBART package if available
library(devtools)
install_github("JingyuHe/XBART@XBCF-RDD")

## Helper functions
tau.pred <- function(x,burnin,num_sweeps)
{
    ## Extract 95% CI and posterior mean from pred.XBCFrd
    ate.posterior <- colMeans(x$tau.adj[,(burnin+1):num_sweeps])
    return(c(quantile(ate.posterior,c(.025,.975)),mean(ate.posterior)))
}

## Setup
set.seed(000)
library(XBART)
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1) # make sure the treatment effect is non-zero

n <- 2000
p <- 2
c <- 0 # Cutoff

## Data
w <- matrix(rnorm(n*p), n, p)
x <- rnorm(n,sd=.5)
z <- x >= c
y <- mu(w, x) + tau(w, x)*z + rnorm(n, 0, 0.1)

## Estimation settings
### Owidth: overlap bandiwdth
### Omin: minimum number of observations inside overlap region for each leaf node
### Opct: If nb is the number of obs in node b, Opct*nb of them have to be inside overlap region
Owidth        <- 0.1
Omin          <- 5
Opct          <- 0.7
num_trees_mod <- 10
num_trees_con <- 10
num_cutpoints <- n
Nmin          <- 20
num_sweeps    <- 20
burnin        <- 10

## XBCF
fit.XBCFrd <- XBCF.rd(y, w, x, c, Owidth = Owidth, Omin = Omin, Opct = Opct,
                      pcat_con = 0, pcat_mod = 0,
                      num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                      num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin)

## Predictions
### Predict tau for a grid of X values with W fixed
grid <- seq(-1,1,by=.1)
pred.grid <- matrix(0,nrow=length(grid),ncol=3)
true.ate <- rep(0,length(grid))
for (i in 1:length(grid))
{
    true.ate[i] <- mean(tau(w,grid[i]))
    pred.grid[i,] <- tau.pred(predict.XBCFrd(fit.XBCFrd, W = w, X = rep(grid[i],n)),
                              burnin,num_sweeps)
}

plot(grid,true.ate,type="l")
matlines(grid,pred.grid,lty=rep(1,3),lwd=c(1,1,1.5),col="cadetblue3")
abline(v=0,lty=2)
### Given the original (w,x) used for fitting the model, we take (w,x=0) as the test set
pred.XBCFrd <- predict.XBCFrd(fit.XBCFrd, W = w, X = rep(0,n))
