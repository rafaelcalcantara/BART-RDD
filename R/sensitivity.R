set.seed(000)
###
library(devtools)
## Install latest version of XBART package if available
install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
## Helper functions
### Extract 95% CI and posterior mean from pred.XBCFrd
tau.pred <- function(x,burnin,num_sweeps)
{
    ate.posterior <- colMeans(x$tau.adj[,(burnin+1):num_sweeps])
    return(mean(ate.posterior))
}
### Minimize MSE for XBCF-RDD for a given grid of Owidth
ate.xbcf <- function(h,n)
{
    pred.list <- vector("list",length(h))
    for (i in 1:length(h))
    {
        ## XBCF
        fit.XBCFrd <- XBCF.rd(y, w, x, c, Owidth = 0.1, Omin = Omin, Opct = Opct[i],
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
### mu and tau functions
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1)
### Plotting function
plot.res <- function(Opct,y,dist)
{
    matplot(x=Opct,y=y,type="b",lty=1,col=c("brown2","darkgoldenrod2","dodgerblue2"),main=dist,pch=16)
    legend("topright",col=c("brown2","darkgoldenrod2","dodgerblue2"),legend=c(100,250,500),lty=1)
}
## Estimation settings
### Owidth: overlap bandiwdth
### Omin: minimum number of observations inside overlap region for each leaf node
### Opct: If nb is the number of obs in node b, Opct*nb of them have to be inside overlap region
Owidth        <- 0.1
Omin          <- 2
Opct          <- c(0.3,0.5,0.7,0.9)
num_trees_mod <- 10
num_trees_con <- 10
Nmin          <- 10
num_sweeps    <- 20
burnin        <- 10
## Example data
n <- c(100,250,500)
p <- 2
c <- 0 # Cutoff
s <- 50
obs <- list("Gaussian"=vector("list",length(Opct)),
            "Uniform"=vector("list",length(Opct)))

results <- list("Gaussian"=matrix(0,length(Opct),length(n)),
                "Uniform"=matrix(0,length(Opct),length(n)))
for (h in 1:2)
{
    if (h==1)
    {
        dist.x <- function(x) rnorm(x,sd=.5)
        d <- "Gaussian"
    } else if (h==2)
    {
        dist.x <- function(x) runif(x,-2,2)
        d <- "Uniform"
    }
    for (i in 1:length(n))
    {
        num_cutpoints <- n[i]
        temp.results <- matrix(0,length(Opct),s)
        for (j in 1:s)
        {
            print(paste0("Distibution of x: ",d))
            print(paste0("Sample size: ",n[i]))
            print(paste0("Simulation ",j))
            w <- matrix(rnorm(n[i]*p), n[i], p)
            x <- dist.x(n[i])
            z <- x >= c
            y <- mu(w, x) + tau(w, x)*z + rnorm(n[i], 0, 0.1)
            true.ate <- mean(tau(w,0))
            ate <- ate.xbcf(Opct,n[i])
            temp.results[,j] <- (ate-true.ate)^2
        }
        results[[h]][,i] <- rowMeans(temp.results)
    }
}
###
png("MSE_opct.png")
par(mfrow=c(1,2))
plot.res(Opct,results$Gaussian,"Gaussian")
plot.res(Opct,results$Uniform,"Uniform")
par(mfrow=c(1,1))
dev.off()
###
minMSE <- do.call("rbind",lapply(results,function(x) apply(x,2,min)))
colnames(minMSE) <- n
minOpct <- do.call("rbind",lapply(results,function(x) apply(x,2,function(x) Opct[which(x==min(x))])))
colnames(minOpct) <- n
