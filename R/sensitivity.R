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
### mu and tau functions
mu <- function(W, X) return(0.1 * rowSums(W) + 1/(1+exp(-5*X)))
tau <- function(W, X) return( sin(mu(W, X)) + 1)
### Plotting function
plot.res <- function(Owidth,y,dist)
{
    matplot(x=Owidth,y=y,type="b",lty=1,col=c("brown2","darkgoldenrod2","dodgerblue2"),main=dist,pch=16)
    legend("topright",col=c("brown2","darkgoldenrod2","dodgerblue2"),legend=c(100,250,500),lty=1)
}
## Estimation settings
### Owidth: overlap bandiwdth
### Omin: minimum number of observations inside overlap region for each leaf node
### Opct: If nb is the number of obs in node b, Opct*nb of them have to be inside overlap region
Owidth        <- c(0.05,0.1,0.15,0.2)
Omin          <- 10
Opct          <- 0.7
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
obs <- list("Gaussian"=vector("list",length(Owidth)),
            "Uniform"=vector("list",length(Owidth)))
###
for(i in 1:2)
{
    for (j in 1:length(Owidth))
    {
        obs[[i]][[j]] <- matrix(0,s,length(n))
    }
}
###
results <- list("Gaussian"=matrix(0,length(Owidth),length(n)),
                "Uniform"=matrix(0,length(Owidth),length(n)))
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
        temp.results <- matrix(0,length(Owidth),s)
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
            ate <- ate.xbcf(Owidth,n[i])
            temp.results[,j] <- (ate-true.ate)^2
            ## Count observations
            for (k in 1:length(Owidth))
            {
                obs[[h]][[k]][j,i] <- sum((x>=-Owidth[k] & x<=Owidth[k]))/n[i]
            }
            ###
        }
        results[[h]][,i] <- rowMeans(temp.results)
    }
}
###
png("MSE_width.png")
par(mfrow=c(1,2))
plot.res(Owidth,results$Gaussian,"Gaussian")
plot.res(Owidth,results$Uniform,"Uniform")
par(mfrow=c(1,1))
dev.off()
###
minMSE <- do.call("rbind",lapply(results,function(x) apply(x,2,min)))
colnames(minMSE) <- n
minWidth <- do.call("rbind",lapply(results,function(x) apply(x,2,function(x) Owidth[which(x==min(x))])))
colnames(minWidth) <- n
###
obs <- lapply(obs,function(x) lapply(x,colMeans))
obs <- 
