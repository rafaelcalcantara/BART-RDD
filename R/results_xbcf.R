library(XBART)
setwd("~/Documents/Git/XBCF-RDD/")
## Weights function
sample.tau.ind <- function(tau,wts)
{
    ## Function to calculate tau for a given sweep
    p <- apply(wts$yhats,2,sample,1)
    gamma <- rbinom(length(p),1,p)
    return(gamma*tau/sum(gamma))
}
## CHECK IF SAMPLE.TAU.SAMPLE IS RIGHT
sample.tau.sample <- function(n,pred,weights)
{
    ## Function to calculate tau for every sweep of a given sample
    ate <- sapply(n,function(x) sample.tau.ind(pred$tau.adj[,x],weights))
    return(colSums(ate))
}
sample.tau <- function(n,pred,weights)
{
    return(mapply(function(x,y) sample.tau.sample(n,x,y),pred,weights))
}
## DGP2
dgp <- readRDS("Data/DGP2.rds")
dgp <- dgp[1:50]
s <- 50
results <- list()
for (i in 1:s)
{
    results[[i]] <- readRDS(paste0("Results/xbcf_dgp2_",i,".rds"))
}
weights <- list()
for (i in 1:s)
{
    weights[[i]] <- readRDS(paste0("Results/weights_dgp2_",i,".rds"))
}
### Obtain ATE posterior
pred <- mapply(function(x,y) predict.XBCFrd(x$fit,y$w,rep(0,500)),results,dgp,SIMPLIFY=F)
post1 <- sapply(pred["tau.adj",],colMeans)
post2 <- mapply(function(x,y) colSums(y$wts*x),pred["tau.adj",],results)
post3 <- mapply(function(x,y) colMeans(x[y==1,]),pred["tau.adj",],tau.sample)
post1.sum <- t(apply(post1,2,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
post2.sum <- t(apply(post2,2,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
post3.sum <- t(apply(post3,2,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
###
mean((post3.sum[,1]-1)^2)
mean(post3.sum[,2]<=1 & 1<=post3.sum[,3])
mean(post3.sum[,3]-post3.sum[,2])
