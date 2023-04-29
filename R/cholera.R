data <- readstata13::read.dta13("R/Merged_1853_1864_data.dta")
data$dist_netw <- data$dist_netw/100
data$dist_2 <- ifelse(data$broad==0,-data$dist_netw,data$dist_netw)
data <- data[,c("deaths","death_ind","dist_netw","dist_2",
                "dist_cent","dist_urinal","no_sewer","old_sewer")]
data <- data[complete.cases(data),]
y <- data$deaths
x <- data$dist_2
w <- data[,c("dist_cent","dist_urinal","no_sewer","old_sewer")]
n <- nrow(data)
c <- 0
##
cct1 <- rdrobust::rdrobust(y,x)
cct2 <- rdrobust::rdrobust(y,x,covs=w)
##
poly <- lm(deaths~broad+dist_netw+dist_netw2+dist_netw3+dist_cent+dist_urinal+no_sewer+old_sewer,data,subset=data$dist_netw<cct1$bws[1,1])
##
library(XBART)
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
                          pcat_con = 0, pcat_mod = 0,
                          num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                          num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                          burnin = burnin, Nmin = Nmin,
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
                   pcat_con = 0, pcat_mod = 0,
                   num_trees_mod = 1, num_trees_con = 1,
                   num_cutpoints = num_cutpoints, num_sweeps = 1,
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
##
Omin          <- 2
Opct          <- 0.7
num_trees_mod <- 10
num_trees_con <- 10
num_cutpoints <- n
Nmin          <- 10
num_sweeps    <- 30
burnin        <- 20
Owidth        <- findOwidth(0.01)
##
res <- matrix(0,50,3)
for (i in 1:50)
{
    Owidth <- i/100
    res[i,] <- pred.ate.xbcf()
}
res
matplot(res,type=c("b","l","l"),lty=2,pch=19,col="blue")
abline(h=0,lty=2)
