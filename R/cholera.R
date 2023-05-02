devtools::install_github("Rafael-C-Alcantara/XBART@XBCF-RDD")
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
                          burnin = burnin, Nmin = Nmin)

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
                   burnin = 0, Nmin = Nmin,
                   p_categorical_con = p_categorical,
                   p_categorical_mod = p_categorical)
    trees_json <- jsonlite::fromJSON(fit$tree_json_mod,simplifyVector=F)
    return(trees_json$trees[["0"]][["0"]][["avail.cutpoints"]])
}
### Find minimum Owidth that leads to more than 1 available cutppoint at root noden
findOwidth <- function(seq)
{
    i <- 0.01
    ## Do this to avoid i s.t. there are no obs in the bandwidth
    while(sum(x>=-i & x<=i)<30) i <- i+0.01
    while(cutpoints(i)==1) i <- i+seq
    return(i)
}
##
data <- readstata13::read.dta13("R/Merged_1853_1864_data.dta")
data$dist_netw <- data$dist_netw/100
data$dist_2 <- ifelse(data$broad==0,-data$dist_netw,data$dist_netw)
vars <- c("log_rentals_1864","dist_netw","dist_2",
          "dist_cent","dist_urinal",
          "dist_square","dist_fire",
          "dist_thea","dist_police",
          "dist_pub","dist_church",
          "dist_bank","dist_vent",
          "dist_pump","dist_pit_fake",
          "no_sewer","old_sewer")
covs <- c("dist_cent","dist_urinal",
          "dist_square","dist_fire",
          "dist_thea","dist_police",
          "dist_pub","dist_church",
          "dist_bank","dist_vent",
          "dist_pump","dist_pit_fake",
          "no_sewer","old_sewer")
data <- data[,vars]
data <- data[complete.cases(data),]
y <- data$log_rentals_1864
x <- data$dist_2
w <- data[,covs]
n <- nrow(data)
c <- 0
##
bins <- aggregate(y,list(bins=cut(x,seq(-1,1,0.05),right=F)),mean)
par(mfrow=c(1,2))
plot(seq(-100,95,5),bins$x,type="n",
     xlab="Distance to Broad St. Pump",ylab="Log rental prices in 1864")
points(seq(-100,-5,5),bins[1:20,2],pch=19)
lines(seq(-100,-5,5),bins[1:20,2])
points(seq(0,95,5),bins[21:40,2],pch=19)
lines(seq(0,95,5),bins[21:40,2])
abline(v=0,lty=2)
hist(100*x,main="",xlab="Distance to Broad St. Pump",ylab="")
par(mfrow=c(1,1))
##
cct1 <- rdrobust::rdrobust(y,x)
cct2 <- rdrobust::rdrobust(y,x,covs=w)
##
Omin          <- 2
Opct          <- 0.7
num_trees_mod <- 20
num_trees_con <- 10
num_cutpoints <- n
Nmin          <- 10
num_sweeps    <- 50
burnin        <- 20
p_categorical <- 2
## Owidth        <- findOwidth(0.01)
Owidth <- 1
##
fit.XBCFrd <- XBCF.rd(y, w, x, c,
                      Owidth = Owidth, Omin = Omin, Opct = Opct,
                      num_trees_mod = num_trees_mod, num_trees_con = num_trees_con,
                      num_cutpoints = num_cutpoints, num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                      random_seed=0)

pred <- predict.XBCFrd(fit.XBCFrd, W = w, X = rep(0,n))
pred.tot <- predict.XBCFrd(fit.XBCFrd, W = w, X = x)
tau <- pred$tau.adj[,(burnin+1):num_sweeps]
mu  <- pred.tot$mu.adj[,(burnin+1):num_sweeps]
ate.xbcf <- colMeans(tau)
##
tau.pred(pred,burnin,num_sweeps)
## Plot individual treatment effects and mu
ind.mu <- t(apply(mu,1,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
ind.mu <- ind.mu[order(ind.mu[,1]),]
ite <- t(apply(tau,1,function(x) c(mean(x),quantile(x,c(0.025,0.975)))))
ite <- ite[order(ite[,1]),]
cart.mu  <- rpart::rpart(mu~.,data=cbind(mu=rowMeans(mu),w))
cart.tau <- rpart::rpart(tau~.,data=cbind(tau=rowMeans(tau),w))
table(ite[,3]<0)
###
matplot(ind.mu,type="l",lty=c(1,2,2),lwd=c(2,1,1),col="cadetblue",
        xlab="House",ylab=expression(mu))
abline(h=0,lty=2)
matplot(ite,type="l",lty=c(1,2,2),lwd=c(2,1,1),col="cadetblue",
        xlab="House",ylab=expression(tau))
abline(h=0,lty=2)
rpart.plot::rpart.plot(cart.mu)
rpart.plot::rpart.plot(cart.tau)
