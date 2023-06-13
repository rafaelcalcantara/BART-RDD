## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("Rafael-C-Alcantara/XBART@XBCF-RDD")
library(XBART)
## Helper functions
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
fit.general <- function(h,y,w,x)
{
    foreach(i=1:length(h),.multicombine=T,.export=c("p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
        {
            fit <- XBCF.rd(y, w, x, c, Owidth = h[i], Omin = Omin, Opct = Opct,
                           num_trees_mod = m, num_trees_con = m,
                           num_cutpoints = n, num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
                           tau_con = 2*var(y)/m,
                           tau_mod = 0.5*var(y)/m, parallel=F,random_seed=0)
            predict.XBCFrd(fit,w,rep(0,n))
        }
}
##
data <- readstata13::read.dta13("Data/Merged_1853_1864_data.dta")
data$dist_netw <- data$dist_netw/100
data$dist_2 <- ifelse(data$broad==0,-data$dist_netw,data$dist_netw)
vars <- c("deaths","dist_netw","dist_2",
          "dist_cent","dist_urinal",
          "no_sewer","old_sewer")
covs <- c("dist_cent","dist_urinal",
          "no_sewer","old_sewer")
data <- data[,vars]
data <- data[complete.cases(data),]
y <- data$deaths
x <- data$dist_2
w <- data[,covs]
n <- nrow(data)
c <- 0
##
bins <- aggregate(y,list(bins=cut(x,seq(-1,1,0.05),right=F)),mean)
par(mfrow=c(1,2))
plot(seq(-100,95,5),bins$x,type="n",
     xlab="Distance to Broad St. Pump",ylab="Deaths")
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
Omin          <- as.integer(0.03*n)
Opct          <- 0.9
## Owidth        <- quantile(abs(x),0.05)
m             <- 10
num_cutpoints <- n
Nmin          <- 10
num_sweeps    <- 50
burnin        <- 20
p_categorical <- 2
h <- quantile(abs(x),seq(0.05,0.25,0.01))
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
fit <- fit.general(h,y,w,x)
stopImplicitCluster()
##
pred <- sapply(fit,function(i) colMeans(i$tau.adj)[(burnin+1):num_sweeps])
pred <- apply(pred,2,function(i) c(mean(i),quantile(i,c(0.025,0.975))))
pred <- t(pred)
matplot(h,pred,lty=c(1,2,2),pch=19,type=c("b","l","l"),col="blue")
abline(v=quantile(abs(x),0.1),lty=2)
abline(v=quantile(abs(x),0.15),lty=2)
##
tau.pred(pred,burnin,num_sweeps)
print(paste0("MSE: ",round(mean((pred.tot$yhats.adj.mean-y)^2),3)))
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
