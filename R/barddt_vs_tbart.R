setwd("~/Git/BART-RDD")
library(doParallel)
no_cores <- 10
# Common parameters
p_categorical <- 0
burnin <- 50
num_sweeps <- 150
# BARDDT parameters
Omin <- 1
Opct <- 0.5
ntrees_con <- 30
ntrees_mod <- 30
# T-BART parameters
ntrees <- 30
## Generate data---------------------------------------------------------------
### Functions
mu0.x <- function(x,k) -0.03*x^5 + k*0.5*x^3 + 0.1*x^2 - 0.1*x + 1.72
# mu0.w <- function(w,k) k*2.9*cos(w)
mu0.w <- function(w,k) k*1.4*cos(w)
tau0.x <- function(x,c) 1/(1-exp(-4-x)) - 1.02
tau0.w <- function(w) sin(3*w)
mu <- function(x,w,k1,k2) mu0.x(x,k1) + mu0.w(w,k2)
tau <- function(x,c,w,ate) tau0.x(x,c) + tau0.w(w) + ate
h.grid <- function(x,c,grid)
{
  abs.x <- sort(abs(x-c))
  out <- rep(0,length(grid))
  names(out) <- grid
  x.right <- sum(c < x)
  x.left <- sum(x < c)
  x.tot <- length(x)
  for(total in grid)
  {
    i <- 1
    sum.right <- sum.left <- 0
    while(sum.right < total | sum.left < total) 
    {
      sum.left <- sum(c-abs.x[i] <= x & x < c)
      sum.right <- sum(c < x & x <= c+abs.x[i])
      if (sum.left == sum(x<c) & sum.right == sum(c<x)) break
      i <- i+1
    }
    out[as.character(total)] <- abs.x[i]
  }
  return(out)
}
### Parameters
n <- 1000
p <- 3 # Dim of w
rho <- 0.2
k1 <- 4 ## Makes mu0.x harder for T-BART
k2 <- 1 ## Makes mu0.w larger
pts_in_window <- 50
s <- 25 ## Sim reps
c <- 0
ate <- 1
sig_error <- 0.2
x <- matrix(rnorm(n*s),n,s)
h <- apply(x,2,function(i) h.grid(i,c,pts_in_window))
test <- sapply(1:s, function(i) c-h[i] <= x[,i] & x[,i] <= c+h[i])
z <- apply(x,2,function(i) as.numeric(i>=c))
w <- lapply(1:s, function(i) matrix(rnorm(n*p,rep(x[,i],p)*rho,sqrt(1-rho^2)),n,p))
cate <- apply(sapply(w,rowMeans), 2, function(i) tau(c,c,i,ate))
cate <- sapply(1:s, function(i) cate[test[,i],i])
y <- sapply(1:s, function(i) mu(x[,i],rowMeans(w[[i]]),k1,k2) + tau(x[,i],c,rowMeans(w[[i]]),ate)*z[,i] + rnorm(n,0,sig_error))
## BARDDT fit------------------------------------------------------------------
fit.barddt <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- y[,i]
  ws <- as.matrix(w[[i]])
  xs <- x[,i]
  Owidth <- h[i]
  sample <- h.grid(xs,c,250)
  train <- c-sample < xs & xs < c+sample
  fit <- XBART::XBCF.rd(ys[train], ws[train,], xs[train], c,
                        Owidth = Owidth, Omin = Omin, Opct = Opct,
                        num_cutpoints = n,
                        num_trees_con = ntrees_con, num_trees_mod = ntrees_mod,
                        num_sweeps = num_sweeps,
                        burnin = burnin,
                        Nmin = 1,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = var(ys[train])/ntrees_con, tau_mod = var(ys[train])/ntrees_mod)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  pred <- XBART::predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}
cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time.barddt <- system.time({
  barddt <- parLapply(cl,1:s,fit.barddt)
})
stopCluster(cl)
## T-BART fit------------------------------------------------------------------
fit.tbart <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- y[,i]
  ws <- as.matrix(w[[i]])
  xs <- x[,i]
  zs <- z[,i]
  Owidth <- h[i]
  fit1 <- XBART::XBART(ys[zs==1], cbind(xs,ws)[zs==1,], num_trees = ntrees,
                       num_cutpoints = sum(zs==1), num_sweeps = num_sweeps,
                       burnin = burnin, p_categorical = p_categorical,
                       tau = var(ys[zs==1])/ntrees, parallel=F)
  fit0 <- XBART::XBART(ys[zs==0], cbind(xs,ws)[zs==0,], num_trees = ntrees,
                       num_cutpoints = sum(zs==0), num_sweeps = num_sweeps,
                       burnin = burnin, p_categorical = p_categorical,
                       tau = var(ys[zs==0])/ntrees, parallel=F)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  test.sample <- cbind(c,ws)[test,]
  pred1 <- XBART::predict.XBART(fit1,test.sample)[,(burnin+1):num_sweeps]
  pred0 <- XBART::predict.XBART(fit0,test.sample)[,(burnin+1):num_sweeps]
  pred1-pred0
}
cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time.tbart <- system.time({
  tbart <- parLapply(cl,1:s,fit.tbart)
})
stopCluster(cl)
## Calculate RMSE--------------------------------------------------------------
barddt <- sapply(barddt,rowMeans)
tbart <- sapply(tbart,rowMeans)
rmse.barddt <- mapply(function(fit,tau) sqrt(mean((fit-tau)^2)),barddt,cate)
rmse.tbart <- mapply(function(fit,tau) sqrt(mean((fit-tau)^2)),tbart,cate)
rmse.ate <- sapply(cate, function(tau) sqrt(mean((tau-ate)^2)))
print("RMSE for BARDDT")
summary(rmse.barddt/rmse.ate)
print("RMSE for T-BART")
summary(rmse.tbart/rmse.ate)
plot(cate[[1]],barddt[[1]],bty="n",pch=19,col="darkblue")
abline(a=0,b=1,col="red")
plot(cate[[1]],tbart[[1]],bty="n",pch=19,col="darkblue")
abline(a=0,b=1,col="red")
plot(rmse.barddt/rmse.ate,rmse.tbart/rmse.ate,bty="n",pch=19,col="darkblue")
abline(a=0,b=1,col="red")
matplot(sapply(w,rowMeans)[test[,1],1],cbind(cate[[1]],tbart[[1]],barddt[[1]]),pch=20,col=c("black","maroon","orange"))
points(sapply(w,rowMeans)[test[,1],1],tbart[[1]],pch=20,col='maroon')
points(sapply(w,rowMeans)[test[,1],1],barddt[[1]],pch=20,col='orange')