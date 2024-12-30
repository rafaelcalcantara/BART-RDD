setwd("~/Git/BART-RDD")
library(doParallel)
set.seed(0)
no_cores <- 11
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
### Functions
mu0.x <- function(x) 0.015*x^5 - 0.006*x^3 + 0.0025*x + 3.5
mu0.w <- function(w) -3*sin(w)
tau0.x <- function(x,c) x
tau0.w <- function(w) sin(w)
mu <- function(x,w) {
  mu.w <- mu0.w(w)
  mu0.x(x) + mu.w
}
tau <- function(x,c,w,ate) {
  tau.w <- tau0.w(w)
  tau0.x(x,c) + tau.w + ate
}
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
## Parameters
N <- c(500,1000,1500)
rho <- c(0.5,0.9)
pts_in_window <- 75
s <- 1000
c <- 0
ate <- 1
sig_error <- 1
ind <- 0
params <- expand.grid(N,rho)
gen.data <- function(ind)
{
  row <- as.numeric(params[ind,])
  n <- row[1]
  Rho <- row[2]
  # u1 <- rnorm(n*s)
  # u2 <- rnorm(n*s,Rho*u1,sqrt(1-Rho^2))
  # u <- pnorm(cbind(u1,u2))
  # x <- matrix(2*qbeta(u[,1],2,4)-1,n,s)
  # x <- matrix(qunif(u[,1],-1,0.2),n,s)
  x <- matrix(rnorm(n*s),n,s)
  h <- apply(x,2,function(i) h.grid(i,c,pts_in_window))
  z <- apply(x,2,function(i) as.numeric(i>=c))
  w <- apply(x,2,function(i) rnorm(n,Rho*i,sqrt(1-Rho^2)))
  # w <- matrix(2*qbeta(u[,2],2,5)-1,n,s)
  # w <- matrix(qunif(u[,2],0.1,1),n,s)
  cate <- apply(w, 2, function(i) tau(c,c,i,ate))
  y <- sapply(1:s, function(i) mu(x[,i],w[,i]) + tau(x[,i],c,w[,i],ate)*z[,i] + rnorm(n,0,sig_error))
  ## Save data
  out <- list(y=y,x=x,z=z,w=w,c=c,h=h,tau.x=cate,tau=ate,n=n,rho=Rho,sig_error=sig_error)
  saveRDS(out,paste0("Data/dgp_",ind,".rds"))
}
##
cl <- makeCluster(no_cores,type="SOCK")
registerDoParallel(cl)
clusterExport(cl,varlist=ls())
time <- system.time({
  out <- parLapply(cl,1:nrow(params),gen.data)
})
stopCluster(cl)
print(time)
print("Data generated")
