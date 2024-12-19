setwd("~/Git/BART-RDD")
library(doParallel)
data <- readRDS("Data/dgp_1.rds")
no_cores <- 11
s <- 1000
n <- 500
c <- 0
ate <- 1
dt <- 2
lvl <- 1
sig <- 0.5
### Functions
mu0.x <- function(x) 3*x^5 - 2.5*x^4 - 1.5*x^3 + 2*x^2 + 3*x + 2
mu0.w <- function(w) 0.5*w^2 - w
tau0.x <- function(x,c) (exp(x)-exp(c))/(1+exp(x))
tau0.w <- function(w,level) {
  if (level==1) out <- sin(3*pi*w)/10 ## 5th degree polynomial on W
  if (level==2) out <- sin(7*pi*w)*(w-0.5)/(5+exp(-2*w)) ## 14th degree polynomial on W
  out <- out-mean(out)
  return(out)
}
mu <- function(x,w) mu0.x(x) + mu0.w(w)
tau <- function(x,c,w,dt,level) tau0.x(x,c) + tau0.w(w,level)*dt
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
load.data <- function(i)
{
  out <- cbind(data$y[,i],data$x[,i],data$w[,i])
}
gen.data <- function(i)
{
  set.seed(i)
  x <- 2*rbeta(n,2,4)-0.75
  h <- h.grid(x,c,10)
  z <- as.numeric(x>=c)
  w <- runif(n)
  cate <- tau(c,c,w,dt,lvl)
  y <- mu(x,w) + tau(x,c,w,dt,lvl)*z + rnorm(n,0,sqrt(sig))
}
###
time.load <- system.time({
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  parLapply(cl,1:s,load.data)
  stopCluster(cl)
})
time.gen <- system.time({
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls()[ls()!="data"])
  parLapply(cl,1:s,gen.data)
  stopCluster(cl)
})
print(time.load)
print(time.gen)