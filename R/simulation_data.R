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
## Parameters
s <- 1000
c <- 0
ate <- 1
delta_tau <- c(2,4,6)
level <- 1
N <- c(500,1000,2500,5000)
sig_error <- 0.5
ind <- 0
params <- expand.grid(delta_tau,level,N,sig_error)
gen.data <- function(ind)
{
  row <- as.numeric(params[ind,])
  dt <- row[1]
  lvl <- row[2]
  n <- row[3]
  sig <- row[4]
  x <- matrix(2*rbeta(n*s,2,4)-0.75,n,s)
  h <- apply(x,2,function(i) h.grid(i,c,10))
  z <- apply(x,2,function(i) as.numeric(i>=c))
  w <- matrix(runif(n*s),n,s)
  cate <- apply(w, 2, function(i) tau(c,c,i,dt,lvl))
  y <- sapply(1:s, function(i) mu(x[,i],w[,i]) + tau(x[,i],c,w[,i],dt,lvl)*z[,i] + rnorm(n,0,sqrt(sig)))
  ## Save data
  out <- list(y=y,x=x,z=z,w=w,c=c,h=h,tau.x=cate,tau=ate,delta_tau=dt,level=lvl,n=n,sig_error=sig)
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