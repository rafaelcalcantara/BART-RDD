setwd("~/Git/BART-RDD")
library(doParallel)
set.seed(0)
no_cores <- 10
if (!dir.exists("Data")) dir.create("Data") ## Create data folder
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
### Functions
mu0.x <- function(x,c) 0.1*(x-c+1)^3
mu0.w <- function(w) (rowMeans(w)+2)^2
mu <- function(x,w,k1,k2,k3,c) k1*mu0.x(x,c)+k2*mu0.w(w)*abs(x+1)^k3*sign(x+1)
tau0 <- function(w,ate) pnorm(rowMeans(w)+2,0,1) + pnorm(rowMeans(w)+1,0,1)/4 + dnorm(2*rowMeans(w),0,1) + dnorm(rowMeans(w)-1,0,1)
tau <- function(w,ate,tau.bar) tau0(w,ate) - tau.bar + ate
# h.grid <- function(x,c,grid)
# {
#   abs.x <- sort(abs(x-c))
#   out <- rep(0,length(grid))
#   names(out) <- grid
#   x.right <- sum(c < x)
#   x.left <- sum(x < c)
#   x.tot <- length(x)
#   for(total in grid)
#   {
#     i <- 1
#     sum.right <- sum.left <- 0
#     while(sum.right < total | sum.left < total) 
#     {
#       sum.left <- sum(c-abs.x[i] <= x & x < c)
#       sum.right <- sum(c < x & x <= c+abs.x[i])
#       if (sum.left == sum(x<c) & sum.right == sum(c<x)) break
#       i <- i+1
#     }
#     out[as.character(total)] <- abs.x[i]
#   }
#   return(out)
# }
## Parameters
n <- 5000
rho <- 0.4
k1 <- c(1,5)
k2 <- c(1,5)
k3 <- c(0,0.5)
p <- c(2,5)
sig_error <- c(0.25,0.5)
s <- 1000
c <- 0
ate <- 0.2
ind <- 0
params <- matrix(c(k1[1],k2[1],k3[1],p[1],sig_error[1],
                   k1[2],k2[1],k3[1],p[1],sig_error[1],
                   k1[1],k2[2],k3[1],p[1],sig_error[1],
                   k1[1],k2[1],k3[2],p[1],sig_error[1],
                   k1[1],k2[1],k3[1],p[2],sig_error[1],
                   k1[1],k2[1],k3[1],p[1],sig_error[2]),
                 6,5,byrow=T)
gen.data <- function(ind)
{
  row <- as.numeric(params[ind,])
  K1 <- row[1]
  K2 <- row[2]
  K3 <- row[3]
  P <- row[4]
  sig <- row[5]
  ## Demeaning tau and adjusting parameters
  x0 <- rnorm(n)
  w0 <- matrix(rnorm(n*P,x0,sqrt(1-rho^2)),n,P)
  K2 <- K2*sd(tau0(w0,ate))/sd(mu(c,w0,K1,1,K3,c))
  tau.bar <- mean(tau0(w0,ate))
  sig <- sig*max(abs(mean(tau(w0,ate,tau.bar))),2*sd(tau(w0,ate,tau.bar)))
  ## Sampling
  x <- matrix(rnorm(n*s),n,s)
  z <- apply(x,2,function(i) as.numeric(i>=c))
  w <- apply(x, 2, function(i) matrix(rnorm(n*P,rep(i,P)*rho,sqrt(1-rho^2)),n,P),simplify=F)
  cate <- sapply(1:s, function(i) tau(w[[i]],ate,tau.bar))
  y <- sapply(1:s, function(i) mu(x[,i],w[[i]],K1,K2,K3,c) + tau(w[[i]],ate,tau.bar)*z[,i] + rnorm(n,0,sig))
  # h <- sapply(1:s, function(i) rdrobust::rdbwselect(y[,i],x[,i],c)$bws[2])
  ## Save data
  out <- list(y=y,x=x,z=z,w=w,c=c,tau.x=cate,tau=ate,n=n,rho=rho,k1=K1,k2=K2,k3=K3,p=P,sig_error=sig)
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
