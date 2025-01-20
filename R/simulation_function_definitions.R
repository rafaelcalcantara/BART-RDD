setwd("~/Git/BART-RDD")
## DGP parameters
n <- 5000
ate <- 0.5
k1 <- 4 ## variability in mu0.x
k2 <- 0.75 ## amplitude of tau0.w relative to ATE (set from 0 to 1)
k3 <- 2 ## sd of mu0.w relative to sd of tau
sig_error <- 0.5 ## relative to sd of tau
p <- 2 # Dim of w
rho <- 0.5
c <- 1
s <- 10 ## Sim reps
pts_in_window <- 50
x.center <- 0 ## In case we want to center x at a different location
print.params <- T ## sanity check
sample <- 1 ## Sample for plotting posterior and monitoring constraints
## DGP functions
### Functions
mu0.x <- function(x,k,ate)
{
  -0.03*x^5 + k*0.5*x^3 + 0.1*x^2 - 0.1*x + 2
}
mu0.w <- function(w,k) 
{
  k*cos(rowMeans(w))
}
tau0.x <- function(x)
{
  0.25/(3+x-min(x))
}
tau0.w <- function(w,k)
{
  k*sin(3*rowMeans(w))
}
mu0 <- function(x,w,k1,k3,ate)
{
  mu0.x(x,k1,ate)*rowMeans(w)
}
tau0 <- function(x,w,ate,k2)
{
  tau0.x(x) + tau0.w(w,k2)
}
mu <- function(x,w,k1,k3,ate,mu.mean)
{
  mu0(x,w,k1,k3,ate) - mu.mean
}
tau <- function(x,w,ate,k2,tau.mean)
{
  tau0(x,w,ate,k2) - tau.mean + ate
}
h.grid <- function(x,c,grid)
{
  abs.x <- sort(abs(x-c))
  out <- rep(0,length(grid))
  names(out) <- grid
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
## Generate data
source("R/simulation_generate_data.R")
## Run simulations
source("R/simulation_estimation.R")
## Monitor split rejection and forced splits
source("R/simulation_monitor_splits.R")