# setwd("~/Downloads/")
## DGP parameters
# n <- 2000
# ate <- 1
# k1 <- 1 ## variability in mu0.x
# k2 <- 1 ## amplitude of tau0.w relative to ATE (set from 0 to 1)
# k3 <- 1 ## sd of mu0.w relative to sd of tau
# sig_error <- 0.2 ## relative to sd of tau
# p <- 2 # Dim of w
# rho <- 0
# c <- 0
# s <- 2 ## Sim reps
# pts_in_window <- 250
# x.center <- 0 ## In case we want to center x at a different location
print.params <- T ## sanity check
## DGP functions
### Functions


mu0.x <- function(x,k,ate) -0.03*x^4 + 0.5*x^3 + 0.5*x^2 - 0.1*x + 1
mu0.w <- function(w,k) k*cos(rowMeans(w))
tau0.x <- function(x,c) rep(1,length(x))
tau0.w <- function(w,k) k*sin(6*rowMeans(w)-1)
mu0 <- function(x,w,k1,k3,ate) k1*mu0.x(x,k1,ate) + k3*rowMeans(w)
tau0 <- function(x,c,w,ate,k2) tau0.w(w,k2) + tau0.x(x,c)
mu <- function(x,w,k1,k3,ate,mu.bar) mu0(x,w,k1,k3,ate) - mu.bar 
tau <- function(x,c,w,ate,k2,tau.bar) tau0(x,c,w,ate,k2) - tau.bar + ate
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
## Generate data
source("simulation_generate_data.R")
## Run simulations
source("simulation_estimation_par.R")