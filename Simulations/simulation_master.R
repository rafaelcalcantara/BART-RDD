library(doParallel)
no_cores <- 10
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0)
{
  ## Running from command line
  cmd.line <- TRUE
  ## DGP parameters
  k1 <- as.numeric(args[1]) # larger gives more weight to x relative to w in mu(x,w)
  k2 <- as.numeric(args[2]) # relative size of sd(tau(w, x = c)) to sd(mu(w, x =c))
  k3 <- as.numeric(args[3]) # separable in x and w or not, 0 is separable
  k4 <- as.numeric(args[4]) # noise relative to sd(mu)
  k5 <- as.numeric(args[5]) # size of the ate as a multiple of the residual standard deviation
  p <- as.numeric(args[6])
  rho <- as.numeric(args[7])
  ## Sample size
  n <- as.numeric(args[8])
  ## Simulation reps
  s <- as.numeric(args[9])
  ## Which models to run
  models <- which(args %in% c("leaf.rdd","tbart","sbart","polynomial"))
  models <- args[models]
} else
{
  ## Running from R
  cmd.line <- FALSE
  ## DGP parameters
  k1 <- 2 # larger gives more weight to x relative to w in mu(x,w)
  k2 <- 0.25 # relative size of sd(tau(w, x = c)) to sd(mu(w, x =c))
  k3 <- 1 # separable in x and w or not, 0 is separable
  k4 <- 0.25 # noise relative to sd(mu)
  k5 <- 1 # size of the ate as a multiple of the residual standard deviation
  p <- 2
  rho <- 0.25
  ## Sample size
  n <- 4000
  ## Simulation reps
  s <- 1
  ## Which models to run
  models <- c("leaf.rdd","tbart","sbart","polynomial")
  #models <- c("polynomial")
}
### Identifier of DGP configuration for names of data and results files
dgp <- paste(c("k1","k2","k3","k4","k5","p","rho"),c(k1,k2,k3,k4,k5,p,rho),collapse="_",sep="_")
print(dgp)
print(class(p))
## Generate data
source("simulation_data.R")
## Estimation
Owidth <- 0.1
source("simulation_estimation.R")