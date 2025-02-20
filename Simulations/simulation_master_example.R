library(doParallel)
no_cores <- 10
c <- 0
Owidth <- 0.1
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
s <- 5
## Which models to run
models <- c("leaf.rdd","tbart","sbart","polynomial")
## Generate data
### Identifier of DGP configuration for names of data and results files
dgp <- paste(c("k1","k2","k3","k4","k5","p","rho"),c(k1,k2,k3,k4,k5,p,rho),collapse="_",sep="_")
source("simulation_data.R")
## Estimation
source("simulation_estimation.R")