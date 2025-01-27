library(doParallel)
c <- 0
Owidth <- 0.1
args <- commandArgs(trailingOnly = TRUE)
## Choose to run simulations or process results
### Toggle on if just processing results from cluster (don't run data and estimation scripts)
if (length(args) > 0)
{
  ## If run from cmd line, only run simulations
  results <- FALSE
} else
{
  ## Choose whether to run simulations or process results when sourcing script from R
  results <- FALSE
}
## Run simulations/process results
if (isFALSE(results))
{
  if (length(args) > 0)
  {
    ## Running from command line
    cmd.line <- TRUE
    no_cores <- 5
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
    ### Identifier of DGP configuration for names of data and results files
    dgp <- paste(c("k1","k2","k3","k4","k5","p","rho"),c(k1,k2,k3,k4,k5,p,rho),collapse="_",sep="_")
    list2env(readRDS(paste0("Data/dgp_",dgp,".rds")),envir=.GlobalEnv)
    ## Estimation
    source("simulation_estimation_local.R")
  } else
  {
    ## Running from R
    cmd.line <- FALSE
    no_cores <- 4
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
    s <- 500
    ## Which models to run
    models <- c("leaf.rdd","tbart","sbart","polynomial")
    ## Generate data
    ### Identifier of DGP configuration for names of data and results files
    dgp <- paste(c("k1","k2","k3","k4","k5","p","rho"),c(k1,k2,k3,k4,k5,p,rho),collapse="_",sep="_")
    source("simulation_data.R")
    ## Estimation
    source("simulation_estimation_local.R")
  }
} else
{
  ## Process results
  source("simulation_results.R")
}