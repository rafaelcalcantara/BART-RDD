library(doParallel)
no_cores <- 10
c <- 0
Owidth <- 0.1
## DGP parameters
k1 <- 2 # larger gives more weight to x relative to w in mu(x,w)
k2 <- 2 # relative size of sd(tau(w, x = c)) to sd(mu(w, x =c))
k3 <- 0 # separable in x and w or not, 0 is separable
k4 <- 0.1 # noise relative to sd(mu)
k5 <- 0.25 # size of the ate as a multiple of the residual standard deviation
p <- 2
rho <- 0
## Sample size
n <- 4000
## Simulation reps
s <- 200
## Which models to run
models <- c('leaf.rdd','tbart','sbart','polynomial')
## Generate data
### Identifier of DGP configuration for names of data and results files
dgp <- 'k1_2_k2_2_k3_0_k4_0.1_k5_0.25_p_2_rho_0'
source("simulation_data.R")
## Estimation
batch <- c(1,10)
for (i in 0:19) { # Running simulations in batches of 10
s0 <- batch[1]+i*10
s1 <- batch[2]+i*10
batch.args <- paste(s0,s1,dgp,n,c,Owidth,collapse=" ")
batch.script <- paste0("nice Rscript --verbose simulation_estimation_cluster.R ",batch.args, " > Logs/", dgp, "/outputFile_batch", i+1, "_", dgp, ".Rout 2> Logs/", dgp, "/errorFile_batch", i+1, "_", dgp, ".Rout"," &")
system(batch.script)
}
