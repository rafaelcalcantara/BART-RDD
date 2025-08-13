library(doParallel)
no_cores <- 10
c <- 0
Owidth <- 0.1
## DGP parameters
k1 <- 1 # larger gives more weight to x relative to w in mu(x,w)
k2 <- 1 # relative size of sd(tau(w, x = c)) to sd(mu(w, x =c))
k3 <- 0 # separable in x and w or not, 0 is separable
k4 <- 0.1 # noise relative to sd(mu)
k5 <- 0 # size of the ate as a multiple of the residual standard deviation
p <- 4
rho <- 0
## Sample size
n <- 4000
## Simulation reps
s <- 100
## Which models to run
models <- c('horseshoe','ridge')
## Generate data
### Identifier of DGP configuration for names of data and results files
dgp <- 'k1_1_k2_1_k3_0_k4_0.1_k5_0_p_4_rho_0'
## Generate the fixed w used in all settings with p features
set.seed(007)
K <- 2*toeplitz(seq(1,0,length.out = p))
w <- MASS::mvrnorm(n,rep(0,p),K)
saveRDS(w, paste0("Data/w_",p,".rds"))
## Create log folders for current DGP
if (dgp %in% list.files("Logs") == FALSE) dir.create(paste0("Logs/",dgp))
## Estimation
source("simulation_estimation_local.R")
