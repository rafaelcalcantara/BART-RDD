library(doParallel)
no_cores <- 10
c <- 0
Owidth <- 0.1
## DGP parameters
k1 <- 1 # larger gives more weight to x relative to w in mu(x,w)
k2 <- 0.25 # relative size of sd(tau(w, x = c)) to sd(mu(w, x =c))
k3 <- 0 # separable in x and w or not, 0 is separable
k4 <- 0.1 # noise relative to sd(mu)
k5 <- 0 # size of the ate as a multiple of the residual standard deviation
p <- 2
rho <- 0
## Sample size
n <- 4000
## Simulation reps
s <- 100
## Which models to run
models <- c('leaf.rdd','tbart','sbart','polynomial')
## Generate data
### Identifier of DGP configuration for names of data and results files
dgp <- 'k1_1_k2_0.25_k3_0_k4_0.1_k5_0_p_2_rho_0'
## Generate the fixed w used in all settings with p features
set.seed(007)
K <- 2*toeplitz(seq(1,0,length.out = p))
w <- MASS::mvrnorm(n,rep(0,p),K)
write.table(w, paste0("Data/w_",p,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")
## Create log folders for current DGP
if (dgp %in% list.files("Logs") == FALSE) dir.create(paste0("Logs/",dgp))
## Estimation
batch <- c(1,1)
for (i in 0:99) { # Running simulations in batches of 10
s0 <- batch[1]+i*1
s1 <- batch[2]+i*1
batch.args <- paste(s0,s1,dgp,n,c,Owidth,paste(k1,k2,k3,k4,k5,p,rho, colapse=" "),collapse=" ")
batch.script <- paste0("nice Rscript --verbose simulation_estimation_cluster.R ",batch.args, " > Logs/", dgp, "/outputFile_batch", i+1, "_", dgp, ".Rout 2> Logs/", dgp, "/errorFile_batch", i+1, "_", dgp, ".Rout"," &")
system(batch.script)
}
