k1 <- c(1,5)
k2 <- c(1,0.25)
k3 <- c(0,1)
k4 <- c(0.1,0.5)
k5 <- c(0,1)
p <- c(2,4)
rho <- c(0,0.5)
dgp <- matrix(0,14,7)
dgp <- matrix(1,7,7)+diag(1,7)
dgp <- rbind(dgp,matrix(2,7,7)+diag(-1,7))
for (i in 1:nrow(dgp))
{
  config <- dgp[i,]
  dgp[i,] <- c(k1[config[1]],k2[config[2]],k3[config[3]],k4[config[4]],
               k5[config[5]],p[config[7]],rho[config[6]])
}
colnames(dgp) <- c("k1","k2","k3","k4","k5","p","rho")
models <- c("leaf.rdd","tbart","sbart","polynomial")
##
txt <- NA
for (i in 1:nrow(dgp))
{
  dgp.id <- paste(colnames(dgp),dgp[i,],collapse="_",sep="_")
  txt[1] <- "library(doParallel)"
  txt[2] <- "no_cores <- 10"
  txt[3] <- "c <- 0"
  txt[4] <- "Owidth <- 0.1"
  txt[5] <- "## DGP parameters"
  txt[6] <- paste0("k1 <- ", dgp[i,1], " # larger gives more weight to x relative to w in mu(x,w)")
  txt[7] <- paste0("k2 <- ", dgp[i,2], " # relative size of sd(tau(w, x = c)) to sd(mu(w, x =c))")
  txt[8] <- paste0("k3 <- ", dgp[i,3], " # separable in x and w or not, 0 is separable")
  txt[9] <- paste0("k4 <- ", dgp[i,4], " # noise relative to sd(mu)")
  txt[10] <- paste0("k5 <- ", dgp[i,5], " # size of the ate as a multiple of the residual standard deviation")
  txt[11] <- paste0("p <- ", dgp[i,6])
  txt[12] <- paste0("rho <- ", dgp[i,7])
  txt[13] <- "## Sample size"
  txt[14] <- paste0("n <- ",n)
  txt[15] <- "## Simulation reps"
  txt[16] <- paste0("s <- ",s)
  txt[17] <- "## Which models to run"
  txt[18] <- paste0("models <- c('",paste(models,collapse="','"),"')")
  txt[19] <- "## Generate data"
  txt[20] <- "### Identifier of DGP configuration for names of data and results files"
  txt[21] <- paste0("dgp <- '", dgp.id,"'")
  txt[22] <- '## Generate the fixed w used in all settings with p features'
  txt[23] <- 'set.seed(007)'
  txt[24] <- 'K <- 2*toeplitz(seq(1,0,length.out = p))'
  txt[25] <- 'w <- MASS::mvrnorm(n,rep(0,p),K)'
  txt[26] <- 'write.table(w, paste0("Data/w_",p,".csv"), row.names = FALSE, col.names = FALSE, sep = ",")'
  txt[27] <- '## Create log folders for current DGP'
  txt[28] <- 'if (dgp %in% list.files("Logs") == FALSE) dir.create(paste0("Logs/",dgp))'
  txt[29] <- "## Estimation"
  if (isFALSE(cluster)){
    txt[30] <- 'source("simulation_estimation_local.R")'
  } else
  {
    txt[30] <- 'batch <- c(1,1)'
    txt[31] <- 'for (i in 0:(s-1)) { # Running simulations in batches of 10' 
    txt[32] <- 's0 <- batch[1]+i*1'
    txt[33] <- 's1 <- batch[2]+i*1'
    txt[34] <- 'batch.args <- paste(s0,s1,dgp,n,c,Owidth,paste(k1,k2,k3,k4,k5,p,rho, colapse=" "),collapse=" ")'
    txt[35] <- 'batch.script <- paste0("nice Rscript --verbose simulation_estimation_cluster.R ",batch.args, " > Logs/", dgp, "/outputFile_batch", i+1, "_", dgp, ".Rout 2> Logs/", dgp, "/errorFile_batch", i+1, "_", dgp, ".Rout"," &")'
    txt[36] <- 'system(batch.script)'
    txt[37] <- '}'
  }
  writeLines(txt,paste0("simulation_master",i,".R"))
}