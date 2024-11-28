## Setup
set.seed(0)
### Parameters
Omin          <- 5
Opct          <- 0.9
ntrees_con    <- 10
ntrees_mod    <- 5
Nmin          <- 5
num_sweeps    <- 150
burnin        <- 50
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- as.matrix(data$w[,i])
  xs <- data$x[,i]
  fit <- XBART::XBCF.rd(ys, ws, xs, c,
                        Owidth = Owidth, Omin = Omin, Opct = Opct,
                        num_trees_mod = ntrees_mod,
                        num_trees_con = ntrees_con,
                        num_cutpoints = n,
                        num_sweeps = num_sweeps,
                        burnin = burnin, Nmin = Nmin,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(ys)/ntrees_con, max_depth = max_depth,
                        tau_mod = 0.5*var(ys)/ntrees_mod)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  pred <- XBART::predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}
##
### BEGIN LOOP
files <- length(list.files("Data"))
s0 <- 1
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  if (paste0("Results/bart_rdd_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/bart_rdd_",i,".rds"))
  } else
  {
    res <- list(results=vector("list",s))
  }
  n <- nrow(data$y)
  s <- ncol(data$y)
  s1 <- s
  c <- data$c
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,s0:s1,fit)
  })
  stopCluster(cl)
  print(time)
  res$results[s0:s1] <- out
  saveRDS(list(results=res$results,time=time/(s1-s0+1)),paste0("Results/bart_rdd_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")