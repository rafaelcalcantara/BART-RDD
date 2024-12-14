## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- as.matrix(data$w[,i])
  xs <- data$x[,i]
  fit <- XBART::XBCF.rd(ys, ws, xs, c,
                        Owidth = Owidth, Omin = Omin, Opct = Opct,
                        num_cutpoints = n,
                        num_trees_con = 10, num_trees_mod = 5,
                        num_sweeps = num_sweeps,
                        burnin = burnin,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(ys)/10, tau_mod = 0.5*var(ys)/5)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  pred <- XBART::predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}
### BEGIN LOOP
files <- length(list.files("Data"))
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  s <- ncol(data$y)
  if (paste0("bart_rdd_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/bart_rdd_",i,".rds"))
  } else
  {
    res <- list(results=vector("list",s))
  }
  n <- data$n
  c <- data$c
  if (n==500)
  {
    Owidth <- Ow[1]
    Omin <- 1
    Opct <- 0.75
  } else if (n==1000)
  {
    Owidth <- Ow[2]
    Omin <- 1
    Opct <- 0.95
  } else if (n==2500)
  {
    Owidth <- Ow[3]
    Omin <- 1
    Opct <- 0.75
  } else
  {
    Owidth <- Ow[4]
    Omin <- 5
    Opct <- 0.95
  }
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,s0:s1,fit)
  })
  stopCluster(cl)
  print(time)
  res$results[s0:s1] <- out
  saveRDS(list(results=res$results,time=time),paste0("Results/bart_rdd_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")