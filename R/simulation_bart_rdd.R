## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- as.matrix(data$w[,i])
  xs <- data$x[,i]
  train <- c-h.grid(xs,c,500) < xs & xs < c+h.grid(xs,c,500)
  fit <- XBART::XBCF.rd(ys[train], ws[train], xs[train], c,
                        Owidth = Owidth, Omin = Omin, Opct = Opct,
                        num_cutpoints = n,
                        num_trees_con = 10, num_trees_mod = 5,
                        num_sweeps = num_sweeps,
                        burnin = burnin,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(ys[train])/10, tau_mod = 0.5*var(ys[train])/5)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  pred <- XBART::predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}
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
    while(sum.right < total/2 | sum.left < total/2) 
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
### BEGIN LOOP
for (i in files)
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
    Opct <- 0.9
  } else if (n==1000)
  {
    Owidth <- Ow[2]
    Omin <- 3
    Opct <- 0.9
  } else if (n==2500)
  {
    Owidth <- Ow[3]
    Omin <- 3
    Opct <- 0.9
  } else
  {
    Owidth <- Ow[4]
    Omin <- 3
    Opct <- 0.9
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