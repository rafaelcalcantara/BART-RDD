## Setup
set.seed(0)
### Functions
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
    while(sum.right < total | sum.left < total) 
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
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys <- data$y[,i]
  ws <- as.matrix(data$w[,i])
  xs <- data$x[,i]
  Owidth <- data$h[i]
  sample <- h.grid(xs,c,1250)
  train <- c-sample < xs & xs < c+sample
  fit <- XBART::XBCF.rd(ys[train], ws[train], xs[train], c,
                        Owidth = Owidth, Omin = Omin, Opct = Opct,
                        num_cutpoints = n,
                        num_trees_con = 20, num_trees_mod = 10,
                        num_sweeps = 250,
                        burnin = 50,
                        Nmin = 1,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(ys[train])/10, tau_mod = 0.5*var(ys[train])/5)
  test <- -Owidth+c<=xs & xs<=Owidth+c
  pred <- XBART::predict.XBCFrd(fit,ws[test,],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}
### Parameters
Omin <- 5
Opct <- 0.75
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
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,s0:s1,fit)
  })
  stopCluster(cl)
  print(time)
  res$results[s0:s1] <- out
  ### DEBUG
  test <- sapply(1:10, function(i) -data$h[i]<=data$x[,i] & data$x[,i] <= data$h[i])
  cate <- unlist(sapply(1:10,function(i) data$tau.x[test[,i],i]))
  w <- unlist(sapply(1:10,function(i) data$w[test[,i],i]))
  pred <- unlist(sapply(out[1:10],rowMeans))
  matplot(w,cbind(pred,cate))
  ###
  saveRDS(list(results=res$results,time=time),paste0("Results/bart_rdd_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")