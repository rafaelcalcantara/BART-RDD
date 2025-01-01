## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  w <- data$w[,i]
  y <- data$y[,i]
  x <- data$x[,i]
  Owidth <- data$h[i]
  deg.x <- 5
  deg.w <- 4
  bw <- rdrobust::rdbwselect(y,x,c=c,covs=w,p=deg.x,q=deg.x+1)$bws[4]
  reg <- subset(data.frame(y=y,x=x,w=w,z=data$z[,i]),x>=-bw & x<=bw)
  model <- lm(y~(poly(x,deg.x,raw=T)+poly(w,deg.w,raw=T))*z,data=reg)
  test.sample <- x>=-Owidth & x<=Owidth
  test1 <- data.frame(x=0,w=w,z=1)[test.sample,]
  test0 <- data.frame(x=0,w=w,z=0)[test.sample,]
  tau <- predict(model,newdata=test1)-predict(model,newdata=test0)
  return(tau)
}
##
### BEGIN LOOP
for (i in files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  s <- ncol(data$y)
  if (paste0("polynomial_full_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/polynomial_full_",i,".rds"))
  } else
  {
    res <- list(results=vector("list",s))
  }
  n <- data$n
  lvl <- data$level
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
  saveRDS(list(results=res$results,time=time),paste0("Results/polynomial_full_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")