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
  w1 <- sin(w)
  w2 <- sin(3*pi*w)
  bw <- rdrobust::rdbwselect(y,x,c=c,covs=cbind(w,w2),p=deg.x,q=deg.x+1)$bws[4]
  reg <- subset(data.frame(y=y,x=x,w1=w1,w2=w2,z=data$z[,i]),x>=-bw & x<=bw)
  model <- lm(y~(log(x-c+1)+w2)*z+poly(x,deg.x,raw=T)+w1,data=reg)
  test.sample <- x>=-Owidth & x<=Owidth
  test1 <- data.frame(x=0,w1=w1,w2=w2,z=1)[test.sample,]
  test0 <- data.frame(x=0,w1=w1,w2=w2,z=0)[test.sample,]
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
  if (paste0("oracle_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/oracle_",i,".rds"))
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
  saveRDS(list(results=res$results,time=time),paste0("Results/oracle_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")