## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  w <- data$w[,i]
  y <- data$y[,i]
  x <- data$x[,i]
  deg <- 5
  if (lvl==1) w2 <- sin(3*pi*w)
  if (lvl==2) w2 <- sin(7*pi*w)*(w-0.5)/(5+exp(-2*w))
  bw <- rdrobust::rdbwselect(y,x,c=c,covs=cbind(w,w2),p=deg,q=deg+1)$bws[4]
  reg <- subset(data.frame(y=y,x=x,w=w,w2=w2,z=data$z[,i]),x>=-bw & x<=bw)
  model <- lm(y~(poly(x,deg,raw=T)+w2)*z+w,data=reg)
  test.sample <- x>=-Owidth & x<=Owidth
  test1 <- data.frame(x=0,w=w,w2=w2,z=1)[test.sample,]
  test0 <- data.frame(x=0,w=w,w2=w2,z=0)[test.sample,]
  tau <- predict(model,newdata=test1)-predict(model,newdata=test0)
  return(tau)
}
##
### BEGIN LOOP
files <- length(list.files("Data"))
s0 <- 1
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  if (paste0("Results/oracle_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/oracle_",i,".rds"))
  } else
  {
    res <- list(results=vector("list",s))
  }
  n <- data$n
  lvl <- data$level
  Owidth <- ifelse(n==500,Ow[2],Ow[1])
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
  saveRDS(list(results=res$results,time=time/(s1-s0+1)),paste0("Results/oracle_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")