## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ws <- data$w[,i]
  y <- data$y[,i]
  x <- data$x[,i]
  bw <- rdrobust::rdbwselect(y,x,c=c,covs=ws,p=5,q=6)$bws[4]
  reg <- subset(data.frame(y=y,x=x,w=ws,w2=sin(15*pi*ws)*(ws-mean(ws))/sqrt(ws) +ws*x+0.2*(x-mean(x))/sd(x),z=data$z[,i]),x>=-bw & x<=bw)
  model <- lm(y~(poly(x,5)+w2)*z+w,data=reg)
  test.sample <- x>=-Owidth & x<=Owidth
  test1 <- data.frame(x=0,w=ws,w2=sin(15*pi*ws)*(1.2*ws-mean(ws)),z=1)[test.sample,]
  test0 <- data.frame(x=0,w=ws,w2=sin(15*pi*ws)*(1.2*ws-mean(ws)),z=0)[test.sample,]
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
  saveRDS(list(results=res$results,time=time/(s1-s0+1)),paste0("Results/oracle_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")