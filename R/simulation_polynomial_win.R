## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  # ifelse(is.list(data$w),ws <- data$w[[i]],ws <- subset(data$w,select=i))
  ws <- data$w[,i]
  # ws <- cbind(ws,sapply(3:9,function(i) as.integer(as.integer(10*data$w[,i]) %% i == 0)))\
  y <- data$y[,i]
  x <- data$x[,i]
  bw <- rdrobust::rdbwselect(y,x,c=c,covs=ws,p=7,q=8)$bws[4]
  # reg <- subset(data.frame(y=y,x=x,w1=as.factor(ws[,1]),w2=as.factor(ws[,2]),
  #                          w3=as.factor(ws[,3]),w4=as.factor(ws[,4]),
  #                          w5=as.factor(ws[,5]),w6=as.factor(ws[,6]),
  #                          w7=as.factor(ws[,7]),w8=as.factor(ws[,8]),
  #                          z=data$z[,i]),x>=-bw & x<=bw)
  reg <- subset(data.frame(y=y,x=x,w=ws,z=data$z[,i]),x>=-bw & x<=bw)
  model <- lm(y~((poly(x,6))*(poly(w,4)))*z,data=reg)
  test.sample <- x>=-Owidth & x<=Owidth
  # test1 <- data.frame(x=0,w1=as.factor(ws[,1]),w2=as.factor(ws[,2]),
  #                     w3=as.factor(ws[,3]),w4=as.factor(ws[,4]),
  #                     w5=as.factor(ws[,5]),w6=as.factor(ws[,6]),
  #                     w7=as.factor(ws[,7]),w8=as.factor(ws[,8]),
  #                     z=1)[test.sample,]
  # test0 <- data.frame(x=0,w1=as.factor(ws[,1]),w2=as.factor(ws[,2]),
  #                     w3=as.factor(ws[,3]),w4=as.factor(ws[,4]),
  #                     w5=as.factor(ws[,5]),w6=as.factor(ws[,6]),
  #                     w7=as.factor(ws[,7]),w8=as.factor(ws[,8]),
  #                     z=0)[test.sample,]
  test1 <- data.frame(x=0,w=ws,z=1)[test.sample,]
  test0 <- data.frame(x=0,w=ws,z=0)[test.sample,]
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
  if (paste0("Results/polynomial_",i,".rds") %in% list.files("Results/") == T)
  {
    res <- readRDS(paste0("Results/polynomial_",i,".rds"))
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
  saveRDS(list(results=res$results,time=time/(s1-s0+1)),paste0("Results/polynomial_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")