## Setup
set.seed(0)
# setwd("../")
if (!dir.exists("Results")) dir.create("Results") ## Create results folder
if (length(list.files("Results")[grep("linear_",list.files("Results"))])!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results")[grep("linear_",list.files("Results"))])
  for (i in files) file.remove(i)
}
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  # ifelse(is.list(data$w),ws <- data$w[[i]],ws <- subset(data$w,select=i))
  ws <- data$ws[,i]
  ws <- cbind(ws,sapply(3:9,function(i) as.integer(as.integer(10*data$w[,i]) %% i == 0)))
  y <- data$y[,i]
  x <- data$x[,i]
  # ws <- as.integer(cut(ws,quantile(ws,seq(0,1,length.out=classes)),include.lowest=T)) %% 2 == 0
  bw <- rdrobust::rdbwselect(y,x,c=c,covs=ws,p=7,q=8)$bws[4]
  reg <- subset(data.frame(y=y,x=x,w1=as.factor(ws[,1]),w2=as.factor(ws[,2]),
                           w3=as.factor(ws[,3]),w4=as.factor(ws[,4]),
                           w5=as.factor(ws[,5]),w6=as.factor(ws[,6]),
                           w7=as.factor(ws[,7]),w8=as.factor(ws[,8]),
                           z=data$z[,i]),x>=-bw & x<=bw)
  model <- lm(y~(x^7+x^6+x^5+x^4+x^3+x^2+x*(w1+w2+w3)+w4+w5+w6+w7+w8)*z,data=reg)
  test.sample <- x>=-Owidth & x<=Owidth
  test1 <- data.frame(x=0,w1=as.factor(ws[,1]),w2=as.factor(ws[,2]),
                      w3=as.factor(ws[,3]),w4=as.factor(ws[,4]),
                      w5=as.factor(ws[,5]),w6=as.factor(ws[,6]),
                      w7=as.factor(ws[,7]),w8=as.factor(ws[,8]),
                      z=1)[test.sample,]
  test0 <- data.frame(x=0,w1=as.factor(ws[,1]),w2=as.factor(ws[,2]),
                      w3=as.factor(ws[,3]),w4=as.factor(ws[,4]),
                      w5=as.factor(ws[,5]),w6=as.factor(ws[,6]),
                      w7=as.factor(ws[,7]),w8=as.factor(ws[,8]),
                      z=0)[test.sample,]
  tau <- predict(model,newdata=test1)-predict(model,newdata=test0)
  return(tau)
}
##
### BEGIN LOOP
files <- length(list.files("Data"))
for (i in 1:files)
{
  print(paste0("DGP: ",i))
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  n <- nrow(data$y)
  s <- ncol(data$y)
  c <- data$c
  classes <- data$classes
  # test <- readRDS(paste0("Data/test_dgp_",i,".rds"))
  # test.w <- test$w
  # test.sample <- cbind(c,test.w)
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s,fit)
  })
  stopCluster(cl)
  print(time)
  saveRDS(list(results=out,time=time),paste0("Results/polynomial_",i,".rds"))
  rm(out)
  gc()
}
##
print("Done!")