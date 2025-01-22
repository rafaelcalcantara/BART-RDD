## Setup
set.seed(0)
### Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  z <- data$z[,i]
  w <- data$w[[i]]
  y <- data$y[,i]
  x <- data$x[,i]
  Owidth <- data$h[i]
  poly.data <- data.frame(y=y,x=x,w1=w[,1],w2=w[,2],z=z)
  poly.data <- subset(poly.data, c-Owidth <= x & x <= c+Owidth)
  deg.x <- 3
  deg.w.mu <- 2
  deg.w.tau <- 3
  poly.fit <- lm(y~(poly(x,degree=deg.x,raw=T)*(poly(w1,degree=deg.w.mu,raw=T)*poly(w2,degree=deg.w.mu,raw=T)))+(poly(w1,degree=deg.w.tau,raw=T)*poly(w2,degree=deg.w.tau,raw=T))*z,data=poly.data)
  test.sample <- x>=-Owidth & x<=Owidth
  test1 <- data.frame(x=0,w1=w[,1],w2=w[,2],z=1)[test.sample,]
  test0 <- data.frame(x=0,w1=w[,1],w2=w[,2],z=0)[test.sample,]
  tau <- predict(poly.fit,newdata=test1)-predict(poly.fit,newdata=test0)
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