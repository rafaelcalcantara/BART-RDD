set.seed(7)
setwd("~/Git/BART-RDD")
library(doParallel)
library(XBART)
library(lattice)
### Parallelization
# no_cores <- detectCores()-1
no_cores <- 125
## Fitting function for prior predictive
fit <- function(i)
{
  print(paste0("Sample: ",i))
  y <- ys[,i]
  y <- y[train]
  x <- x[train]
  w <- w[train]
  fit <- XBART::XBCF.rd(y, w, x, c,
                        Owidth = hs, Omin = Om, Opct = Op,
                        num_cutpoints = n,
                        num_sweeps = num_sweeps,
                        burnin = burnin,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical)
  test <- -hs+c<=x & x<=hs+c
  pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
  pred <- rowMeans(pred$tau.adj[,(burnin+1):num_sweeps])
  sqrt(mean((pred-cate[test])^2))
}
### Prior predictive
#### Setup
mu.prior <- function(x,w) 10*w + 6.5*x^5 - 2.6*x^3 + 1.25*x + 0.5
tau.prior <- function(x,w,c) log(x-c+1) + 2*sin(3*pi*w)-mean(2*sin(3*pi*w))
ate <- 2
c <- 0
s <- no_cores ## no of samples of th synthetic DGP
N <- c(500,1000,2500,5000)
Omin <- c(1,5)
Opct <- c(0.75,0.9)
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
x.list <- list(`500`=NA,`1000`=NA,`2500`=NA,`5000`=NA)
h.list <- list(`500`=NA,`1000`=NA,`2500`=NA,`5000`=NA)
for (n in N) {
  x.list[[as.character(n)]] <- 2*rbeta(n,2,4)-0.75
  h.list[[as.character(n)]] <- h.grid(x.list[[as.character(n)]],c,c(10,30,50))
}
#### Loop
params <- c("N","ATE","Omin","Opct","h")
num_sweeps    <- 150
burnin        <- 50
p_categorical <- 0
for (n in N)
{
  ## Data for calibrating the prior to the simulation DGPs
  x <- x.list[[as.character(n)]]
  z <- as.numeric(x>=c)
  w <- runif(n)
  cate <- tau.prior(c,w,c)
  h <- h.list[[as.character(n)]]
  Ey <- mu.prior(x,w) + z*ate + tau.prior(x,w,c)*z
  ys <-  Ey + matrix(rnorm(n*s,0,sqrt(0.5)),n,s)
  for (Om in Omin)
  {
    for (Op in Opct)
    {
      for (hs in h)
      {
        print(paste(params,c(n,ate,Om,Op,hs),sep=": "))
        train <- c-h.grid(x,c,250) < x & x < c+h.grid(x,c,250)
        ## Fit the model
        cl <- makeCluster(no_cores,type="SOCK")
        registerDoParallel(cl)
        clusterExport(cl,varlist=ls())
        time <- system.time({
          out <- parSapply(cl,1:s,fit)
        })
        stopCluster(cl)
        print(time)
        assign(paste(params,c(n,ate,Om,Op,hs),sep=".",collapse="_"),out)
      }
    }
  }
}
save.image("prior_predictive.RData")
load("prior_predictive.RData")
## Results
### RMSE
objects <- mget(ls())
sel <- grep(paste0("ATE.",ate),names(objects))
p <- do.call("rbind",strsplit(names(objects[sel]),"_"))
p <- apply(p,2, function(i) sub("N.|ATE.|Omin.|Opct.|h.","",i))
rmse <- sapply(objects[sel],mean)
rmse <- cbind(p,RMSE=rmse)
colnames(rmse) <- c(params,"RMSE")
rmse <- apply(rmse,2,as.numeric)
rmse <- data.frame(rmse)
### Plots
# pdf("Figures/prior_predictive.pdf")
par(mfrow=c(2,2),bty="L")
###
for (n in N)
{
  rmse.plot <- reshape(subset(rmse,N==n,select=c("Omin","Opct","h","RMSE")),direction="wide",timevar = "Omin",idvar=c("h","Opct"))
  rmse.plot <- reshape(rmse.plot,direction="wide",timevar = "Opct",idvar="h")
  matplot(x=rmse.plot[,1],rmse.plot[,-1],type="b",lty=2,col=1:length(Omin),
          pch=c(15,16),
          cex.axis=0.75,cex.lab=0.75,cex=0.75,
          ylab="RMSE",xlab="h",main=paste0("N=",n))
  # if (n==1000)
  # {
  #   legend("top",title=expression(N[Omin]),legend=Omin,bty="n",cex=0.65,ncol=2,lty=2,
  #          col=1:length(Omin),lwd=2)
  #   legend("topright",title=bquote(alpha),legend=Opct,bty="n",cex=0.65,pt.cex=0.85,ncol=2,pch=15:17)
  # }
}
# dev.off()
for (n in N) print(subset(rmse[rmse$N==n,],RMSE==min(RMSE)))