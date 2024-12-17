set.seed(7)
setwd("../")
library(doParallel)
library(XBART)
library(lattice)
### Parallelization
no_cores <- detectCores()-1
# no_cores <- 125
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
                        num_trees_con = 10,
                        num_trees_mod = 5,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(y)/30,
                        tau_mod = 0.5*var(y)/10)
  test <- -hs+c<=x & x<=hs+c
  pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
  pred <- rowMeans(pred$tau.adj[,(burnin+1):num_sweeps])
  sqrt(mean((pred-cate[test])^2))
}
### Prior predictive
#### Setup
mu.prior <- function(x,w) w + 1/(1+exp(-5*x)) + (1-abs(x-c))*sin(x)/10
tau.prior <- function(x,w,tau.bar) tau.bar - log(x+1)/50 + (w - mean(w))
ate <- 0.5
c <- 0
s <- no_cores*5 ## no of samples of th synthetic DGP
N <- c(500,1000,2500,5000)
Omin <- c(1,3,5)
Opct <- c(0.6,0.75,0.95)
h.grid <- function(x,c,grid)
{
  abs.x <- sort(abs(x))
  out <- rep(0,length(grid))
  names(out) <- grid
  for(total in grid)
  {
    i <- 1
    sum.right <- sum.left <- 0
    while(sum.right < total/2 & sum.left < total/2) 
    {
      sum.left <- sum(c-abs.x[i] < x & x < c)
      sum.right <- sum(c < x & x < c+abs.x[i])
      i <- i+1
    }
    out[as.character(total)] <- abs.x[i]
  }
  return(out)
}
h.list <- list(`500`=h.grid(2*rbeta(500,2,4)-0.75,c,2:12*10),
               `1000`=h.grid(2*rbeta(1000,2,4)-0.75,c,2:12*10),
               `2500`=h.grid(2*rbeta(2500,2,4)-0.75,c,2:12*10),
               `5000`=h.grid(2*rbeta(5000,2,4)-0.75,c,2:12*10))
#### Loop
params <- c("N","ATE","Omin","Opct","h")
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 0
for (n in N)
{
  ## Data for calibrating the prior to the simulation DGPs
  x <- 2*rbeta(n,2,4)-0.75
  z <- as.numeric(x>=c)
  w <- runif(n)
  cate <- tau.prior(c,w,ate)
  h <- h.list[[as.character(n)]]
  Ey <- mu.prior(x,w) + tau.prior(x,w,ate)*z
  ys <-  Ey + matrix(rnorm(n*s,0,0.1*sd(Ey)),n,s)
  for (Om in Omin)
  {
    for (Op in Opct)
    {
      for (hs in h)
      {
        print(paste(params,c(n,ate,Om,Op,hs),sep=": "))
        train <- -5*hs < x & x < 5*hs
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
  matplot(x=rmse.plot[,1],rmse.plot[,-1],type="b",lty=2,col=rep(1:length(Omin),length(Opct)),
          pch=c(rep(15,length(Opct)),rep(16,length(Opct)),rep(17,length(Opct))),
          cex.axis=0.75,cex.lab=0.75,cex=0.75,
          ylab="RMSE",xlab="h",main=paste0("N=",n))
  if (n==1000)
  {
    legend("top",title=expression(N[Omin]),legend=Omin,bty="n",cex=0.65,ncol=2,lty=2,
           col=1:length(Omin),lwd=2)
    legend("topright",title=bquote(alpha),legend=Opct,bty="n",cex=0.65,pt.cex=0.85,ncol=2,pch=15:17)
  }
}
# dev.off()
for (n in N) print(subset(rmse[rmse$N==n,],RMSE==min(RMSE)))