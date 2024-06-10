set.seed(7)
setwd("~/../Git/BART-RDD/")
library(doParallel)
library(XBART)
library(lattice)
### Parallelization
no_cores <- detectCores()-1
## Fitting function for prior predictive
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ntrees        <- 5
  Nmin          <- 5
  num_sweeps    <- 120
  burnin        <- 20
  p_categorical <- 2
  y <- ys[,i]
  fit <- XBART::XBCF.rd(y, w, x, c,
                        Owidth = hs, Omin = Om, Opct = Op,
                        num_trees_mod = ntrees,
                        num_trees_con = ntrees,
                        num_cutpoints = obs,
                        num_sweeps = num_sweeps,
                        burnin = burnin, Nmin = Nmin,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(y)/ntrees,
                        tau_mod = 0.5*var(y)/ntrees)
  test <- -hs+c<=x & x<=hs+c
  pred <- XBART::predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
  mean(colMeans(pred$tau.adj[,(burnin+1):num_sweeps]))
}
### Function to plot data
plotit <- function(X,Y,z,c,om) {
  plot(X,Y,col=z+1,pch=19,bty="n",main=bquote(omega==.(om)))
  abline(v=c, lty=2)
}
### Prior predictive
#### Setup
mu.prior <- function(x,w,omega) 0.3*rowSums(w) + (1+2*x-x^2+0.5*x^3)/5+(1-abs(x-c))*sin(omega*x)/10
tau.prior <- function(x,tau.bar) tau.bar - log(x+1)/50
omega <- 5
ate <- 0.5
N <- 1000
c <- 0
s <- 20 ## no of samples of th synthetic DGP
Omin <- c(1,5,10)
Opct <- seq(0.6,0.9,0.1)
h <- seq(0.05,0.2,length=4)
#### Loop
params <- c("N","ATE","Omega","Omin","Opct","h")
for (i in N)
{
  ## Data for calibrating the prior to the simulation DGPs
  obs <- i
  x <- 2*rbeta(i,2,4)-0.75
  z <- as.numeric(x>=c)
  w <- cbind(runif(i,-0.1,0.1),rnorm(i,0,0.2),rbinom(i,1,0.4)-0.4,rbinom(i,1,dnorm(x,c,0.5)))
  w[,4] <- w[,4]-mean(w[,4])
  for (j in 1:length(ate))
  {
    taubar <- ate[j]
    for (k in 1:length(omega))
    {
      ome <- omega[k]
      Ey <- mu.prior(x,w,ome) + tau.prior(x,taubar)*z
      ys <-  Ey + matrix(rnorm(i*s,0,0.5),i,s)
      for (l in 1:length(Omin))
      {
        Om <- Omin[l]
        for (m in 1:length(Opct))
        {
          Op <- Opct[m]
          for (n in 1:length(h))
          {
            hs <- h[n]
            print(paste(params,c(i,taubar,ome,Om,Op,hs),sep=": "))
            ## Fit the model
            cl <- makeCluster(no_cores,type="SOCK")
            registerDoParallel(cl)
            clusterExport(cl,varlist=ls())
            time <- system.time({
              out <- parSapply(cl,1:s,fit)
            })
            stopCluster(cl)
            print(time)
            assign(paste(params,c(i,taubar,ome,Om,Op,hs),sep=".",collapse="_"),out)
          }
        }
      }
    }
  }
}
save.image("prior_predictive.RData")
load("prior_predictive.RData")
## Results
### RMSE
rmse <- matrix(0,length(grep("N.",ls())),7)
colnames(rmse) <- c("N","ATE","Omega","Omin","Opct","h","RMSE")
for(i in 1:length(ate))
{
  objects <- mget(ls())
  sel <- grep(paste0("ATE.",ate[i]),names(objects))
  p <- do.call("rbind",strsplit(names(objects[sel]),"_"))
  p <- apply(p,2, function(i) sub("N.|ATE.|Omega.|Omin.|Opct.|h.","",i))
  temp <- do.call("cbind",objects[sel])
  temp <- sqrt(colMeans((temp-ate[i])^2))
  temp <- cbind(p,temp)
  rmse[(nrow(temp)*(i-1)+1):(nrow(temp)*i),] <- temp
}
rmse <- apply(rmse,2,as.numeric)
rmse <- data.frame(rmse)
### Var
var <- matrix(0,length(grep("N.",ls())),7)
colnames(var) <- c("N","ATE","Omega","Omin","Opct","h","Variance")
for(i in 1:length(ate))
{
  objects <- mget(ls())
  sel <- grep(paste0("ATE.",ate[i]),names(objects))
  p <- do.call("rbind",strsplit(names(objects[sel]),"_"))
  p <- apply(p,2, function(i) sub("N.|ATE.|Omega.|Omin.|Opct.|h.","",i))
  temp <- do.call("cbind",objects[sel])
  temp <- apply(temp,2,var)
  temp <- cbind(p,temp)
  var[(nrow(temp)*(i-1)+1):(nrow(temp)*i),] <- temp
}
var <- apply(var,2,as.numeric)
var <- data.frame(var)
### Plots
boxplot(rmse[,"RMSE"]~rmse[,"h"],xlab="h",ylab="RMSE")
boxplot(var[,"Variance"]~rmse[,"h"],xlab="h",ylab="Variance")
####
plot.res <- function(n,ate,omin,opct)
{
  p <- subset(rmse, N==n & ATE==ate & Omin==omin & Opct==opct, select=c("Omega","h","RMSE"))
  p <- reshape(p,timevar="h",idvar="Omega",direction="wide")
  matplot(omega,p[,-1],type="b",lty=1,pch=19,col=1:length(h),bty="n",ylab="RMSE",
          xlab=bquote(omega), ylim = c(min(p[,-1]),max(p[,-1])+0.1),
          main=paste(params[c(1,2,4,5)],c(n,ate,omin,opct),sep=": ",collapse="; "))
  legend("topleft",legend=h,title="h",lty=1,pch=19,col=1:length(h),ncol=2,cex=0.7)
}
plot.new()
par(mfrow=c(1,2))
plot.res(1000,1,5,0.7)
plot.res(1000,1,5,0.9)
par(mfrow=c(1,1))