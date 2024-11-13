set.seed(7)
setwd("../")
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
  # pred <- XBART::predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
  pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
  mean(colMeans(pred$tau.adj[,(burnin+1):num_sweeps]))
}
### Prior predictive
#### Setup
mu.prior <- function(x,w,omega) w + 1/(1+exp(-5*x)) + (1-abs(x-c))*sin(omega*x)/10
tau.prior <- function(x,tau.bar) tau.bar - log(x+1)/50
omega <- 0
ate <- 0.4
N <- 1000
c <- 0
s <- 20 ## no of samples of th synthetic DGP
Omin <- c(5,10,20)
Opct <- seq(0.6,0.9,length=3)
h <- seq(0.05,0.2,length=4)
classes <- 5
p <- 0.4
#### Loop
params <- c("N","ATE","Omega","Omin","Opct","h")
for (i in N)
{
  ## Data for calibrating the prior to the simulation DGPs
  obs <- i
  x <- 2*rbeta(i,2,4)-0.75
  z <- as.numeric(x>=c)
  w <- rbinom(i,classes,p)+1
  ## plot(x,mu.prior(x,w,omega[3])+tau.prior(x,0.25)*z+rnorm(N),col=z+1,pch=19)
  for (j in 1:length(ate))
  {
    taubar <- ate[j]
    for (k in 1:length(omega))
    {
      ome <- omega[k]
      Ey <- mu.prior(x,w,ome) + tau.prior(x,taubar)*z
      ys <-  Ey + matrix(rnorm(i*s,0,1),i,s)
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
# par(bty="L",mfrow=c(2,2))
# boxplot(RMSE~h,data=subset(rmse,Omega==0))
# boxplot(RMSE~h,data=subset(rmse,Omega==5))
# boxplot(RMSE~Omin,data=subset(rmse,Omega==0))
# boxplot(RMSE~Omin,data=subset(rmse,Omega==5))
# boxplot(RMSE~Opct,data=subset(rmse,Omega==0))
# boxplot(RMSE~Opct,data=subset(rmse,Omega==5))
# par(mfrow=c(2,2))
# plot(aggregate(RMSE~h,data=rmse,mean),type="b",lty=2)
# plot(aggregate(RMSE~Omin,data=rmse,mean),type="b",lty=2)
# plot(aggregate(RMSE~Opct,data=rmse,mean),type="b",lty=2)
# boxplot(Variance~h,data=var)
# boxplot(Variance~Omin,data=var)
# boxplot(Variance~Opct,data=var)
# ###
rmse.plot <- reshape(subset(rmse,select=c("Omin","Opct","h","RMSE")),direction="wide",timevar = "Omin",idvar=c("h","Opct"))
rmse.plot <- reshape(rmse.plot,direction="wide",timevar = "Opct",idvar="h")
plot.labels <- substring(colnames(rmse.plot)[-1],6)
plot.labels <- sapply(strsplit(plot.labels,"\\."), function(i) paste(i[1],paste0(i[2:3],collapse="."),sep="; "))
col <- rainbow(ncol(rmse.plot)-1)
pdf("Figures/prior_predictive.pdf")
par(mfrow=c(1,1),bty="L")
matplot(x=rmse.plot[,1],rmse.plot[,-1],type="b",lty=2,col=rep(1:3,3),
        pch=c(15,15,15,16,16,16,17,17,17),
        cex.axis=0.75,cex.lab=0.75,cex=0.75,
        ylab="RMSE",xlab="h")
legend("top",title=expression(N[Omin]),legend=Omin,bty="n",cex=0.65,ncol=3,lty=2,
       col=c(1,3,2),lwd=2)
legend("topright",title=bquote(alpha),legend=Opct,bty="n",cex=0.65,pt.cex=0.85,ncol=3,pch=c(0:2))
dev.off()