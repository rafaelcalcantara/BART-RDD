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
  y <- ys[,i]
  fit <- XBART::XBCF.rd(y, w, x, c,
                        Owidth = hs, Omin = Om, Opct = Op,
                        num_trees_mod = ntc,
                        num_trees_con = ntc,
                        num_cutpoints = n,
                        num_sweeps = num_sweeps,
                        burnin = burnin, Nmin = nm,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = tc*var(y)/ntc,
                        tau_mod = tm*var(y)/ntm)
  test <- -hs+c<=x & x<=hs+c
  pred <- XBART::predict.XBCFrd(fit,w[test],rep(c,sum(test)))
  pred <- rowMeans(pred$tau.adj[,(burnin+1):num_sweeps])
  sqrt(mean((pred-cate[test])^2))
}
### Prior predictive
#### Setup
mu.prior <- function(x,w) w + 1/(1+exp(-5*x)) + (1-abs(x-c))*sin(x)/10
tau.prior <- function(x,w,tau.bar) tau.bar - log(x+1)/50 + 0.01*(w - mean(w))
ate <- 0.4
c <- 0
s <- 11 ## no of samples of th synthetic DGP
N <- c(500)
Omin <- c(1,5,10)
Opct <- seq(0.6,0.9,length=2)
tau.sig <- c(0.5,1,2)
ntrees <- c(5,10,15)
Nmin <- c(5,10,15)
#### Loop
params <- c("N","ATE","Omin","Opct","h","tcon","tmod","Ntcon","Ntmod","Nmin")
num_sweeps    <- 150
burnin        <- 50
p_categorical <- 0
for (n in N)
{
  ## Data for calibrating the prior to the simulation DGPs
  x <- 2*rbeta(n,2,4)-0.75
  z <- as.numeric(x>=c)
  w <- runif(n)
  cate <- tau.prior(c,w,ate)
  h <- quantile(abs(x),c(0.05,0.1,0.15,0.2))
  h <- 0.068
  Ey <- mu.prior(x,w) + tau.prior(x,w,ate)*z
  ys <-  Ey + matrix(rnorm(n*s,0,0.1*sd(tau.prior(x,w,ate))),n,s)
  for (tc in tau.sig)
  {
    for (tm in tau.sig)
    {
      for (ntc in ntrees)
      {
        for (ntm in ntrees)
        {
          for (nm in Nmin)
          {
            for (Om in Omin)
            {
              for (Op in Opct)
              {
                for (hs in h)
                {
                  print(paste(params,c(n,ate,Om,Op,hs,tc,tm,ntc,ntm,nm),sep=": "))
                  ## Fit the model
                  cl <- makeCluster(no_cores,type="SOCK")
                  registerDoParallel(cl)
                  clusterExport(cl,varlist=ls())
                  time <- system.time({
                    out <- parSapply(cl,1:s,fit)
                  })
                  stopCluster(cl)
                  print(time)
                  assign(paste(params,c(n,ate,Om,Op,hs,tc,tm,ntc,ntm,nm),sep=".",collapse="_"),out)
                }
              }
            }
          }
        }
      }
    }
  }
}
save.image("prior_predictive.RData")
# load("prior_predictive.RData")
## Results
### RMSE
objects <- mget(ls())
sel <- grep(paste0("ATE.",ate),names(objects))
p <- do.call("rbind",strsplit(names(objects[sel]),"_"))
p <- apply(p,2, function(i) sub("N.|ATE.|Omin.|Opct.|h.|tcon.|tmod.|Ntcon.|Ntmod.|Nmin.","",i))
rmse <- sapply(objects[sel],mean)
rmse <- cbind(p,RMSE=rmse)
colnames(rmse) <- c(params,"RMSE")
rmse <- apply(rmse,2,as.numeric)
rmse <- data.frame(rmse)
### Plots
rmse.plot <- reshape(subset(rmse,select=c("Omin","Opct","h","RMSE")),direction="wide",timevar = "Omin",idvar=c("h","Opct"))
rmse.plot <- reshape(rmse.plot,direction="wide",timevar = "Opct",idvar="h")
plot.labels <- substring(colnames(rmse.plot)[-1],6)
plot.labels <- sapply(strsplit(plot.labels,"\\."), function(i) paste(i[1],paste0(i[2:3],collapse="."),sep="; "))
col <- rainbow(ncol(rmse.plot)-1)
# pdf("Figures/prior_predictive.pdf")
par(mfrow=c(1,1),bty="L")
matplot(x=rmse.plot[,1],rmse.plot[,-1],type="b",lty=2,col=rep(1:3,3),
        pch=c(15,15,15,16,16,16,17,17,17),
        cex.axis=0.75,cex.lab=0.75,cex=0.75,
        ylab="RMSE",xlab="h")
legend("top",title=expression(N[Omin]),legend=Omin,bty="n",cex=0.65,ncol=3,lty=2,
       col=1:3,lwd=2)
legend("topright",title=bquote(alpha),legend=Opct,bty="n",cex=0.65,pt.cex=0.85,ncol=3,pch=c(0:2))
# dev.off()