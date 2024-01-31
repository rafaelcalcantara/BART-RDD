## Setup
setwd("~/../Git/BART-RDD/")
set.seed(7)
library(parallel)
library(foreach)
library(doParallel)
library(XBART)
## Parallelization
no_cores <- detectCores() - 2
##
s1            <- 5*no_cores
s2            <- 1
n             <- 500
c             <- 0
Omin          <- 1
Opct          <- 0.9
ntrees        <- 5
Nmin          <- 5
num_sweeps    <- 550
burnin        <- 50
p_categorical <- 1
### DGP
x1 <- matrix(2*rbeta(n*s2,2,4)-0.75,n,s2)
x2 <- matrix(rnorm(n*s2,0,1),n,s2)
x3 <- matrix(2*rbeta(n*s2,2,4)-0.5,n,s2)
str <- 1
w <- matrix(rbinom(n*s2,3,0.7),n,s2)+1
z1 <- apply(x1,2,function(i) as.numeric(i>=c))
z2 <- apply(x2,2,function(i) as.numeric(i>=c))
z3 <- apply(x3,2,function(i) as.numeric(i>=c))
mu <- function(x,w) 0.3*w + 0.1*w*x + 0.01*w*(x+0.05)^2 + 0.05*w*x^3
tau <- function(x,w) 0.1 + 0.05*(w-mean(w)) - 0.02*w*x + 0.03*w*(x+0.01)^2 - 0.01*w*x^3
Ey1 <- Ey2 <- Ey3 <- matrix(0,n,s2)
for (i in 1:s2) 
{
  Ey1[,i] <- mu(x1[,i],w[,i]) + tau(x1[,i],w[,i])*z1[,i]
  Ey2[,i] <- mu(x2[,i],w[,i]) + tau(x2[,i],w[,i])*z2[,i]
  Ey3[,i] <- mu(x3[,i],w[,i]) + tau(x3[,i],w[,i])*z3[,i]
}
e1 <- matrix(rnorm(n*s2,0,sd(Ey1)/str),n,s2)
e2 <- matrix(rnorm(n*s2,0,sd(Ey2)/str),n,s2)
e3 <- matrix(rnorm(n*s2,0,sd(Ey3)/str),n,s2)
## Plotting data and functions
pdf("Figures/prior_true_1.pdf")
par(mfrow=c(3,3))
plot(x1[,1],mu(x1[,1],w[,1]),xlab="X",ylab=expression(mu),bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
plot(x1[,1],tau(x1[,1],w[,1]),xlab="X",ylab=expression(tau),bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
plot(x1[,1],Ey1[,1],xlab="X",ylab="E[Y|X]",bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
##
plot(x2[,1],mu(x2[,1],w[,1]),xlab="X",ylab=expression(mu),bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
plot(x2[,1],tau(x2[,1],w[,1]),xlab="X",ylab=expression(tau),bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
plot(x2[,1],Ey2[,1],xlab="X",ylab="E[Y|X]",bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
##
plot(x3[,1],mu(x3[,1],w[,1]),xlab="X",ylab=expression(mu),bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
plot(x3[,1],tau(x3[,1],w[,1]),xlab="X",ylab=expression(tau),bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
plot(x3[,1],Ey3[,1],xlab="X",ylab="E[Y|X]",bty="n",col=w[,1],pch=20)
abline(v=0,lty=2)
dev.off()
###
y1 <- Ey1 + e1
y2 <- Ey2 + e2
y3 <- Ey3 + e3
pdf("Figures/prior_true_2.pdf")
par(mfrow=c(1,3))
plot(x1[,1],y1[,1],bty="n",pch=20,col=z1[,1]+1,xlab="X",ylab="Y")
abline(v=0,lty=2)
plot(x2[,1],y2[,1],bty="n",pch=20,col=z2[,1]+1,xlab="X",ylab="Y")
abline(v=0,lty=2)
plot(x3[,1],y3[,1],bty="n",pch=20,col=z3[,1]+1,xlab="X",ylab="Y")
abline(v=0,lty=2)
dev.off()
## Prior predictive
pars <- runif(6,0.05,0.3)
mu.prior <- function(x) as.vector(cbind(1,x,x^2)%*%pars[1:3])
tau.prior <- function(x) as.vector(cbind(1,x,x^2)%*%pars[4:6])
ys1 <- ys2 <- ys3 <- vector("list",s2)
for (i in 1:s2)
{
  ys1[[i]] <- mu.prior(x1[,i]) + tau.prior(x1[,i])*z1[,i] + matrix(rnorm(n*s1,0,sd(mu.prior(x1[,i]) + tau.prior(x1[,i])*z1[,i])),n,s1)
  ys2[[i]] <- mu.prior(x2[,i]) + tau.prior(x2[,i])*z2[,i] + matrix(rnorm(n*s1,0,sd(mu.prior(x2[,i]) + tau.prior(x2[,i])*z2[,i])),n,s1)
  ys3[[i]] <- mu.prior(x3[,i]) + tau.prior(x3[,i])*z3[,i] + matrix(rnorm(n*s1,0,sd(mu.prior(x3[,i]) + tau.prior(x3[,i])*z3[,i])),n,s1)
}
pdf("Figures/prior_synthetic.pdf")
par(mfrow=c(1,3))
plot(x1[,1],ys1[[1]][,1],bty="n",pch=20,col=z1[,1]+1,xlab="X",ylab=expression(Y[prior]))
abline(v=0,lty=2)
plot(x2[,1],ys2[[1]][,1],bty="n",pch=20,col=z2[,1]+1,xlab="X",ylab=expression(Y[prior]))
abline(v=0,lty=2)
plot(x3[,1],ys3[[1]][,1],bty="n",pch=20,col=z3[,1]+1,xlab="X",ylab=expression(Y[prior]))
abline(v=0,lty=2)
dev.off()
## Functions
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ifelse(is.matrix(y),ys <- y[,i],ys <- y)
  ifelse(is.matrix(w0),ws <- w0[,i],ws <- w0)
  ifelse(is.matrix(x),xs <- x[,i],xs <- x)
  ifelse(length(h)>1,hs <- h[i],hs <- h)
  ifelse(length(Omin)>1,os <- Omin[i],os <- Omin)
  ifelse(length(Opct)>1,ops <- Opct[i],ops <- Opct)
  fit <- XBART::XBCF.rd(ys, ws, xs, c,
                        Owidth = hs, Omin = os, Opct = ops,
                        num_trees_mod = ntrees,
                        num_trees_con = ntrees,
                        num_cutpoints = n,
                        num_sweeps = num_sweeps,
                        burnin = burnin, Nmin = Nmin,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(ys)/ntrees,
                        tau_mod = 0.5*var(ys)/ntrees)
  test <- -hs+c<=xs & xs<=hs+c
  pred <- XBART::predict.XBCFrd(fit,ws[test],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}
####
h1 <- round(rbind(0.1*apply(x1,2,sd),0.25*apply(x1,2,sd),0.5*apply(x1,2,sd),apply(x1,2,sd)),2)
h2 <- round(rbind(0.1*apply(x2,2,sd),0.25*apply(x2,2,sd),0.5*apply(x2,2,sd),apply(x2,2,sd)),2)
h3 <- round(rbind(0.1*apply(x3,2,sd),0.25*apply(x3,2,sd),0.5*apply(x3,2,sd),apply(x3,2,sd)),2)
####
rmse1 <- rmse2 <- rmse3 <- matrix(0,ncol(h1),nrow(h1))
for (i in 1:s2)
{
  y <- ys1[[i]]
  x <- x1[,i]
  w0 <- w[,i]
  ate <- tau.prior(0)
  for (j in 1:nrow(h1))
  {
    h <- h1[j,i]
    print(paste0("DGP: ",i))
    cl <- makeCluster(no_cores,type="SOCK")
    registerDoParallel(cl)
    clusterExport(cl,varlist=ls())
    time <- system.time({
      out <- parLapply(cl,1:s1,fit)
    })
    stopCluster(cl)
    rmse1[i,j] <- sqrt(mean((colMeans(sapply(out,colMeans))-ate)^2))
    print(time)
    gc()
  }
}
for (i in 1:s2)
{
  y <- ys2[[i]]
  x <- x2[,i]
  w0 <- w[,i]
  ate <- tau.prior(0)
  for (j in 1:nrow(h2))
  {
    h <- h2[j,i]
    print(paste0("DGP: ",i))
    cl <- makeCluster(no_cores,type="SOCK")
    registerDoParallel(cl)
    clusterExport(cl,varlist=ls())
    time <- system.time({
      out <- parLapply(cl,1:s1,fit)
    })
    stopCluster(cl)
    rmse2[i,j] <- sqrt(mean((colMeans(sapply(out,colMeans))-ate)^2))
    print(time)
    gc()
  }
}
for (i in 1:s2)
{
  y <- ys3[[i]]
  x <- x3[,i]
  w0 <- w[,i]
  ate <- tau.prior(0)
  for (j in 1:nrow(h3))
  {
    h <- h3[j,i]
    print(paste0("DGP: ",i))
    cl <- makeCluster(no_cores,type="SOCK")
    registerDoParallel(cl)
    clusterExport(cl,varlist=ls())
    time <- system.time({
      out <- parLapply(cl,1:s1,fit)
    })
    stopCluster(cl)
    rmse3[i,j] <- sqrt(mean((colMeans(sapply(out,colMeans))-ate)^2))
    print(time)
    gc()
  }
}
###
rmse <- data.frame(DGP1=rmse1[1,],DGP2=rmse2[1,],DGP3=rmse3[1,],row.names = paste("h",1:4,sep="_"))
rmse <- round(rmse,3)
h.tab <- data.frame(DGP1=h1,DGP2=h2,DGP3=h3,row.names = paste("h",1:4,sep="_"))
rmse.tab <- rbind(rep(c("h","RMSE"),3),
                  c(h.tab[1,1],rmse[1,1],h.tab[1,2],rmse[1,2],h.tab[1,3],rmse[1,3]),
                  c(h.tab[2,1],rmse[2,1],h.tab[2,2],rmse[2,2],h.tab[2,3],rmse[2,3]),
                  c(h.tab[3,1],rmse[3,1],h.tab[3,2],rmse[3,2],h.tab[3,3],rmse[3,3]),
                  c(h.tab[4,1],rmse[4,1],h.tab[4,2],rmse[4,2],h.tab[4,3],rmse[4,3]))
colnames(rmse.tab) <- c("DGP1","","DGP2","","DGP3","")
###
y <- y1[,1]
x <- x1[,1]
w0 <- w[,1]
h <- h1[which(rmse1[1,]==min(rmse1[1,])),1]
Omin <- 1:no_cores
Opct.grid <- c(0.7,0.9)
est.omin <- vector("list",2)
het.omin <- matrix(0,no_cores,2)
for (i in 1:length(Opct.grid))
{
  Opct <- Opct.grid[i]
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:length(Omin),fit)
  })
  stopCluster(cl)
  est.omin[[i]] <- t(apply(sapply(out,colMeans),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  het.omin[,i] <- apply(sapply(out,rowMeans),2,function(i) length(table(i)))
}
##
Opct <- seq(0.7,0.9,length=no_cores)
Omin.grid <- c(1,10)
est.opct <- vector("list",2)
het.opct <- matrix(0,no_cores,2)
for (i in 1:length(Omin.grid))
{
  Omin <- Omin.grid[i]
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:length(Opct),fit)
  })
  stopCluster(cl)
  est.opct[[i]] <- t(apply(sapply(out,colMeans),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  het.opct[,i] <- apply(sapply(out,rowMeans),2,function(i) length(table(i)))
}
##
pdf("Figures/prior_omin1.pdf")
par(mfrow=c(2,2))
matplot(est.omin[[1]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.7")))
matplot(est.omin[[2]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.9")))
matplot(Opct,est.opct[[1]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 1")))
matplot(Opct,est.opct[[2]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 10")))
dev.off()
pdf("Figures/prior_het_omin1.pdf")
par(mfrow=c(2,2))
matplot(het.omin[,1],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.7")))
matplot(het.omin[,2],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.9")))
matplot(Opct,het.opct[,1],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 1")))
matplot(Opct,het.opct[,2],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 10")))
dev.off()
###
y <- y2[,1]
x <- x2[,1]
w0 <- w[,1]
h <- h2[which(rmse2[1,]==min(rmse2[1,])),1]
Omin <- 1:no_cores
Opct.grid <- c(0.7,0.9)
est.omin <- vector("list",2)
het.omin <- matrix(0,no_cores,2)
for (i in 1:length(Opct.grid))
{
  Opct <- Opct.grid[i]
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:length(Omin),fit)
  })
  stopCluster(cl)
  est.omin[[i]] <- t(apply(sapply(out,colMeans),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  het.omin[,i] <- apply(sapply(out,rowMeans),2,function(i) length(table(i)))
}
##
Opct <- seq(0.7,0.9,length=no_cores)
Omin.grid <- c(1,10)
est.opct <- vector("list",2)
het.opct <- matrix(0,no_cores,2)
for (i in 1:length(Omin.grid))
{
  Omin <- Omin.grid[i]
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:length(Opct),fit)
  })
  stopCluster(cl)
  est.opct[[i]] <- t(apply(sapply(out,colMeans),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  het.opct[,i] <- apply(sapply(out,rowMeans),2,function(i) length(table(i)))
}
##
pdf("Figures/prior_omin2.pdf")
par(mfrow=c(2,2))
matplot(est.omin[[1]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.7")))
matplot(est.omin[[2]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.9")))
matplot(Opct,est.opct[[1]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 1")))
matplot(Opct,est.opct[[2]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 10")))
dev.off()
pdf("Figures/prior_het_omin2.pdf")
par(mfrow=c(2,2))
matplot(het.omin[,1],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.7")))
matplot(het.omin[,2],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.9")))
matplot(Opct,het.opct[,1],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 1")))
matplot(Opct,het.opct[,2],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 10")))
dev.off()
###
y <- y3[,1]
x <- x3[,1]
w0 <- w[,1]
h <- h3[which(rmse3[1,]==min(rmse3[1,])),1]
Omin <- 1:no_cores
Opct.grid <- c(0.7,0.9)
est.omin <- vector("list",2)
het.omin <- matrix(0,no_cores,2)
for (i in 1:length(Opct.grid))
{
  Opct <- Opct.grid[i]
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:length(Omin),fit)
  })
  stopCluster(cl)
  est.omin[[i]] <- t(apply(sapply(out,colMeans),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  het.omin[,i] <- apply(sapply(out,rowMeans),2,function(i) length(table(i)))
}
##
Opct <- seq(0.7,0.9,length=no_cores)
Omin.grid <- c(1,10)
est.opct <- vector("list",2)
het.opct <- matrix(0,no_cores,2)
for (i in 1:length(Omin.grid))
{
  Omin <- Omin.grid[i]
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:length(Opct),fit)
  })
  stopCluster(cl)
  est.opct[[i]] <- t(apply(sapply(out,colMeans),2,function(i) c(mean(i),quantile(i,c(0.025,0.975)))))
  het.opct[,i] <- apply(sapply(out,rowMeans),2,function(i) length(table(i)))
}
##
pdf("Figures/prior_omin3.pdf")
par(mfrow=c(2,2))
matplot(est.omin[[1]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.7")))
matplot(est.omin[[2]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.9")))
matplot(Opct,est.opct[[1]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 1")))
matplot(Opct,est.opct[[2]],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 10")))
dev.off()
pdf("Figures/prior_het_omin3.pdf")
par(mfrow=c(2,2))
matplot(het.omin[,1],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.7")))
matplot(het.omin[,2],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(N[Omin]),main=expression(paste(alpha," = 0.9")))
matplot(Opct,het.opct[,1],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 1")))
matplot(Opct,het.opct[,2],col=1,type=c("b","l","l"),pch=20,lty=c(1,2,2),bty="n",ylab="",xlab=expression(alpha),main=expression(paste(N[Omin]," = 10")))
dev.off()
###
stargazer::stargazer(rmse.tab,summary = F,rownames = F,
                     title="Candidate $h$ values and RMSE for each (X,W) sample",
                     label="tab:rmse.prior",
                     font.size = "small")