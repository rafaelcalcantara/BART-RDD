## Setup
set.seed(7)
library(parallel)
library(foreach)
library(doParallel)
# devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
library(rdrobust)
library(bayesrdd)
library(rpart)
library(rpart.plot)
library(xtable)
## Parallelization
no_cores <- detectCores() - 2
##
s1            <- no_cores
Omin          <- 1
Opct          <- 0.9
ntrees        <- 5
Nmin          <- 5
num_sweeps    <- 120
burnin        <- 20
p_categorical <- 5
### Read data
setwd("~/../Git/BART-RDD/")
data <- read.csv("Application/gpa.csv")
y <- data$nextGPA
x <- data$X
w <- data[,4:11]
c <- 0
n <- nrow(data)
z <- as.numeric(x>c)
### Summary statistics
#### Tables
sum.stat <- function(dat) t(apply(dat,2,function(i) c(Mean=mean(i),SD=sd(i),Min=min(i),Max=max(i),Cor=cor(dat[,1],i))))
sum0 <- sum.stat(cbind(y,x,w))
sum1 <- by(cbind(y,x,w[,1:5]),data$loc_campus1,sum.stat)[[2]]
sum2 <- by(cbind(y,x,w[,1:5]),data$loc_campus2,sum.stat)[[2]]
sum3 <- by(cbind(y,x,w[,1:5]),data$loc_campus3,sum.stat)[[2]]
sum.tab <- lapply(1:7, function(i) round(rbind(sum0[i,],sum1[i,],sum2[i,],sum3[i,]),3))
sum.tab <- lapply(1:7, function(i) cbind(Sample=c("Full","Campus 1","Campus 2","Campus 3"),sum.tab[[i]]))
sum.tab <- do.call("rbind",sum.tab)
sum.tab <- cbind(unlist(lapply(c("Y","X",names(w)[1:5]), function(i) c(i,rep("",3)))),sum.tab)
obs <- c(n,table(data$loc_campus1)[2],table(data$loc_campus2)[2],table(data$loc_campus3)[2])
stargazer::stargazer(sum.tab,summary=F,rownames = F,
                     label = "tab:sum.stat", title = "Summary statistics",
                     table.placement = "!htpb", font.size = "small",
                     notes = c("\\tiny Sample size:",paste(c("\\tiny Total:","Campus 1:","Campus 2:","Campus 3:"),obs,collapse="; ")))
#### Plots
##### W
l0 <- lapply(1:ncol(w), function(i) predict(loess(w~x,data=data.frame(w=w[x<c & x>-0.5,i],x=x[x<c & x>-0.5])),se=T)[1:2])
l0 <- lapply(l0,function(i) cbind(i$fit,i$fit-1.96*i$se.fit,i$fit+1.96*i$se.fit))
l1 <- lapply(1:ncol(w), function(i) predict(loess(w~x,data=data.frame(w=w[x>c & x<0.5,i],x=x[x>c & x<0.5])),se=T)[1:2])
l1 <- lapply(l1,function(i) cbind(i$fit,i$fit-1.96*i$se.fit,i$fit+1.96*i$se.fit))
##### Y
yhat0 <- predict(loess(y~x,data=data.frame(y=y[x<c & x>-0.5],x=x[x<c & x>-0.5])),se=T)[1:2]
yhat0 <- cbind(yhat0$fit,yhat0$fit-1.96*yhat0$se.fit,yhat0$fit+1.96*yhat0$se.fit)
yhat1 <- predict(loess(y~x,data=data.frame(y=y[x>c & x<0.5],x=x[x>c & x<0.5])),se=T)[1:2]
yhat1 <- cbind(yhat1$fit,yhat1$fit-1.96*yhat1$se.fit,yhat1$fit+1.96*yhat1$se.fit)
#####
xplot <- x[x>-0.5 & x<0.5]
pdf("Figures/application_summary.pdf")
par(mfrow=c(3,3))
matplot(xplot,rbind(yhat0,yhat1),xlab="X",ylab="Y",bty="n",type="n")
matlines(sort(xplot[xplot<c]),yhat0[order(xplot[xplot<c]),],col=1,lty=c(1,2,2))
matlines(sort(xplot[xplot>c]),yhat1[order(xplot[xplot>c]),],col=1,lty=c(1,2,2))
abline(v=c,lty=2)
for (i in 1:ncol(w))
{
  matplot(xplot,rbind(l0[[i]],l1[[i]]),xlab="X",ylab=colnames(w)[i],bty="n",type="n")
  matlines(sort(xplot[xplot<c]),l0[[i]][order(xplot[xplot<c]),],col=1,lty=c(1,2,2))
  matlines(sort(xplot[xplot>c]),l1[[i]][order(xplot[xplot>c]),],col=1,lty=c(1,2,2))
  abline(v=c,lty=2)
}
dev.off()
### Prior elicitation
pars <- runif(6,0.05,0.3)
mu.prior <- function(x) as.vector(cbind(1,x,x^2)%*%pars[1:3])
tau.prior <- function(x) as.vector(cbind(1,x,x^2)%*%pars[4:6])
ys <- matrix(0,nrow(data),s1)
for (i in 1:s1)
{
  ys[,i] <- mu.prior(x) + tau.prior(x)*z + rnorm(n,0,sd(mu.prior(x) + tau.prior(x)*z)*0.5)
  
}
h.grid <- sapply(c(0.1,0.25,0.5,0.75,1),function(i) i*sd(x))
fit <- function(i)
{
  print(paste0("Sample: ",i))
  ys1 <- ys[,i]
  fit <- XBART::XBCF.rd(ys1, w, x, c,
                        Owidth = hs, Omin = Omin, Opct = Opct,
                        num_trees_mod = ntrees,
                        num_trees_con = ntrees,
                        num_cutpoints = n,
                        num_sweeps = num_sweeps,
                        burnin = burnin, Nmin = Nmin,
                        p_categorical_con = p_categorical,
                        p_categorical_mod = p_categorical,
                        tau_con = 2*var(ys1)/ntrees,
                        tau_mod = 0.5*var(ys1)/ntrees)
  test <- -hs+c<=x & x<=hs+c
  pred <- XBART::predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
  pred$tau.adj[,(burnin+1):num_sweeps]
}
rmse <- rep(0,length(h.grid))
ate <- tau.prior(0)
for (i in 1:length(h.grid))
{
  print(paste0("h: ",h.grid[i]))
  hs <- h.grid[i]
  cl <- makeCluster(no_cores,type="SOCK")
  registerDoParallel(cl)
  clusterExport(cl,varlist=ls())
  time <- system.time({
    out <- parLapply(cl,1:s1,fit)
  })
  stopCluster(cl)
  rmse[i] <- sqrt(mean((colMeans(sapply(out,colMeans))-ate)^2))
  print(time)
  gc()
}
h <- h.grid[which(h.grid==min(h.grid))]
rmse.tab <- cbind(h=h.grid,RMSE=rmse)
rmse.tab <- round(rmse.tab,3)
table(-h<x & x<h)
stargazer::stargazer(rmse.tab,summary = F,rownames = F,font.size = "small",
                     label="tab:rmse.data",title="Results from prior elicitation")
### Summary for identification strip
# h <-h.grid[1]
sum.stat.h <- function(dat,id0,id1) 
{
  out <- round(cbind(t(apply(dat[id0,],2,function(i) c(Mean=mean(i),SD=sd(i)))),
                     t(apply(dat[id1,],2,function(i) c(Mean=mean(i),SD=sd(i))))), 3)
  return(out)
}
id0 <- -h < x & x < 0
id1 <- 0 < x & x< h
sum0 <- sum.stat.h(cbind(y,x,w),id0,id1)
sum1 <- sum.stat.h(cbind(y,x,w[,1:5])[data$loc_campus1==1,],id0[data$loc_campus1==1],id1[data$loc_campus1==1])
sum2 <- sum.stat.h(cbind(y,x,w[,1:5])[data$loc_campus2==1,],id0[data$loc_campus2==1],id1[data$loc_campus2==1])
sum3 <- sum.stat.h(cbind(y,x,w[,1:5])[data$loc_campus3==1,],id0[data$loc_campus3==1],id1[data$loc_campus3==1])
sum.tab <- lapply(1:7, function(i) rbind(sum0[i,],sum1[i,],sum2[i,],sum3[i,]))
sum.tab <- lapply(1:7, function(i) cbind(Sample=c("Full","Campus 1","Campus 2","Campus 3"),sum.tab[[i]]))
sum.tab <- do.call("rbind",sum.tab)
sum.tab <- cbind(unlist(lapply(c("Y","X",names(w)[1:5]), function(i) c(i,rep("",3)))),sum.tab)
sum.tab <- rbind(colnames(sum.tab),sum.tab)
colnames(sum.tab) <- c("","","(1)","","(2)","")
obs0 <- c(sum(id0),table(data$loc_campus1[id0])[2],table(data$loc_campus2[id0])[2],table(data$loc_campus3[id0])[2])
obs1 <- c(sum(id1),table(data$loc_campus1[id1])[2],table(data$loc_campus2[id1])[2],table(data$loc_campus3[id1])[2])
stargazer::stargazer(sum.tab,summary=F,rownames = F,
                     label = "tab:sum.stat.h", title = "Summary statistics - identification strip",
                     table.placement = "!htpb", font.size = "small",
                     notes = c("\\tiny Sample size (control/treatment):",
                               paste(c("\\tiny Total:","Campus 1:","Campus 2:","Campus 3:"),
                                     paste(obs0,obs1,sep="/"),collapse="; ")))
### Fitting the model
burnin <- 50
num_sweeps <- 550
fit <- XBART::XBCF.rd(y, w, x, c,
                      Owidth = h, Omin = Omin, Opct = Opct,
                      num_trees_mod = ntrees,
                      num_trees_con = ntrees,
                      num_cutpoints = n,
                      num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical,
                      tau_con = 2*var(y)/ntrees,
                      tau_mod = 0.5*var(y)/ntrees,
                      parallel = T, nthread = no_cores)
test <- -h+c<=x & x<=h+c
pred <- XBART::predict.XBCFrd(fit,w[test,],rep(c,sum(test)))
pred <- pred$tau.adj[,(burnin+1):num_sweeps]
saveRDS(pred,"Results/bart_rdd_posterior.rds")
pred <- readRDS("Results/bart_rdd_posterior.rds")
### Other estimators
#### Local polynomial
llr <- rdrobust::rdrobust(y,x,c,covs=w[,1:7],masspoints = F)
#### BCF
fit <- XBART::XBCF.discrete(y=y, Z=z, X_con = as.matrix(cbind(x,w)), X_mod = as.matrix(cbind(x,w)),
                      num_trees_mod = ntrees,
                      num_trees_con = ntrees,
                      num_cutpoints = n,
                      num_sweeps = num_sweeps,
                      burnin = burnin, Nmin = Nmin,
                      p_categorical_con = p_categorical,
                      p_categorical_mod = p_categorical,
                      tau_con = 2*var(y)/ntrees,
                      tau_mod = 0.5*var(y)/ntrees,
                      parallel = T, nthread = no_cores)
pred.bcf <- XBART::predict.XBCFdiscrete(fit,X_con = as.matrix(cbind(0,w)), X_mod = as.matrix(cbind(0,w)),Z=z,pihat=z,burnin=burnin)
pred.bcf <- pred.bcf$tau.adj[,(burnin+1):num_sweeps]
saveRDS(pred.bcf,"Results/bcf_posterior.rds")
pred.bcf <- readRDS("Results/bcf_posterior.rds")
#### SBART
fit <- XBART::XBART(y, as.matrix(cbind(x,w,z)), num_trees = ntrees,
                    num_cutpoints = n, num_sweeps = num_sweeps,
                    burnin = burnin, Nmin = Nmin,
                    p_categorical = p_categorical,
                    tau = var(y)/ntrees, parallel=T, nthread = no_cores)
pred1 <- XBART::predict.XBART(fit,cbind(rep(0,n),as.matrix(w),rep(1,n)))[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit,cbind(rep(0,n),as.matrix(w),rep(0,n)))[,(burnin+1):num_sweeps]
pred.sbart <- pred1-pred0
saveRDS(pred.sbart,"Results/sbart_posterior.rds")
pred.sbart <- readRDS("Results/sbart_posterior.rds")
#### TBART
fit1 <- XBART::XBART(y[z==1], as.matrix(cbind(x,w))[z==1,], num_trees = ntrees,
                     num_cutpoints = sum(z==1), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==1])/ntrees, parallel=T, nthread = no_cores)
fit0 <- XBART::XBART(y[z==0], as.matrix(cbind(x,w))[z==0,], num_trees = ntrees,
                     num_cutpoints = sum(z==0), num_sweeps = num_sweeps,
                     burnin = burnin, Nmin = Nmin,
                     p_categorical = p_categorical,
                     tau = var(y[z==0])/ntrees, parallel=T, nthread = no_cores)
pred1 <- XBART::predict.XBART(fit1,cbind(rep(0,n),as.matrix(w)))[,(burnin+1):num_sweeps]
pred0 <- XBART::predict.XBART(fit0,cbind(rep(0,n),as.matrix(w)))[,(burnin+1):num_sweeps]
pred.tbart <- pred1-pred0
saveRDS(pred.tbart,"Results/tbart_posterior.rds")
pred.tbart <- readRDS("Results/tbart_posterior.rds")
#### CGS
cgs.sample <- -1 < x & x < 1
p <- c(.9,.1)
mz <- c(4,3)
mztau <- c(3,2)
mw <- rep(5,3)
lamstmean0_ <- rep(1,6)
lamstsd0_ <- 5*rep(1,6)
s2mean0_ <- .3
s2sd0_ <- 1
pred.cgs <- bayesrddest(y = y[cgs.sample],
                        z = x[cgs.sample],
                        V = as.matrix(w)[,4:8],
                        W = as.matrix(w)[,1:3],
                        tau = c,
                        p = p,
                        mw = mw,
                        mz = mz,
                        mztau = mztau,
                        beta0_ = rep(0,15),
                        lamstmean0_ = lamstmean0_,
                        lamstsd0_ = lamstsd0_,
                        s2mean0_ = s2mean0_,
                        s2sd0_ = s2sd0_,
                        distribution = "gaussian",
                        hetero = FALSE,
                        n0=100, m=1000)
saveRDS(pred.cgs,"Results/pred_cgs.rds")
pred.cgs <- readRDS("Results/pred_cgs.rds")
### ATE
ate.post <- colMeans(pred)
ate.sum <- cbind(Mean=mean(ate.post), SD=sd(ate.post),
                 `2.5%`=quantile(ate.post,0.025),`97.5%`=quantile(ate.post,0.975),
                 Median=median(ate.post),Min=min(ate.post),Max=max(ate.post))
ate.sum <- round(ate.sum,3)
stargazer::stargazer(ate.sum,summary=F,rownames = F,
                     label="tab:ate.results",
                     title="BART-RDD posterior summary for the ATE",
                     font.size = "small")
ate.bcf <- colMeans(pred.bcf)
ate.sbart <- colMeans(pred.sbart)
ate.tbart <- colMeans(pred.tbart)
ate.other <- cbind(BCF=c(round(mean(ate.bcf),3),paste0("(",round(quantile(ate.bcf,0.025),3),",",round(quantile(ate.bcf,0.975),3),")")),
                   `S-BART`=c(round(mean(ate.sbart),3),paste0("(",round(quantile(ate.sbart,0.025),3),",",round(quantile(ate.sbart,0.975),3),")")),
                   `T-BART`=c(round(mean(ate.tbart),3),paste0("(",round(quantile(ate.tbart,0.025),3),",",round(quantile(ate.tbart,0.975),3),")")),
                   LLR=c(round(llr$Estimate[2],3),paste0("(",round(llr$ci[2,1],3),",",round(llr$ci[2,2],3),")")),
                   CGS=c(round(mean(pred.cgs$atem),3),paste0("(",round(quantile(pred.cgs$atem,0.025),3),",",round(quantile(pred.cgs$atem,0.975),3),")")))
stargazer::stargazer(ate.other, summary = F, rownames = F,
                     label="tab:ate.others",
                     title="ATE point estimate and 95% confidence interval for different estimators",
                     font.size = "small")
###
cate <- rpart(y~.,data.frame(y=rowMeans(pred),w[test,]),control = rpart.control(cp=0.03))
cate.campus1 <- rpart(y~.,data.frame(y=rowMeans(pred)[data$loc_campus1[test]==1],w[test & data$loc_campus1==1,]),control = rpart.control(cp=0.1))
cate.campus2 <- rpart(y~.,data.frame(y=rowMeans(pred)[data$loc_campus2[test]==1],w[test & data$loc_campus2==1,]),control = rpart.control(cp=0.047))
cate.campus3 <- rpart(y~.,data.frame(y=rowMeans(pred)[data$loc_campus3[test]==1],w[test & data$loc_campus3==1,]),control = rpart.control(cp=0.04))
pdf("Figures/cate_gpa.pdf")
par(mfrow=c(2,2))
rpart.plot(cate,main="Full sample")
rpart.plot(cate.campus1,main="Campus 1")
rpart.plot(cate.campus2,main="Campus 2")
rpart.plot(cate.campus3,main="Campus 3")
dev.off()
####
cate.hspct.leq.37.campus1 <- do.call("cbind",by(pred,data$hsgrade_pct[test] < 37
                                           & data$loc_campus1[test] == 1, colMeans))
cate.hspct.leq.36.campus2 <- do.call("cbind",by(pred,data$hsgrade_pct[test] < 36
                                                & data$loc_campus2[test] == 1, colMeans))
cate.age.leq.19.campus3 <- do.call("cbind",by(pred,data$age_at_entry[test] < 19
                                                & data$loc_campus3[test] == 1, colMeans))
cate.hspct.geq.37.campus1 <- do.call("cbind",by(pred,data$hsgrade_pct[test] >= 37
                                                & data$loc_campus1[test] == 1, colMeans))
cate.hspct.geq.36.campus2 <- do.call("cbind",by(pred,data$hsgrade_pct[test] >= 36
                                                & data$loc_campus2[test] == 1, colMeans))
cate.age.geq.19.campus3 <- do.call("cbind",by(pred,data$age_at_entry[test] >= 19
                                                & data$loc_campus3[test] == 1, colMeans))
cate.campus <- do.call("cbind",by(pred,as.factor(data$loc_campus1+2*data$loc_campus2+3*data$loc_campus3)[test],colMeans))
####
pct1 <- mean((cate.hspct.leq.37.campus1[,2]-cate.hspct.geq.37.campus1[,2])>0)
pct2 <- mean((cate.hspct.leq.36.campus2[,2]-cate.hspct.geq.36.campus2[,2])>0)
pct3 <- mean((cate.age.leq.19.campus3[,2]-cate.age.geq.19.campus3[,2])>0)
pct4 <- mean((cate.campus[,2]-cate.campus[,1])>0)
pct5 <- mean((cate.campus[,3]-cate.campus[,1])>0)
pct6 <- mean((cate.campus[,3]-cate.campus[,2])>0)
pct1;pct2;pct3;pct4;pct5;pct6
####
d1 <- density(cate.hspct.leq.37.campus1[,2]-cate.hspct.geq.37.campus1[,2])
d2 <- density(cate.hspct.leq.36.campus2[,2]-cate.hspct.geq.36.campus2[,2])
d3 <- density(cate.age.leq.19.campus3[,2]-cate.age.geq.19.campus3[,2])
d4 <- density(cate.campus[,2]-cate.campus[,1])
d5 <- density(cate.campus[,3]-cate.campus[,1])
d6 <- density(cate.campus[,3]-cate.campus[,2])
pdf("Figures/cate_difference.pdf")
par(mfrow=c(2,2))
plot(d1,bty="n",main="Campus 1 - hsgrade_pct > 37",ylab="",xlab="Difference in subgroup average treatment effect")
polygon(c(0,d1$x[d1$x>0]),c(0,d1$y[d1$x>0]),col=rgb(96/256,96/256,96/256,0.4),border=1)
plot(d2,bty="n",main="Campus 2 - hsgrade_pct > 36",ylab="",xlab="Difference in subgroup average treatment effect")
polygon(c(0,d2$x[d2$x>0]),c(0,d2$y[d2$x>0]),col=rgb(96/256,96/256,96/256,0.4),border=1)
plot(d3,bty="n",main="Campus 3 - age_at_entry > 19",ylab="",xlab="Difference in subgroup average treatment effect")
polygon(c(0,d3$x[d3$x>0]),c(0,d3$y[d3$x>0]),col=rgb(96/256,96/256,96/256,0.4),border=1)
plot(d4,type="n",bty="n",
     xlim=c(min(d4$x,d5$x,d6$x),max(d4$x,d5$x,d6$x)),
     ylim=c(min(d4$y,d5$y,d6$y),max(d4$y,d5$y,d6$y)),
     main="",
     ylab="",xlab="Difference in campus average treatment effect")
lines(d4,col=rgb(1,0,0))
lines(d5,col=rgb(0,1,0))
lines(d6,col=rgb(0,0,1))
polygon(c(0,d4$x[d4$x>0]),c(0,d4$y[d4$x>0]),col=rgb(1,0,0,0.1),border=rgb(1,0,0))
polygon(c(0,d5$x[d5$x>0]),c(0,d5$y[d5$x>0]),col=rgb(0,1,0,0.1),border=rgb(0,1,0))
polygon(c(0,d6$x[d6$x>0]),c(0,d6$y[d6$x>0]),col=rgb(0,0,1,0.1),border=rgb(0,0,1))
legend("topleft",legend=c("2 - 1","3 - 1","3 - 2"),
       col=c(rgb(1,0,0),rgb(0,1,0),rgb(0,0,1)),lty=1,cex=0.75)
dev.off()
####
cate.hspct.geq.43.campus3 <- do.call("cbind",by(pred,data$hsgrade_pct[test] >= 43
                                           & data$loc_campus3[test] != 1, colMeans))
####
plot(density(cate.hspct.leq.43.campus1[,2]-cate.hspct.geq.43.campus3[,2]))
###
## sample <- -0.3<=x & x<=0.3 ## For XBCF and loess plot
## ## Plot data
## png("Figures/gpa_data.png")
## plot(x,y,pch=21,bg="azure",cex=0.5,xlab="1st Year GPA",ylab="2nd Year GPA",bty="n")
## dev.off()
## ###
## l0 <- loess(y~x,data=data.frame(y=y[x<c & sample],x=x[x<c & sample]))
## l0 <- predict(l0,se=T)
## l0 <- data.frame(Est=l0$fit, LI=l0$fit-1.96*l0$se.fit,
##                  UI=l0$fit+1.96*l0$se.fit)
## l1 <- loess(y~x,data=data.frame(y=y[x>=c & sample],x=x[x>=c & sample]))
## l1 <- predict(l1,se=T)
## l1 <- data.frame(Est=l1$fit, LI=l1$fit-1.96*l1$se.fit,
##                  UI=l1$fit+1.96*l1$se.fit)
## ###
## png("Figures/gpa.png")
## matplot(x[sample],rbind(l0,l1),type="n",
##         ylab="2nd Year GPA",
##         xlab="1st Year GPA",bty="n")
## matlines(sort(x[x<c & sample]),l0[order(x[x<c & sample]),],col="black",lty=1,lwd=c(1.5,1,1))
## matlines(sort(x[x>=c & sample]),l1[order(x[x>=c & sample]),],col="black",lty=1,lwd=c(1.5,1,1))
## abline(v=c,lty=2)
## dev.off()
## ## Descriptive statistics
## ### Full sample
## stats <- apply(w,2,function(x) c(mean(x),sd(x)))
## stats <- t(stats)
## rownames(stats) <- c("High school grade percentile",
##                      "Credits attempted in first year",
##                      "Age at entry",
##                      "Male",
##                      "Born in North America",
##                      "Campus 1","Campus 2","Campus 3")
## colnames(stats) <- c("Mean","Std. Dev")
## stats <- rbind(`Next Term GPA`=c(mean(y),sd(y)),
##                `Distance from cutoff`=c(mean(x),sd(x)),
##                `Treatment assignment`= c(mean(x>0),sd(x>0)),
##                stats)
## stats <- round(stats,2)
## ### XBCF window
## sample <- -0.1<=x & x<=0.1
## stats.xbcf <- apply(w[sample,],2,function(x) c(mean(x),sd(x)))
## stats.xbcf <- t(stats.xbcf)
## rownames(stats.xbcf) <- c("High school grade percentile",
##                      "Credits attempted in first year",
##                      "Age at entry",
##                      "Male",
##                      "Born in North America",
##                      "Campus 1","Campus 2","Campus 3")
## colnames(stats.xbcf) <- c("Mean","Std. Dev")
## stats.xbcf <- rbind(`Next Term GPA`=c(mean(y[sample]),sd(y[sample])),
##                `Distance from cutoff`=c(mean(x[sample]),sd(x[sample])),
##                `Treatment assignment`= c(mean(x[sample]>0),sd(x[sample]>0)),
##                stats.xbcf)
## stats.xbcf <- round(stats.xbcf,2)
## ### CKT window
## bw <- rdbwselect(y,x,c=0,covs=w)$bws[1]
## ckt.sample <- -bw<=x & x<=bw
## stats.ckt <- apply(w[ckt.sample,],2,function(x) c(mean(x),sd(x)))
## stats.ckt <- t(stats.ckt)
## rownames(stats.ckt) <- c("High school grade percentile",
##                           "Credits attempted in first year",
##                           "Age at entry",
##                           "Male",
##                           "Born in North America",
##                           "Campus 1","Campus 2","Campus 3")
## colnames(stats.ckt) <- c("Mean","Std. Dev")
## stats.ckt <- rbind(`Next Term GPA`=c(mean(y[ckt.sample]),sd(y[ckt.sample])),
##                `Distance from cutoff`=c(mean(x[ckt.sample]),sd(x[ckt.sample])),
##                `Treatment assignment`= c(mean(x[ckt.sample]>0),sd(x[ckt.sample]>0)),
##                stats.ckt)
## stats.ckt <- round(stats.ckt,2)
## ### Combining tables
## desc.stats <- cbind(stats,stats.xbcf,stats.ckt)
## desc.stats <- cbind(desc.stats[,1:2],"",
##                     desc.stats[,3:4],"",desc.stats[,5:6])
## print(xtable(desc.stats,caption="Descriptive statistics",
##              label="tab:desc.stats",align=c("l",rep("c",8))))
## ## Plot high school grade percentile against score
## png("Figures/grades_score.png")
## plot(x,w$hsgrade_pct,pch=21,bg="azure",main="",
##      xlab="1st year GPA",ylab="High school grade percentile",
##      cex=0.5,bty="n")
## dev.off()
## ## Plot high school grade percentile distribution
## d1 <- density(w$hsgrade_pct)
## d2 <- density(w[sample,]$hsgrade_pct)
## png("Figures/grades_dist.png")
## plot(d1,type="n",ylim=c(min(d1$y,d2$y),max(d1$y,d2$y)),bty="n",
##      xlab="High school grade percentile",ylab="",main="")
## polygon(d1,border="blue",col=rgb(0,0,1,0.25))
## polygon(d2,border="green",col=rgb(0,1,0,0.25))
## legend("topright",col=c("blue","green"),legend=c("Full","XBCF"),
##        title="Sample",lty=1)
## dev.off()
## ## 1st Year GPA distribution
## tab <- c(tail(table(x[x<=0]),5),head(table(x[x>=0]),5))
## names(tab)[1:4] <- substr(names(tab)[1:4],1,6)
## names(tab)[6:10] <- substr(names(tab)[6:10],1,5)
## png("Figures/gpa_cutoff_dist.png")
## barplot(tab,cex.names=0.75,cex.axis=0.75)
## dev.off()
## RDRobust estimation
cct1 <- rdrobust(y,x,c,all=T)
cct2 <- rdrobust(y,x,c,covs=w,all=T)
## XBCF estimation
fit.xbcf <- function(h,y,w,x,p_cat,pred.band,Omin)
{
    t0 <- Sys.time()
    fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
                   num_trees_mod = m, num_trees_con = m,
                   num_cutpoints = n, num_sweeps = num_sweeps,
                   burnin = burnin, Nmin = Nmin,
                   p_categorical_con = p_cat, p_categorical_mod = p_cat,
                   tau_con = 2*var(y)/m,
                   tau_mod = 0.5*var(y)/m, parallel=F)
    test <- -pred.band <= x & x<= pred.band
    pred <- predict.XBCFrd(fit,w[test,],rep(0,sum(test)))
    pred <- pred$tau.adj[,(burnin+1):num_sweeps]
    post <- colMeans(pred,na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    print(paste0("h=",h,"; Omin=",Omin,"; ATE=",mean(post)))
    return(list(ate.post=post,pred=pred,Owidth=h,time=dt))
}
###
n             <- length(y)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 150
burnin        <- 50
p_categorical <- ncol(w)
num_cutpoints <- n
###
### Parallelization
num_cores <- detectCores() - 1
registerDoParallel(num_cores)
## ### With controls
## fit.cont <- function(s)
## {
##     foreach(i=0:s,.multicombine=T,.export=c("y","w","x","p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
##         {
##             if (i<3)
##             {
##                 print(paste0("Model: h=",0.08+i/100,", Omin=",1))
##                 fit.xbcf(0.08+i/100,y,w,x,p_categorical,0.1,1)
##             } else
##             {
##                 print(paste0("Model: h=",0.08+(i-3)/100,", Omin=",5))
##                 fit.xbcf(0.08+(i-3)/100,y,w,x,p_categorical,0.1,5)
##             }
##         }
## }
## xbcf.cont <- fit.cont(5)
## ### Saving results to avoid running it all again
## save(xbcf.cont,file="Tables/application1.RData")
## rm(xbcf.cont)
## gc()
## print("Finished estimations")
load("Tables/application1.RData")
ate.post.cont1 <- xbcf.cont[[1]]$ate.post
ate.post.cont2 <- xbcf.cont[[2]]$ate.post
ate.post.cont3 <- xbcf.cont[[3]]$ate.post
ate.post.cont4 <- xbcf.cont[[4]]$ate.post
ate.post.cont5 <- xbcf.cont[[5]]$ate.post
ate.post.cont6 <- xbcf.cont[[6]]$ate.post
## ### Without controls
## fit.no.cont <- function(s)
## {
##     foreach(i=0:s,.multicombine=T,.export=c("y","w","x","p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
##         {
##             if (i<3)
##             {
##                 print(paste0("Model: h=",0.08+i/100,", Omin=",1))
##                 fit.xbcf(0.08+i/100,y,NULL,x,0,0.1,1)
##             } else
##             {
##                 print(paste0("Model: h=",0.08+(i-3)/100,", Omin=",5))
##                 fit.xbcf(0.08+(i-3)/100,y,NULL,x,0,0.1,5)
##             }
##         }
## }
## xbcf.no.cont <- fit.no.cont(5)
## ### Saving results to avoid running it all again
## save(xbcf.no.cont,file="Tables/application2.RData")
## rm(xbcf.no.cont)
## gc()
## print("Finished estimation")
load("Tables/application2.RData")
ate.post.no.cont1 <- xbcf.no.cont[[1]]$ate.post
ate.post.no.cont2 <- xbcf.no.cont[[2]]$ate.post
ate.post.no.cont3 <- xbcf.no.cont[[3]]$ate.post
ate.post.no.cont4 <- xbcf.no.cont[[4]]$ate.post
ate.post.no.cont5 <- xbcf.no.cont[[5]]$ate.post
ate.post.no.cont6 <- xbcf.no.cont[[6]]$ate.post
## Main results
r1 <- paste(round(quantile(ate.post.no.cont3,c(.025,.975)),2),collapse=",")
r2 <- paste(round(quantile(ate.post.cont3,c(.025,.975)),2),collapse=",")
r3 <- paste(round(cct1$ci[3,],2),collapse=",")
r4 <- paste(round(cct2$ci[3,],2),collapse=",")
test <- -0.1<=x & x<=0.1
###
results <- matrix(NA,5,6)
results[1,] <- c("","Controls","hat{tau}","95% CI","h","N")
results[2,] <- c("XBCF","No",round(mean(ate.post.no.cont1),2),
                 paste0("[",r1,"]"),
                 0.1,sum(test))
results[3,] <- c("","Yes",round(mean(ate.post.cont1),2),
                 paste0("[",r2,"]"),
                 0.1,sum(test))
results[4,] <- c("CCT","No",round(cct1$coef[1],2),
                 paste0("[",r3,"]"),
                 round(cct1$bws[1,1],2),
                 sum(cct1$N_h))
results[5,] <- c("","Yes",round(cct2$coef[1],2),
                 paste0("[",r4,"]"),
                 round(cct2$bws[1,1],2),
                 sum(cct2$N_h))
print(xtable(results,caption="RD Estimates",label="tab:gpa.res"),
      include.rownames=F,include.colnames=F)
## Alternative specifications for BART-RDD
r1 <- paste(round(quantile(ate.post.no.cont1,c(.025,.975)),2),collapse=",")
r2 <- paste(round(quantile(ate.post.no.cont2,c(.025,.975)),2),collapse=",")
r3 <- paste(round(quantile(ate.post.no.cont4,c(.025,.975)),2),collapse=",")
r4 <- paste(round(quantile(ate.post.no.cont5,c(.025,.975)),2),collapse=",")
r5 <- paste(round(quantile(ate.post.no.cont6,c(.025,.975)),2),collapse=",")
r6 <- paste(round(quantile(ate.post.cont1,c(.025,.975)),2),collapse=",")
r7 <- paste(round(quantile(ate.post.cont2,c(.025,.975)),2),collapse=",")
r8 <- paste(round(quantile(ate.post.cont4,c(.025,.975)),2),collapse=",")
r9 <- paste(round(quantile(ate.post.cont5,c(.025,.975)),2),collapse=",")
r10 <- paste(round(quantile(ate.post.cont6,c(.025,.975)),2),collapse=",")
###
results.check <- matrix(NA,11,5)
results.check[1,] <- c("Controls","hat{tau}","95% CI","h","N")
results.check[2,] <- c("No",round(mean(ate.post.no.cont1),2),
                 paste0("[",r1,"]"),
                 0.08,1)
results.check[3,] <- c("No",round(mean(ate.post.no.cont2),2),
                 paste0("[",r2,"]"),
                 0.09,1)
results.check[4,] <- c("No",round(mean(ate.post.no.cont4),2),
                 paste0("[",r3,"]"),
                 0.08,5)
results.check[5,] <- c("No",round(mean(ate.post.no.cont5),2),
                 paste0("[",r4,"]"),
                 0.09,5)
results.check[6,] <- c("No",round(mean(ate.post.no.cont6),2),
                 paste0("[",r5,"]"),
                 0.10,5)
results.check[7,] <- c("Yes",round(mean(ate.post.cont1),2),
                 paste0("[",r6,"]"),
                 0.08,1)
results.check[8,] <- c("Yes",round(mean(ate.post.cont2),2),
                 paste0("[",r7,"]"),
                 0.09,1)
results.check[9,] <- c("Yes",round(mean(ate.post.cont4),2),
                 paste0("[",r8,"]"),
                 0.08,5)
results.check[10,] <- c("Yes",round(mean(ate.post.cont5),2),
                 paste0("[",r9,"]"),
                 0.09,5)
results.check[11,] <- c("Yes",round(mean(ate.post.cont6),2),
                 paste0("[",r10,"]"),
                 0.10,5)
results.check <- results.check[,c(1,5,4,2,3)]
print(xtable(results.check,caption="RD Estimates - Sensitivity analysis",label="tab:gpa.robust"),
      include.rownames=F,include.colnames=F)
## CATE
cart1 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[1]]$pred),
                              w[test,]))
cart2 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[2]]$pred),
                              w[test,]))
cart3 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[3]]$pred),
                              w[test,]))
cart4 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[4]]$pred),
                              w[test,]))
cart5 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[5]]$pred),
                              w[test,]))
cart6 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[6]]$pred),
                              w[test,]))
### Plots
cart <- list(cart1,cart2,cart3,cart4,cart5,cart6)
for (i in 1:length(cart))
{
    pdf(paste0("Figures/gpa_cart_",i,".pdf"))
    rpart.plot(cart[[i]])
    dev.off()
}
#### Outcomes per age, grade and campus
bp.age <- data.frame(y=y,Treated=x>=0,Age=w$age_at_entry>=19)[test,]
boxplot(bp.age$y~bp.age$Treated+bp.age$Age,col=c("red","blue"),
        axes=F,xlab="",ylab="")
axis(1,c(1.5,3.5),
     labels=c("Below 19 years old","Above 19 years old"))
axis(2)
legend("topright",legend=c("Treatment","Control"),
       fill=c("blue","red"))
#####
bp.camp <- data.frame(y=y,Treated=x>=0,Camp=ifelse(w$loc_campus1==1,1,ifelse(w$loc_campus2,2,3)))[test,]
boxplot(bp.camp$y~bp.camp$Treated+bp.camp$Camp,col=c("red","blue"),
        axes=F,xlab="",ylab="")
axis(1,c(1.5,3.5,5.5),
     labels=c("Campus 1","Campus 2","Campus 3"))
axis(2)
legend("topright",legend=c("Treatment","Control"),
       fill=c("blue","red"))
#####
bp.hsg <- data.frame(y=y,Treated=x>=0,Grade=w$hsgrade_pct>=34)[test,]
boxplot(bp.hsg$y~bp.hsg$Treated+bp.hsg$Grade,col=c("red","blue"),
        axes=F,xlab="",ylab="")
axis(1,c(1.5,3.5),labels=c("Below 34","Above 34"))
axis(2)
legend("topright",legend=c("Treatment","Control"),
       fill=c("blue","red"))
#### Posterior for subgroups
s1 <- w[test,]$age_at_entry<19 & w[test,]$loc_campus3==1
s2 <- w[test,]$age_at_entry<19 & w[test,]$loc_campus3==0
s3 <- w[test,]$age_at_entry<19 & w[test,]$hsgrade_pct<34 & w[test,]$loc_campus3==1
s4 <- w[test,]$age_at_entry<19 & w[test,]$hsgrade_pct<34 & w[test,]$loc_campus3==0
#####
p1 <- colMeans(xbcf.cont[[3]]$pred[s1,])
p2 <- colMeans(xbcf.cont[[3]]$pred[s2,])
p3 <- colMeans(xbcf.cont[[3]]$pred[s3,])
p4 <- colMeans(xbcf.cont[[3]]$pred[s4,])
#####
png("Figures/post_age.png")
boxplot(by(xbcf.cont[[2]]$pred,w[test,]$age_at_entry,colMeans),axes=F,xlab="Age at entry",ylab=expression(tau))
axis(1,at=1:5,labels=17:21)
axis(2)
dev.off()
##
grade.sd <- by(xbcf.cont[[3]]$pred,w[test,]$hsgrade_pct,function(x) sd(colMeans(x)))
grade.mean <- by(xbcf.cont[[3]]$pred,w[test,]$hsgrade_pct,function(x) mean(colMeans(x)))
grade.sd.u <- grade.mean + grade.sd
grade.sd.l <- grade.mean - grade.sd
grade.plot <- cbind(grade.sd.l,grade.mean,grade.sd.u)
png("Figures/post_grade.png")
matplot(grade.plot,type=c("l","b","l"),lty=c(2,1,2),pch=19,axes=F,col="black",xlab="High school grade percentile",ylab=expression(tau))
axis(1,at=c(1,12,25,38,50,62,74),labels=rownames(grade.plot)[c(1,12,25,38,50,62,74)])
axis(2)
dev.off()
##
den <- density(p1-p2)
png("Figures/post_diff_camp_1.png")
plot(den,bty="n",xlab=expression(Delta),ylab="Density",main="")
polygon(c(den$x[den$x>=0],0),c(den$y[den$x>=0],0),col="black",
        density=25,angle=45)
dev.off()
den <- density(p3-p4)
png("Figures/post_diff_camp_2.png")
plot(den,bty="n",xlab=expression(Delta),ylab="Density",main="")
polygon(c(den$x[den$x>=0],0),c(den$y[den$x>=0],0),col="black",
        density=25,angle=45)
dev.off()
###
stopImplicitCluster()
