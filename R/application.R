## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("rafaelcalcantara/XBART@XBCF-RDD")
library(XBART)
library(rdrobust)
library(rpart)
library(rpart.plot)
library(xtable)
###
setwd("~/Documents/Git/XBCF-RDD/")
## data <- readRDS("Application/cost_sharing.rds")
## list2env(data,globalenv())
data <- read.csv("Application/gpa.csv")
y <- data$nextGPA
x <- data$X
w <- data[,4:11]
c <- 0
sample <- -0.3<=x & x<=0.3 ## For XBCF and loess plot
### Original data is children-level information, here it is
#### aggregated per distance to the cutoff, so variables are
#### averages for children t days away from the cutoff
## Plot data
png("Figures/gpa_data.png")
plot(x,y,pch=21,bg="azure",cex=0.5,xlab="1st Year GPA",ylab="2nd Year GPA")
dev.off()
###
l0 <- loess(y~x,data=data.frame(y=y[x<c & sample],x=x[x<c & sample]))
l0 <- predict(l0,se=T)
l0 <- data.frame(Est=l0$fit, LI=l0$fit-1.96*l0$se.fit,
                 UI=l0$fit+1.96*l0$se.fit)
l1 <- loess(y~x,data=data.frame(y=y[x>=c & sample],x=x[x>=c & sample]))
l1 <- predict(l1,se=T)
l1 <- data.frame(Est=l1$fit, LI=l1$fit-1.96*l1$se.fit,
                 UI=l1$fit+1.96*l1$se.fit)
###
png("Figures/gpa.png")
matplot(x[sample],rbind(l0,l1),type="n",
        ylab="2nd Year GPA",
        xlab="1st Year GPA")
matlines(sort(x[x<c & sample]),l0[order(x[x<c & sample]),],col="black",lty=1,lwd=c(1.5,1,1))
matlines(sort(x[x>=c & sample]),l1[order(x[x>=c & sample]),],col="black",lty=1,lwd=c(1.5,1,1))
abline(v=c,lty=2)
dev.off()
## Descriptive statistics
### Full sample
stats <- apply(w,2,function(x) c(mean(x),sd(x)))
stats <- t(stats)
rownames(stats) <- c("High school grade percentile",
                     "Credits attempted in first year",
                     "Age at entry",
                     "Male",
                     "Born in North America",
                     "Campus 1","Campus 2","Campus 3")
colnames(stats) <- c("Mean","Std. Dev")
stats <- rbind(`Next Term GPA`=c(mean(y),sd(y)),
               `Distance from cutoff`=c(mean(x),sd(x)),
               `Treatment assignment`= c(mean(x>0),sd(x>0)),
               stats)
stats <- round(stats,2)
### XBCF subsample
stats.subsample <- apply(w[sample,],2,function(x) c(mean(x),sd(x)))
stats.subsample <- t(stats.subsample)
rownames(stats.subsample) <- c("High school grade percentile",
                     "Credits attempted in first year",
                     "Age at entry",
                     "Male",
                     "Born in North America",
                     "Campus 1","Campus 2","Campus 3")
colnames(stats.subsample) <- c("Mean","Std. Dev")
stats.subsample <- rbind(`Next Term GPA`=c(mean(y[sample]),sd(y[sample])),
               `Distance from cutoff`=c(mean(x[sample]),sd(x[sample])),
               `Treatment assignment`= c(mean(x[sample]>0),sd(x[sample]>0)),
               stats.subsample)
stats.subsample <- round(stats.subsample,2)
### CKT window
bw <- rdbwselect(y,x,c=0,covs=w)$bws[1]
ckt.sample <- -bw<=x & x<=bw
stats.ckt <- apply(w[ckt.sample,],2,function(x) c(mean(x),sd(x)))
stats.ckt <- t(stats.ckt)
rownames(stats.ckt) <- c("High school grade percentile",
                          "Credits attempted in first year",
                          "Age at entry",
                          "Male",
                          "Born in North America",
                          "Campus 1","Campus 2","Campus 3")
colnames(stats.ckt) <- c("Mean","Std. Dev")
stats.ckt <- rbind(`Next Term GPA`=c(mean(y[ckt.sample]),sd(y[ckt.sample])),
               `Distance from cutoff`=c(mean(x[ckt.sample]),sd(x[ckt.sample])),
               `Treatment assignment`= c(mean(x[ckt.sample]>0),sd(x[ckt.sample]>0)),
               stats.ckt)
stats.ckt <- round(stats.ckt,2)
### Combining tables
desc.stats <- cbind(stats,stats.subsample,stats.ckt)
desc.stats <- cbind(desc.stats[,1:2],"",
                    desc.stats[,3:4],"",desc.stats[,5:6])
print(xtable(desc.stats,caption="Descriptive statistics",
             label="tab:desc.stats",align=c("l",rep("c",8))))
## Plot high school grade percentile against score
png("Figures/grades_score.png")
plot(x,w$hsgrade_pct,pch=21,bg="azure",main="",
     xlab="1st year GPA",ylab="High school grade percentile",
     cex=0.5,bty="n")
dev.off()
## Plot high school grade percentile distribution
d1 <- density(w$hsgrade_pct)
d2 <- density(w[sample,]$hsgrade_pct)
png("Figures/grades_dist.png")
plot(d1,type="n",ylim=c(min(d1$y,d2$y),max(d1$y,d2$y)),bty="n",
     xlab="High school grade percentile",ylab="",main="")
polygon(d1,border="blue",col=rgb(0,0,1,0.25))
polygon(d2,border="green",col=rgb(0,1,0,0.25))
legend("topright",col=c("blue","green"),legend=c("Full","XBCF"),
       title="Sample",lty=1)
dev.off()
## RDRobust estimation
cct1 <- rdrobust(y,x,c,all=T)
cct2 <- rdrobust(y,x,c,covs=w,all=T)
## XBCF estimation
fit.xbcf <- function(h,y,w,x,p_cat)
{
    t0 <- Sys.time()
    fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
                           num_trees_mod = m, num_trees_con = m,
                           num_cutpoints = n, num_sweeps = num_sweeps,
                           burnin = burnin, Nmin = Nmin,
                           p_categorical_con = p_cat, p_categorical_mod = p_cat,
                           tau_con = 2*var(y)/m,
                   tau_mod = 0.5*var(y)/m, parallel=T,
                   nthread = no_cores)
    pred <- predict.XBCFrd(fit,w,rep(0,n))
    post <- colMeans(pred$tau.adj[,(burnin+1):num_sweeps],na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(ate.post=post,pred=fit,Owidth=h,time=dt))
}
###
n             <- sum(sample)
Omin          <- as.integer(0.03*n)
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 120
burnin        <- 20
p_categorical <- ncol(w)
num_cutpoints <- n
###
### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
###
xbcf1 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],NULL,x[sample],0)
ate.post1 <- xbcf1$ate.post
xbcf2 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post2 <- xbcf2$ate.post
## Comparing estimates
r1 <- paste(round(quantile(ate.post1,c(.025,.975)),2),collapse=",")
r2 <- paste(round(quantile(ate.post2,c(.025,.975)),2),collapse=",")
r3 <- paste(round(cct1$ci[3,],2),collapse=",")
r4 <- paste(round(cct2$ci[3,],2),collapse=",")
###
results <- matrix(NA,5,6)
results[1,] <- c("","Controls","hat{tau}","95% CI","h","N")
results[2,] <- c("XBCF","No",round(mean(ate.post1),2),
                 paste0("[",r1,"]"),
                 round(xbcf1$Owidth,2),
                 n)
results[3,] <- c("","Yes",round(mean(ate.post2),2),
                 paste0("[",r2,"]"),
                 round(xbcf2$Owidth,2),
                 n)
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
## Robustness
### Different windows
xbcf3 <- fit.xbcf(quantile(abs(x[sample]),0.1,na.rm=T),y[sample],NULL,x[sample],0)
ate.post3 <- xbcf3$ate.post
xbcf4 <- fit.xbcf(quantile(abs(x[sample]),0.1,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post4 <- xbcf4$ate.post
xbcf5 <- fit.xbcf(quantile(abs(x[sample]),0.15,na.rm=T),y[sample],NULL,x[sample],0)
ate.post5 <- xbcf5$ate.post
xbcf6 <- fit.xbcf(quantile(abs(x[sample]),0.15,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post6 <- xbcf6$ate.post
xbcf7 <- fit.xbcf(quantile(abs(x[sample]),0.2,na.rm=T),y[sample],NULL,x[sample],0)
ate.post7 <- xbcf7$ate.post
xbcf8 <- fit.xbcf(quantile(abs(x[sample]),0.2,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post8 <- xbcf8$ate.post
###
r3 <- paste(round(quantile(ate.post3,c(.025,.975)),3),collapse=",")
r4 <- paste(round(quantile(ate.post4,c(.045,.975)),3),collapse=",")
r5 <- paste(round(quantile(ate.post5,c(.025,.975)),3),collapse=",")
r6 <- paste(round(quantile(ate.post6,c(.065,.975)),3),collapse=",")
r7 <- paste(round(quantile(ate.post7,c(.025,.975)),3),collapse=",")
r8 <- paste(round(quantile(ate.post8,c(.085,.975)),3),collapse=",")
###
robustness <- matrix(NA,7,4)
robustness[1,] <- c("Controls","hat{tau}","95% CI","h")
robustness[2,] <- c("No",round(mean(ate.post3),2),
                    paste0("[",r3,"]"),
                    round(xbcf3$Owidth,2))
robustness[3,] <- c("Yes",round(mean(ate.post4),2),
                    paste0("[",r4,"]"),
                    round(xbcf4$Owidth,2))
robustness[4,] <- c("No",round(mean(ate.post5),2),
                    paste0("[",r5,"]"),
                    round(xbcf5$Owidth,2))
robustness[5,] <- c("Yes",round(mean(ate.post6),2),
                    paste0("[",r6,"]"),
                    round(xbcf6$Owidth,2))
robustness[6,] <- c("No",round(mean(ate.post7),2),
                    paste0("[",r7,"]"),
                    round(xbcf7$Owidth,2))
robustness[7,] <- c("Yes",round(mean(ate.post8),2),
                    paste0("[",r8,"]"),
                    round(xbcf8$Owidth,2))
print(xtable(robustness,caption="XBCF Estimates - Robustness to h",label="tab:gpa.robust.h"),
      include.rownames=F,include.colnames=F)
### Different subsamples
sample        <- -1<=x & x<=1
n             <- sum(sample)
n3            <- n
Omin          <- as.integer(0.03*n)
num_cutpoints <- n
xbcf3 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],NULL,x[sample],0)
ate.post3 <- xbcf3$ate.post
xbcf4 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post4 <- xbcf4$ate.post
####
sample        <- -0.5<=x & x<=0.5
n             <- sum(sample)
n5            <- n
Omin          <- as.integer(0.03*n)
num_cutpoints <- n
xbcf5 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],NULL,x[sample],0)
ate.post5 <- xbcf5$ate.post
xbcf6 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post6 <- xbcf6$ate.post
####
sample        <- -0.1<=x & x<=0.1
n             <- sum(sample)
n7            <- n
Omin          <- as.integer(0.03*n)
num_cutpoints <- n
xbcf7 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],NULL,x[sample],0)
ate.post7 <- xbcf7$ate.post
xbcf8 <- fit.xbcf(quantile(abs(x[sample]),0.125,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post8 <- xbcf8$ate.post
###
r3 <- paste(round(quantile(ate.post3,c(.025,.975)),3),collapse=",")
r4 <- paste(round(quantile(ate.post4,c(.045,.975)),3),collapse=",")
r5 <- paste(round(quantile(ate.post5,c(.025,.975)),3),collapse=",")
r6 <- paste(round(quantile(ate.post6,c(.065,.975)),3),collapse=",")
r7 <- paste(round(quantile(ate.post7,c(.025,.975)),3),collapse=",")
r8 <- paste(round(quantile(ate.post8,c(.085,.975)),3),collapse=",")
###
robustness <- matrix(NA,7,5)
robustness[1,] <- c("Controls","hat{tau}","95% CI","h","N")
robustness[2,] <- c("No",round(mean(ate.post3),2),
                    paste0("[",r3,"]"),
                    round(xbcf3$Owidth,2),n3)
robustness[3,] <- c("Yes",round(mean(ate.post4),2),
                    paste0("[",r4,"]"),
                    round(xbcf4$Owidth,2),n3)
robustness[4,] <- c("No",round(mean(ate.post5),2),
                    paste0("[",r5,"]"),
                    round(xbcf5$Owidth,2),n5)
robustness[5,] <- c("Yes",round(mean(ate.post6),2),
                    paste0("[",r6,"]"),
                    round(xbcf6$Owidth,2),n5)
robustness[6,] <- c("No",round(mean(ate.post7),2),
                    paste0("[",r7,"]"),
                    round(xbcf7$Owidth,2),n7)
robustness[7,] <- c("Yes",round(mean(ate.post8),2),
                    paste0("[",r8,"]"),
                    round(xbcf8$Owidth,2),n7)
print(xtable(robustness,caption="XBCF Estimates - Robustness to sample choice",label="tab:gpa.robust.n"),
      include.rownames=F,include.colnames=F)
## Heterogeneous effects
num_sweeps <- 1020
sample        <- -0.3<=x & x<=0.3
n             <- sum(sample)
Omin          <- as.integer(0.03*n)
num_cutpoints <- n
xbcf.het <- fit.xbcf(quantile(abs(x[sample]),0.2,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
pred <- predict.XBCFrd(xbcf.het$pred,w[sample,],rep(0,sum(sample)))
pred <- pred$tau.adj[,(burnin+1):num_sweeps]
summary(colMeans(pred))
xbcf.het$pred$importance_treatment
cart <- rpart(y~.,data.frame(y=rowMeans(pred),w[sample,]))
### Plotting
png("Figures/cart_tau.png")
rpart.plot(cart)
dev.off()
### Subgroups
s1 <- w[sample,]$hsgrade_pct >= 35
s2 <- w[sample,]$hsgrade_pct >= 35 & w[sample,]$totcredits_year1 >= 4.8
s3 <- w[sample,]$hsgrade_pct < 35 & w[sample,]$totcredits_year1 >= 4.8
s4 <- w[sample,]$hsgrade_pct >= 35 & w[sample,]$totcredits_year1 < 4.8
s5 <- w[sample,]$hsgrade_pct < 35 & w[sample,]$totcredits_year1 < 4.8
s6 <- w[sample,]$hsgrade_pct >= 35 & w[sample,]$male==1
s7 <- w[sample,]$hsgrade_pct < 35 & w[sample,]$male==1
s6 <- w[sample,]$hsgrade_pct >= 35 & w[sample,]$male==0
s7 <- w[sample,]$hsgrade_pct < 35 & w[sample,]$male==0
###
sapply(by(pred,s1,colMeans),function(i) c(mean(i),quantile(i,c(0.025,0.975))))
sapply(by(pred,s2,colMeans),function(i) c(mean(i),quantile(i,c(0.025,0.975))))
sapply(by(pred,s3,colMeans),function(i) c(mean(i),quantile(i,c(0.025,0.975))))
sapply(by(pred,s4,colMeans),function(i) c(mean(i),quantile(i,c(0.025,0.975))))
### Grade pct < 35 seems to be the only one making a big difference
png("Figures/cate_posterior.png")
plot(density(colMeans(pred[s1,])),
     xlab="", ylab="", main="",col="red",axes=F)
axis(1)
axis(2)
lines(density(colMeans(pred[!s1,])),col="blue")
polygon(density(colMeans(pred[s1,])),col=rgb(1,0,0,0.25))
polygon(density(colMeans(pred[!s1,])),col=rgb(0,0,1,0.25))
legend("topleft",col=c("red","blue"),lty=1,
       legend=c("Above 35","Below 35"),
       title="High School\nGrade Percentile",cex=0.75)
abline(v=0,lty=2)
dev.off()
###
den <- density(colMeans(pred[!s1,])-colMeans(pred[s1,]))
png("Figures/cate_difference.png")
plot(den,bty="n",xlab=expression(Delta),ylab="Density",main="")
polygon(c(den$x[den$x>=0],0),c(den$y[den$x>=0],0),col="black",
        density=25,angle=45)
dev.off()
###
stopImplicitCluster()
