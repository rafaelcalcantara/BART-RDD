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
plot(x,y,pch=21,bg="azure",cex=0.5)
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
    post <- colMeans(pred$tau.adj,na.rm=T)
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
## no_cores <- detectCores() - 1
no_cores <- 1
registerDoParallel(no_cores)
###
xbcf1 <- fit.xbcf(quantile(abs(x),0.125,na.rm=T),y[sample],NULL,x[sample],0)
ate.post1 <- xbcf1$ate.post
xbcf2 <- fit.xbcf(quantile(abs(x),0.125,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
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
## Heterogeneous effects
### Covariates are: share of male children, household income per capita, share of children born in Taipei, and birth year
pred <- predict.XBCFrd(xbcf2$pred,w[sample,],x[sample])
pred <- pred$tau.adj[,(burnin+1):num_sweeps]
cart <- rpart(y~.,data.frame(y=rowMeans(pred),w[sample,]),control=rpart.control(cp=0.01))
rpart.plot(cart)
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
xbcf3 <- fit.xbcf(quantile(abs(x[sample]),0.1,na.rm=T),y[sample],NULL,x[sample],0)
ate.post3 <- xbcf3$ate.post
xbcf4 <- fit.xbcf(quantile(abs(x[sample]),0.1,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post4 <- xbcf4$ate.post
####
sample        <- -0.5<=x & x<=0.5
n             <- sum(sample)
n5            <- n
Omin          <- as.integer(0.03*n)
num_cutpoints <- n
xbcf5 <- fit.xbcf(quantile(abs(x[sample]),0.15,na.rm=T),y[sample],NULL,x[sample],0)
ate.post5 <- xbcf5$ate.post
xbcf6 <- fit.xbcf(quantile(abs(x[sample]),0.15,na.rm=T),y[sample],w[sample,],x[sample],p_categorical)
ate.post6 <- xbcf6$ate.post
####
sample        <- -0.1<=x & x<=0.1
n             <- sum(sample)
n7            <- n
Omin          <- as.integer(0.03*n)
num_cutpoints <- n
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
###
stopImplicitCluster()
#### Old
## p <- data.frame(Method = c("CCT1","CCT2","XBCF1","XBCF2"),
##                 Est = c(cct1$coef[1],cct2$coef[1],
##                         mean(ate.post1),mean(ate.post2)),
##                 LI = c(cct1$ci[1,1],cct2$ci[1,1],
##                        quantile(ate.post1,0.025),
##                        quantile(ate.post2,0.025)),
##                 UI = c(cct1$ci[1,2],cct2$ci[1,2],
##                        quantile(ate.post1,0.975),
##                        quantile(ate.post2,0.975)),
##                 Position = c(1,3,5,7))
## p1 <- reshape(p,direction="long",idvar=c("Method","Position"),varying=c("Est","LI","UI"),times=c("Est","LI","UI"),v.name="value")
## cols <- c("firebrick","dodgerblue")
## plot(p1$Position,p1$value,"n",xlab="",ylab="",axes=F)
## axis(2)
## axis(1,at=c(1,2,6,8),labels=c("","CCT","XBCF",""))
## segments(x0=p$Position,y0=p$LI,y1=p$UI,col=cols,lwd=1.5)
## segments(x0=p$Position-0.05,x1=p$Position+0.05,y0=p$LI,col=cols,lwd=1.5)
## segments(x0=p$Position-0.05,x1=p$Position+0.05,y0=p$UI,col=cols,lwd=1.5)
## points(x=p$Position,y=p$Est,col=cols,pch=19)
## legend("topright",legend=c("Yes","No"),cex=0.7,
##        col=c(cols[2],cols[1]),lty=1,lwd=1.5,pch=19,
##        title="Includes covariates",horiz=T)
