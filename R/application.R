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
plot(x,y,pch=21,bg="azure",cex=0.5,xlab="1st Year GPA",ylab="2nd Year GPA",bty="n")
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
        xlab="1st Year GPA",bty="n")
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
### XBCF window
sample <- -0.01<=x & x<=0.01
stats.xbcf <- apply(w[sample,],2,function(x) c(mean(x),sd(x)))
stats.xbcf <- t(stats.xbcf)
rownames(stats.xbcf) <- c("High school grade percentile",
                     "Credits attempted in first year",
                     "Age at entry",
                     "Male",
                     "Born in North America",
                     "Campus 1","Campus 2","Campus 3")
colnames(stats.xbcf) <- c("Mean","Std. Dev")
stats.xbcf <- rbind(`Next Term GPA`=c(mean(y[sample]),sd(y[sample])),
               `Distance from cutoff`=c(mean(x[sample]),sd(x[sample])),
               `Treatment assignment`= c(mean(x[sample]>0),sd(x[sample]>0)),
               stats.xbcf)
stats.xbcf <- round(stats.xbcf,2)
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
desc.stats <- cbind(stats,stats.xbcf,stats.ckt)
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
fit.xbcf <- function(h,y,w,x,p_cat,pred.band,Omin)
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
    test <- -pred.band <= x & x<= pred.band
    pred <- predict.XBCFrd(fit,w[test,],rep(0,sum(test)))
    pred <- pred$tau.adj[,(burnin+1):num_sweeps]
    post <- colMeans(pred,na.rm=T)
    t1 <- Sys.time()
    dt <- difftime(t1,t0)
    print(paste0("Elapsed time: ",round(dt,2)," seconds"))
    return(list(fit=fit,ate.post=post,pred=pred,Owidth=h,time=dt))
}
###
n             <- length(y)
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
### With controls
xbcf.cont.1 <- fit.xbcf(0.01,y,w,x,p_categorical,0.01,10)
ate.post.cont1 <- xbcf.cont.1$ate.post
xbcf.cont.2 <- fit.xbcf(0.02,y,w,x,p_categorical,0.01,10)
ate.post.cont2 <- xbcf.cont.2$ate.post
xbcf.cont.3 <- fit.xbcf(0.03,y,w,x,p_categorical,0.01,10)
ate.post.cont3 <- xbcf.cont.3$ate.post
xbcf.cont.4 <- fit.xbcf(0.04,y,w,x,p_categorical,0.01,10)
ate.post.cont4 <- xbcf.cont.4$ate.post
xbcf.cont.5 <- fit.xbcf(0.01,y,w,x,p_categorical,0.01,50)
ate.post.cont5 <- xbcf.cont.5$ate.post
xbcf.cont.6 <- fit.xbcf(0.02,y,w,x,p_categorical,0.01,50)
ate.post.cont6 <- xbcf.cont.6$ate.post
xbcf.cont.7 <- fit.xbcf(0.03,y,w,x,p_categorical,0.01,50)
ate.post.cont7 <- xbcf.cont.7$ate.post
xbcf.cont.8 <- fit.xbcf(0.04,y,w,x,p_categorical,0.01,50)
ate.post.cont8 <- xbcf.cont.8$ate.post
### Without controls
xbcf.no.cont.1 <- fit.xbcf(0.01,y,NULL,x,0,0.01,10)
ate.post.no.cont1 <- xbcf.no.cont.1$ate.post
xbcf.no.cont.2 <- fit.xbcf(0.02,y,NULL,x,0,0.01,10)
ate.post.no.cont2 <- xbcf.no.cont.2$ate.post
xbcf.no.cont.3 <- fit.xbcf(0.03,y,NULL,x,0,0.01,10)
ate.post.no.cont3 <- xbcf.no.cont.3$ate.post
xbcf.no.cont.4 <- fit.xbcf(0.04,y,NULL,x,0,0.01,10)
ate.post.no.cont4 <- xbcf.no.cont.4$ate.post
xbcf.no.cont.5 <- fit.xbcf(0.01,y,NULL,x,0,0.01,50)
ate.post.no.cont5 <- xbcf.no.cont.5$ate.post
xbcf.no.cont.6 <- fit.xbcf(0.02,y,NULL,x,0,0.01,50)
ate.post.no.cont6 <- xbcf.no.cont.6$ate.post
xbcf.no.cont.7 <- fit.xbcf(0.03,y,NULL,x,0,0.01,50)
ate.post.no.cont7 <- xbcf.no.cont.7$ate.post
xbcf.no.cont.8 <- fit.xbcf(0.04,y,NULL,x,0,0.01,50)
ate.post.no.cont8 <- xbcf.no.cont.8$ate.post
## Main results
r1 <- paste(round(quantile(ate.post.no.cont1,c(.025,.975)),2),collapse=",")
r2 <- paste(round(quantile(ate.post.cont1,c(.025,.975)),2),collapse=",")
r3 <- paste(round(cct1$ci[3,],2),collapse=",")
r4 <- paste(round(cct2$ci[3,],2),collapse=",")
###
results <- matrix(NA,5,6)
results[1,] <- c("","Controls","hat{tau}","95% CI","h","N")
results[2,] <- c("XBCF","No",round(mean(ate.post.no.cont1),2),
                 paste0("[",r1,"]"),
                 0.01,10)
results[3,] <- c("","Yes",round(mean(ate.post.cont1),2),
                 paste0("[",r2,"]"),
                 0.01,10)
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
## Alternative specifications for XBCF
r1 <- paste(round(quantile(ate.post.no.cont1,c(.025,.975)),2),collapse=",")
r2 <- paste(round(quantile(ate.post.no.cont2,c(.025,.975)),2),collapse=",")
r3 <- paste(round(quantile(ate.post.no.cont3,c(.025,.975)),2),collapse=",")
r4 <- paste(round(quantile(ate.post.no.cont4,c(.025,.975)),2),collapse=",")
r5 <- paste(round(quantile(ate.post.no.cont5,c(.025,.975)),2),collapse=",")
r6 <- paste(round(quantile(ate.post.no.cont6,c(.025,.975)),2),collapse=",")
r7 <- paste(round(quantile(ate.post.no.cont7,c(.025,.975)),2),collapse=",")
r8 <- paste(round(quantile(ate.post.no.cont8,c(.025,.975)),2),collapse=",")
r9 <- paste(round(quantile(ate.post.cont1,c(.025,.975)),2),collapse=",")
r10 <- paste(round(quantile(ate.post.cont2,c(.025,.975)),2),collapse=",")
r11 <- paste(round(quantile(ate.post.cont3,c(.025,.975)),2),collapse=",")
r12 <- paste(round(quantile(ate.post.cont4,c(.025,.975)),2),collapse=",")
r13 <- paste(round(quantile(ate.post.cont5,c(.025,.975)),2),collapse=",")
r14 <- paste(round(quantile(ate.post.cont6,c(.025,.975)),2),collapse=",")
r15 <- paste(round(quantile(ate.post.cont7,c(.025,.975)),2),collapse=",")
r16 <- paste(round(quantile(ate.post.cont8,c(.025,.975)),2),collapse=",")
###
results.check <- matrix(NA,17,5)
results.check[1,] <- c("Controls","hat{tau}","95% CI","h","N")
results.check[2,] <- c("No",round(mean(ate.post.no.cont1),2),
                 paste0("[",r1,"]"),
                 0.01,10)
results.check[3,] <- c("No",round(mean(ate.post.no.cont2),2),
                 paste0("[",r2,"]"),
                 0.02,10)
results.check[4,] <- c("No",round(mean(ate.post.no.cont3),2),
                 paste0("[",r3,"]"),
                 0.03,10)
results.check[5,] <- c("No",round(mean(ate.post.no.cont4),2),
                 paste0("[",r4,"]"),
                 0.04,10)
results.check[6,] <- c("No",round(mean(ate.post.no.cont5),2),
                 paste0("[",r5,"]"),
                 0.01,50)
results.check[7,] <- c("No",round(mean(ate.post.no.cont6),2),
                 paste0("[",r6,"]"),
                 0.02,50)
results.check[8,] <- c("No",round(mean(ate.post.no.cont7),2),
                 paste0("[",r7,"]"),
                 0.03,50)
results.check[9,] <- c("No",round(mean(ate.post.no.cont8),2),
                 paste0("[",r8,"]"),
                 0.04,50)
results.check[10,] <- c("Yes",round(mean(ate.post.cont1),2),
                 paste0("[",r9,"]"),
                 0.01,10)
results.check[11,] <- c("Yes",round(mean(ate.post.cont2),2),
                 paste0("[",r10,"]"),
                 0.02,10)
results.check[12,] <- c("Yes",round(mean(ate.post.cont3),2),
                 paste0("[",r11,"]"),
                 0.03,10)
results.check[13,] <- c("Yes",round(mean(ate.post.cont4),2),
                 paste0("[",r12,"]"),
                 0.04,10)
results.check[14,] <- c("Yes",round(mean(ate.post.cont5),2),
                 paste0("[",r13,"]"),
                 0.01,50)
results.check[15,] <- c("Yes",round(mean(ate.post.cont6),2),
                 paste0("[",r14,"]"),
                 0.02,50)
results.check[16,] <- c("Yes",round(mean(ate.post.cont7),2),
                 paste0("[",r15,"]"),
                 0.03,50)
results.check[17,] <- c("Yes",round(mean(ate.post.cont8),2),
                 paste0("[",r16,"]"),
                 0.04,50)
print(xtable(results.check,caption="RD Estimates - Sensitivity analysis",label="tab:gpa.robust"),
      include.rownames=F,include.colnames=F)
## CATE
test <- -0.01<=x & x<=0.01
cart1 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.1$pred),
                              w[test,]))
cart2 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.2$pred),
                              w[test,]))
cart3 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.3$pred),
                              w[test,]))
cart4 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.4$pred),
                              w[test,]))
cart5 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.5$pred),
                              w[test,]))
cart6 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.6$pred),
                              w[test,]))
cart7 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.7$pred),
                              w[test,]))
cart8 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont.8$pred),
                              w[test,]))
### Plots
cart <- list(cart1,cart2,cart3,cart4,cart5,cart6,cart7,cart8)
for (i in length(cart))
{
    png(paste0("Figures/gpa_cart_",i,".png"))
    rpart.plot(cart[[i]])
    dev.off()
}
###
stopImplicitCluster()
save.image("Tables/application.RData")
