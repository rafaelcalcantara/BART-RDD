## Setup
library(parallel)
library(foreach)
library(doParallel)
devtools::install_github("JingyuHe/XBART@XBCF-RDD")
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
sample <- -0.1<=x & x<=0.1
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
## 1st Year GPA distribution
tab <- c(tail(table(x[x<=0]),5),head(table(x[x>=0]),5))
names(tab)[1:4] <- substr(names(tab)[1:4],1,6)
names(tab)[6:10] <- substr(names(tab)[6:10],1,5)
png("Figures/gpa_cutoff_dist.png")
barplot(tab,cex.names=0.75,cex.axis=0.75)
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
##     foreach(i=1:s,.multicombine=T,.export=c("y","w","x","p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
##         {
##             if (i<5)
##             {
##                 print(paste0("Model: h=",0.08+i/100,", Omin=",10))
##                 fit.xbcf(0.08+i/100,y,w,x,p_categorical,0.1,10)
##             } else
##             {
##                 print(paste0("Model: h=",0.08+(i-4)/100,", Omin=",50))
##                 fit.xbcf(0.08+(i-4)/100,y,w,x,p_categorical,0.1,50)
##             }
##         }
## }
## xbcf.cont <- fit.cont(8)
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
ate.post.cont7 <- xbcf.cont[[7]]$ate.post
ate.post.cont8 <- xbcf.cont[[8]]$ate.post
## ### Without controls
## fit.no.cont <- function(s)
## {
##     foreach(i=1:s,.multicombine=T,.export=c("y","w","x","p","c","Omin","Opct","m","n","num_sweeps","burnin","Nmin","p_categorical")) %dopar%
##         {
##             if (i<5)
##             {
##                 print(paste0("Model: h=",0.08+i/100,", Omin=",10))
##                 fit.xbcf(0.08+i/100,y,NULL,x,0,0.1,10)
##             } else
##             {
##                 print(paste0("Model: h=",0.08+(i-4)/100,", Omin=",50))
##                 fit.xbcf(0.08+(i-4)/100,y,NULL,x,0,0.1,50)
##             }
##         }
## }
## xbcf.no.cont <- fit.no.cont(8)
## ### Saving results to avoid running it all again
## save(xbcf.no.cont,file="Tables/application2.RData")
## rm(xbcf.no.cont)
## gc()
load("Tables/application2.RData")
ate.post.no.cont1 <- xbcf.no.cont[[1]]$ate.post
ate.post.no.cont2 <- xbcf.no.cont[[2]]$ate.post
ate.post.no.cont3 <- xbcf.no.cont[[3]]$ate.post
ate.post.no.cont4 <- xbcf.no.cont[[4]]$ate.post
ate.post.no.cont5 <- xbcf.no.cont[[5]]$ate.post
ate.post.no.cont6 <- xbcf.no.cont[[6]]$ate.post
ate.post.no.cont7 <- xbcf.no.cont[[7]]$ate.post
ate.post.no.cont8 <- xbcf.no.cont[[8]]$ate.post
## Main results
r1 <- paste(round(quantile(ate.post.no.cont2,c(.025,.975)),2),collapse=",")
r2 <- paste(round(quantile(ate.post.cont2,c(.025,.975)),2),collapse=",")
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
sample <- sapply(seq(0.09,0.12,0.01),function(i) sum(-i<=x & x<=i))
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
                 0.09,10)
results.check[3,] <- c("No",round(mean(ate.post.no.cont2),2),
                 paste0("[",r2,"]"),
                 0.10,10)
results.check[4,] <- c("No",round(mean(ate.post.no.cont3),2),
                 paste0("[",r3,"]"),
                 0.11,10)
results.check[5,] <- c("No",round(mean(ate.post.no.cont4),2),
                 paste0("[",r4,"]"),
                 0.12,10)
results.check[6,] <- c("No",round(mean(ate.post.no.cont5),2),
                 paste0("[",r5,"]"),
                 0.09,50)
results.check[7,] <- c("No",round(mean(ate.post.no.cont6),2),
                 paste0("[",r6,"]"),
                 0.10,50)
results.check[8,] <- c("No",round(mean(ate.post.no.cont7),2),
                 paste0("[",r7,"]"),
                 0.11,50)
results.check[9,] <- c("No",round(mean(ate.post.no.cont8),2),
                 paste0("[",r8,"]"),
                 0.12,50)
results.check[10,] <- c("Yes",round(mean(ate.post.cont1),2),
                 paste0("[",r9,"]"),
                 0.09,10)
results.check[11,] <- c("Yes",round(mean(ate.post.cont2),2),
                 paste0("[",r10,"]"),
                 0.10,10)
results.check[12,] <- c("Yes",round(mean(ate.post.cont3),2),
                 paste0("[",r11,"]"),
                 0.11,10)
results.check[13,] <- c("Yes",round(mean(ate.post.cont4),2),
                 paste0("[",r12,"]"),
                 0.12,10)
results.check[14,] <- c("Yes",round(mean(ate.post.cont5),2),
                 paste0("[",r13,"]"),
                 0.09,50)
results.check[15,] <- c("Yes",round(mean(ate.post.cont6),2),
                 paste0("[",r14,"]"),
                 0.10,50)
results.check[16,] <- c("Yes",round(mean(ate.post.cont7),2),
                 paste0("[",r15,"]"),
                 0.11,50)
results.check[17,] <- c("Yes",round(mean(ate.post.cont8),2),
                 paste0("[",r16,"]"),
                 0.12,50)
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
cart7 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[7]]$pred),
                              w[test,]))
cart8 <- rpart(y~.,data.frame(y=rowMeans(xbcf.cont[[8]]$pred),
                              w[test,]))
### Plots
cart <- list(cart1,cart2,cart3,cart4,cart5,cart6,cart7,cart8)
for (i in 1:length(cart))
{
    png(paste0("Figures/gpa_cart_",i,".png"))
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
s3 <- w[test,]$age_at_entry<19 & w[test,]$hsgrade_pct<34 & w[test,]$loc_campus1==1
s4 <- w[test,]$age_at_entry<19 & w[test,]$hsgrade_pct<34 & w[test,]$loc_campus3==0
#####
p1 <- colMeans(xbcf.cont[[2]]$pred[s1,])
p2 <- colMeans(xbcf.cont[[2]]$pred[s2,])
p3 <- colMeans(xbcf.cont[[2]]$pred[s3,])
p4 <- colMeans(xbcf.cont[[2]]$pred[s4,])
#####
png("Figures/post_age.png")
boxplot(by(xbcf.cont[[2]]$pred,w[test,]$age_at_entry,colMeans),axes=F,xlab="Age at entry",ylab=expression(tau))
axis(1,at=1:5,labels=17:21)
axis(2)
dev.off()
##
grade.sd <- by(xbcf.cont[[2]]$pred,w[test,]$hsgrade_pct,function(x) sd(colMeans(x)))
grade.mean <- by(xbcf.cont[[2]]$pred,w[test,]$hsgrade_pct,function(x) mean(colMeans(x)))
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
