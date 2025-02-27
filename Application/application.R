# Setup--------------------------------------------------------------
library(stochtree)
library(rpart)
library(rpart.plot)
library(xtable)
library(MASS)
run.model <- FALSE ## toggle to actually fit BARDDT
## Read data
data <- read.csv("gpa.csv")
y <- data$nextGPA
x <- data$X
x <- x/sd(x)
w <- data[,4:11]
c <- 0
n <- nrow(data)
z <- as.numeric(x>c)
h <- 0.1
id0 <- -h < x & x < 0
id1 <- 0 < x & x< h
if ("Figures" %in% list.files() == FALSE) dir.create("Figures")
if ("Tables" %in% list.files() == FALSE) dir.create("Tables")
if ("Results" %in% list.files() == FALSE) dir.create("Results")
## Summary statistics-----------------------------------------------
### Tables
#### Full sample
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
sum.tab.tex <- stargazer::stargazer(sum.tab,summary=F,rownames = F,
                                    label = "tab:sum.stat", title = "Summary statistics",
                                    table.placement = "!htpb", font.size = "small",
                                    notes = c("\\tiny Sample size:",paste(c("\\tiny Total:","Campus 1:","Campus 2:","Campus 3:"),obs,collapse="; ")))
writeLines(sum.tab.tex,"Tables/application_sum_stat.tex")
#### Test sample
sum.stat.h <- function(dat,id0,id1) 
{
  out <- round(cbind(t(apply(dat[id0,],2,function(i) c(Mean=mean(i),SD=sd(i)))),
                     t(apply(dat[id1,],2,function(i) c(Mean=mean(i),SD=sd(i))))), 3)
  return(out)
}
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
sum.tab.test.tex <- stargazer::stargazer(sum.tab,summary=F,rownames = F,
                                         label = "tab:sum.stat.h", title = "Summary statistics - identification strip",
                                         table.placement = "!htpb", font.size = "small",
                                         notes = c("\\tiny Sample size (control/treatment):",
                                                   paste(c("\\tiny Total:","Campus 1:","Campus 2:","Campus 3:"),
                                                         paste(obs0,obs1,sep="/"),collapse="; ")))
writeLines(sum.tab.test.tex,"Tables/application_sum_stat_test.tex")
### Plots
#### W
l0 <- lapply(1:ncol(w), function(i) predict(loess(w~x,data=data.frame(w=w[id0,i],x=x[id0])),se=T)[1:2])
l0 <- lapply(l0,function(i) cbind(i$fit,i$fit-1.96*i$se.fit,i$fit+1.96*i$se.fit))
l1 <- lapply(1:ncol(w), function(i) predict(loess(w~x,data=data.frame(w=w[id1,i],x=x[id1])),se=T)[1:2])
l1 <- lapply(l1,function(i) cbind(i$fit,i$fit-1.96*i$se.fit,i$fit+1.96*i$se.fit))
#### Y
yhat0 <- predict(loess(y~x,data=data.frame(y=y[id0],x=x[id0])),se=T)[1:2]
yhat0 <- cbind(yhat0$fit,yhat0$fit-1.96*yhat0$se.fit,yhat0$fit+1.96*yhat0$se.fit)
yhat1 <- predict(loess(y~x,data=data.frame(y=y[id1],x=x[id1])),se=T)[1:2]
yhat1 <- cbind(yhat1$fit,yhat1$fit-1.96*yhat1$se.fit,yhat1$fit+1.96*yhat1$se.fit)
####
xplot <- x[id0 | id1]
pdf("Figures/application_summary.pdf",width=9,height=12)
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
## Fitting the model-----------------------------
test <- id0 | id1
w$male <- factor(w$male,ordered=FALSE)
w$bpl_north_america <- factor(w$bpl_north_america,ordered=FALSE)
w$loc_campus1 <- factor(w$loc_campus1,ordered=FALSE)
w$loc_campus2 <- factor(w$loc_campus2,ordered=FALSE)
w$loc_campus3 <- factor(w$loc_campus3,ordered=FALSE)
if (isTRUE(run.model))
{
  barddt.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
  barddt.mean.parmlist <- list(num_trees=150, min_samples_leaf=20, alpha=0.95, beta=2,
                               max_depth=20, sample_sigma2_leaf=FALSE, sigma2_leaf_init = diag(rep(0.1/150,4)))
  B <- cbind(z*x,(1-z)*x, z,rep(1,n))
  B1 <- cbind(rep(c,n), rep(0,n), rep(1,n), rep(1,n))
  B0 <- cbind(rep(0,n), rep(c,n), rep(0,n), rep(1,n))
  time.fit <- system.time({
    barddt.fit = stochtree::bart(X_train= cbind(x,w), y_train=y,
                                 W_train = B, mean_forest_params=barddt.mean.parmlist,
                                 general_params=barddt.global.parmlist,
                                 num_mcmc=10000,num_gfr=30)
  })
  time.pred <- system.time({
    B1 <- B1[test,]
    B0 <- B0[test,]
    xmat_test <- cbind(x=rep(0,n),w)[test,]
    pred1 <- predict(barddt.fit,xmat_test,B1)$y_hat
    pred0 <- predict(barddt.fit,xmat_test,B0)$y_hat
    post <- pred1-pred0
  })
  saveRDS(post,"Results/bart_rdd_posterior.rds")
}
pred <- readRDS("Results/bart_rdd_posterior.rds")
###
cate <- rpart(y~.,data.frame(y=rowMeans(pred),w[test,]),control = rpart.control(cp=0.015))
## Define separate colors for left and rightmost nodes
plot.cart <- function(rpart.obj)
{
  rpart.frame <- rpart.obj$frame
  left <- which.min(rpart.frame$yval)
  right <- which.max(rpart.frame$yval)
  nodes <- rep(NA,nrow(rpart.frame))
  for (i in 1:length(nodes))
  {
    if (rpart.frame$yval[i]==rpart.frame$yval[right]) nodes[i] <- "gold2"
    else if (rpart.frame$yval[i]==rpart.frame$yval[left]) nodes[i] <- "tomato3"
    else nodes[i] <- "lightblue3"
  }
  return(nodes)
}
## Plot CART tree
pdf("Figures/cate_gpa.pdf")
par(mfrow=c(1,1))
rpart.plot(cate,main="",box.col=plot.cart(cate))
dev.off()
## Define function to produce KD estimates of the joint distribution of two subgroups
cate.kde <- function(rpart.obj,pred)
{
  rpart.frame <- rpart.obj$frame
  left <- rpart.obj$where==which.min(rpart.frame$yval)
  right <- rpart.obj$where==which.max(rpart.frame$yval)
  ## Calculate CATE posterior for groups A and B
  cate.a <- do.call("cbind",by(pred,left, colMeans))
  cate.b <- do.call("cbind",by(pred,right, colMeans))
  cate.a <- cate.a[,2]
  cate.b <- cate.b[,2]
  ## Estimate kernel density
  denshat <- MASS::kde2d(cate.a, cate.b, n=200)
  return(denshat)
}
pdf("Figures/cate_difference.pdf")
contour(cate.kde(cate,pred),bty='n',xlab="Group A",ylab="Group B")
abline(a=0,b=1)
dev.off()