# Setup--------------------------------------------------------------
if ("Figures" %in% list.files() == FALSE) dir.create("Figures")
library(stochtree)
library(rpart)
library(rpart.plot)
library(xtable)
library(MASS)
library(doParallel)
library(foreach)
library(bayeslm)
run.barddt <- FALSE ## toggle to fit BARDDT
run.sbart <- FALSE ## toggle to fit S-BART
run.tbart <- FALSE ## toggle to fit T-BART
run.llr <- FALSE ## toggle to fit LLR
## Read data-------------------------------------------------------------------
data <- read.csv("gpa.csv")
y <- data$nextGPA
x <- data$X
x <- x/sd(x)
w <- data[,4:11]
c <- 0
n <- nrow(data)
z <- as.numeric(x>c)
Owidth <- 0.1
id0 <- -Owidth < x
id1 <- x < Owidth
test <- -Owidth < x & x < Owidth
ntest <- sum(test)
if ("Figures" %in% list.files() == FALSE) dir.create("Figures")
if ("Tables" %in% list.files() == FALSE) dir.create("Tables")
if ("Results" %in% list.files() == FALSE) dir.create("Results")
## Summary statistics----------------------------------------------------------
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
## Fitting the models----------------------------------------------------------
w$totcredits_year1 <- factor(w$totcredits_year1,ordered=TRUE)
w$male <- factor(w$male,ordered=FALSE)
w$bpl_north_america <- factor(w$bpl_north_america,ordered=FALSE)
w$loc_campus1 <- factor(w$loc_campus1,ordered=FALSE)
w$loc_campus2 <- factor(w$loc_campus2,ordered=FALSE)
w$loc_campus3 <- factor(w$loc_campus3,ordered=FALSE)
## BARDDT
if (isTRUE(run.barddt))
{
  ## We will sample multiple chains sequentially
  num_chains <- 20
  num_gfr <- 2
  num_burnin <- 0
  num_mcmc <- 500
  bart_models <- list()
  ## Define basis functions for training and testing
  B <- cbind(z*x,(1-z)*x, z,rep(1,n))
  B1 <- cbind(rep(c,n), rep(0,n), rep(1,n), rep(1,n))
  B0 <- cbind(rep(0,n), rep(c,n), rep(0,n), rep(1,n))
  B1 <- B1[test,]
  B0 <- B0[test,]
  B_test <- rbind(B1,B0)
  xmat_test <- cbind(x=rep(0,n),w)[test,]
  xmat_test <- rbind(xmat_test,xmat_test)
  ### We combine the basis for Z=1 and Z=0 to feed it to the BART call and get the Y(z) predictions instantaneously
  ### Then we separate the posterior matrix between each Z and calculate the CATE prediction
  ## Sampling trees in parallel
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  start_time <- Sys.time()
  bart_model_outputs <- foreach (i = 1:num_chains) %dopar% {
    random_seed <- i
    ## Lists to define BARDDT parameters
    barddt.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
    barddt.mean.parmlist <- list(num_trees=50, min_samples_leaf=20, alpha=0.95, beta=2,
                                 max_depth=20, sample_sigma2_leaf=FALSE, sigma2_leaf_init = diag(rep(0.1/50,4)))
    bart_model <- stochtree::bart(
      X_train = cbind(x,w), leaf_basis_train = B, y_train = y,
      X_test = xmat_test, leaf_basis_test = B_test,
      num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc,
      general_params = barddt.global.parmlist, mean_forest_params = barddt.mean.parmlist
    )
    bart_model <- bart_model$y_hat_test[1:ntest,]-bart_model$y_hat_test[(ntest+1):(2*ntest),]
  }
  stopCluster(cl)
  ## Combine CATE predictions
  pred <- do.call("cbind",bart_model_outputs)

  end_time <- Sys.time()

  print(end_time - start_time)
  ## Save the results
  saveRDS(pred,"bart_rdd_posterior.rds")
} else
{
  pred <- readRDS("bart_rdd_posterior.rds")
}
## S-BART
if (isTRUE(run.sbart))
{
  ## We will sample multiple chains sequentially
  num_chains <- 20
  num_gfr <- 2
  num_burnin <- 0
  num_mcmc <- 500
  bart_models <- list()
  ## Define test set
  xmat_test.1 <- cbind(x=c,z=1,w)[test,]
  xmat_test.0 <- cbind(x=c,z=0,w)[test,]
  ## Sampling trees in parallel
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  start_time <- Sys.time()
  bart_model_outputs <- foreach (i = 1:num_chains) %dopar% {
    random_seed <- i
    ## Lists to define BARDDT parameters
    barddt.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
    barddt.mean.parmlist <- list(num_trees=50, min_samples_leaf=20, alpha=0.95, beta=2,
                                 max_depth=20, sample_sigma2_leaf=TRUE)
    bart_model <- stochtree::bart(
      X_train = cbind(x,z,w), y_train = y,
      num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc,
      general_params = barddt.global.parmlist, mean_forest_params = barddt.mean.parmlist
    )
    pred1 <- predict(bart_model,xmat_test.1)$y_hat
    pred0 <- predict(bart_model,xmat_test.0)$y_hat
    post <- pred1-pred0
  }
  stopCluster(cl)
  ## Combine CATE predictions
  pred.sbart <- do.call("cbind",bart_model_outputs)

  end_time <- Sys.time()

  print(end_time - start_time)
  ## Save the results
  saveRDS(pred.sbart,"sbart_posterior.rds")
} else
{
  pred.sbart <- readRDS("sbart_posterior.rds")
}
## T-BART
if (isTRUE(run.tbart))
{
  ## We will sample multiple chains sequentially
  num_chains <- 20
  num_gfr <- 2
  num_burnin <- 0
  num_mcmc <- 500
  bart_models <- list()
  ## Define test set
  xmat_test <- cbind(x=c,w)[test,]
  ## Sampling trees in parallel
  ncores <- 5
  cl <- makeCluster(ncores)
  registerDoParallel(cl)

  start_time <- Sys.time()
  bart_model_outputs <- foreach (i = 1:num_chains) %dopar% {
    random_seed <- i
    ## Lists to define BARDDT parameters
    barddt.global.parmlist <- list(standardize=T,sample_sigma_global=TRUE,sigma2_global_init=0.1)
    barddt.mean.parmlist <- list(num_trees=50, min_samples_leaf=20, alpha=0.95, beta=2,
                                 max_depth=20, sample_sigma2_leaf=TRUE)
    ## X>c regression
    bart_model <- stochtree::bart(
      X_train = cbind(x,w)[x>c,], y_train = y[x>c],
      num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc,
      general_params = barddt.global.parmlist, mean_forest_params = barddt.mean.parmlist
    )
    pred1 <- predict(bart_model,xmat_test)$y_hat
    ## X<=c regression
    bart_model <- stochtree::bart(
      X_train = cbind(x,w)[x<=c,], y_train = y[x<=c],
      num_gfr = num_gfr, num_burnin = num_burnin, num_mcmc = num_mcmc,
      general_params = barddt.global.parmlist, mean_forest_params = barddt.mean.parmlist
    )
    pred0 <- predict(bart_model,xmat_test)$y_hat
    ## Posterior
    post <- pred1-pred0
  }
  stopCluster(cl)
  ## Combine CATE predictions
  pred.tbart <- do.call("cbind",bart_model_outputs)

  end_time <- Sys.time()

  print(end_time - start_time)
  ## Save the results
  saveRDS(pred.tbart,"tbart_posterior.rds")
} else
{
  pred.tbart <- readRDS("tbart_posterior.rds")
}
## LLR
if (isTRUE(run.llr))
{
  dfw <- data.frame(w=w)
  fmla <- as.formula(paste0("y~(",names(dfw)[1],"+",paste(paste("as.factor(",names(dfw)[-1],")",sep=""),collapse="+"),")*x*z"))
  df <- data.frame(x=x,w=w,y=y,z=z)
  df$z <- as.factor(df$z)
  h <- rdrobust::rdbwselect(df$y,df$x,c,covs=df$w)$bws[4]
  df.train <- subset(df,c-h<=x & x<=c+h)
  #######
  tt = terms(as.formula(fmla))
  TT = delete.response(tt)
  mf <- model.frame(TT, data = df.train)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  #######
  poly.fit <- bayeslm(fmla,df.train,penalize=rep(1,ncol(X)-1),singular=TRUE,prior="horseshoe")
  ######
  df.test <- df[test,]
  xtest.a <- df.test
  xtest.a$x <- c
  xtest.a$z <- "1"
  xtest.b <- df.test
  xtest.b$x <- c
  xtest.b$z <- "0"
  xtest <- rbind(xtest.a,xtest.b,df.test)
  ####
  tt = terms(as.formula(fmla))
  TT = delete.response(tt)
  mf <- model.frame(TT, data = xtest)
  X <- model.matrix(attr(mf, "terms"), data = mf)
  X1 <- X[1:(nrow(X)/3),]
  X0 <- X[(nrow(X)/3+1):(2*nrow(X)/3),]
  pred.llr <- apply(poly.fit$beta,1, function(i) (X1-X0) %*% i)
  ## Save results
  saveRDS(pred.llr,"llr.rds")
} else
{
  pred.llr <- readRDS("llr.rds")
}
# Figure 7---------------------------------------------------------------------
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
# Figure 8---------------------------------------------------------------------
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
## Produce KD estimates
kdens.barddt <- cate.kde(cate,pred)
kdens.sbart <- cate.kde(cate,pred.sbart)
kdens.tbart <- cate.kde(cate,pred.tbart)
kdens.llr <- cate.kde(cate,pred.llr)

colors <- c("magenta3","dodgerblue","green3")

contour.min <- 1

ylim <- quantile(c(kdens.barddt$y,kdens.tbart$y,kdens.sbart$y),c(0.3,0.9))
xlim <- quantile(c(kdens.barddt$x,kdens.tbart$x,kdens.sbart$x),c(0.3,0.9))
drawlabels <- FALSE
lwd <- 1.2
lvl <- c(1,15,50)

pdf("Figures/contour.pdf",width=9,height=4.5)
layout(matrix(c(1,2,1,3,1,4),ncol=3),heights = c(0.25,1.75))
par(mar=c(0.1,4,0.05,0.1))
plot.new()
legend("center",col=c("black",colors),legend = c("BARDDT","T-BART","S-BART","Horseshoe"),ncol=4,lty=1,lwd=1.5)
par(mar=c(4,4,0.1,0.1))
contour(kdens.barddt,bty='n',xlab="Group A",ylab="Group B",col="black",levels=lvl,ylim=ylim,xlim=xlim,drawlabels=drawlabels,lwd=lwd+0.25)
contour(cate.kde(cate,pred.tbart),bty='n',xlab="Group A",ylab="Group B",add=TRUE,col=colors[1],levels=lvl,drawlabels=drawlabels,lwd=lwd)
abline(a=0,b=1)
contour(kdens.barddt,bty='n',xlab="Group A",ylab="Group B",col="black",levels=lvl,ylim=ylim,xlim=xlim,drawlabels=drawlabels,lwd=lwd+0.25)
contour(cate.kde(cate,pred.sbart),bty='n',xlab="Group A",ylab="Group B",add=TRUE,col=colors[2],levels=lvl,drawlabels=drawlabels,lwd=lwd)
abline(a=0,b=1)
contour(kdens.barddt,bty='n',xlab="Group A",ylab="Group B",col="black",levels=lvl,ylim=ylim,xlim=xlim,drawlabels=drawlabels,lwd=lwd+0.25)
contour(cate.kde(cate,pred.llr),bty='n',xlab="Group A",ylab="Group B",add=TRUE,col=colors[3],levels=lvl,drawlabels=drawlabels,lwd=lwd)
abline(a=0,b=1)
dev.off()
# ### BARDDT
# pdf("Figures/cate_difference.pdf")
# contour(cate.kde(cate,pred),bty='n',xlab="Group A",ylab="Group B")
# abline(a=0,b=1)
# dev.off()
# ### S-BART
# pdf("Figures/cate_difference_sbart.pdf")
# contour(cate.kde(cate,pred.sbart),bty='n',xlab="Group A",ylab="Group B")
# abline(a=0,b=1)
# dev.off()
# ### T-BART
# pdf("Figures/cate_difference_tbart.pdf")
# contour(cate.kde(cate,pred.tbart),bty='n',xlab="Group A",ylab="Group B")
# abline(a=0,b=1)
# dev.off()
# ### LLR
# pdf("Figures/cate_difference_llr.pdf")
# contour(cate.kde(cate,pred.llr),bty='n',xlab="Group A",ylab="Group B")
# abline(a=0,b=1)
# dev.off()
