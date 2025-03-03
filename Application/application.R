# Setup--------------------------------------------------------------
library(stochtree)
library(rpart)
library(rpart.plot)
library(xtable)
library(MASS)
library(doParallel)
library(foreach)
run.model <- TRUE ## toggle to actually fit BARDDT
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
id0 <- -h < x
id1 <- x < h
test <- -h < x & x < h
ntest <- sum(test)
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
## Fitting the model-----------------------------
w$totcredits_year1 <- factor(w$totcredits_year1,ordered=TRUE)
w$male <- factor(w$male,ordered=FALSE)
w$bpl_north_america <- factor(w$bpl_north_america,ordered=FALSE)
w$loc_campus1 <- factor(w$loc_campus1,ordered=FALSE)
w$loc_campus2 <- factor(w$loc_campus2,ordered=FALSE)
w$loc_campus3 <- factor(w$loc_campus3,ordered=FALSE)
if (isTRUE(run.model))
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
  pred <- readRDS("Results/bart_rdd_posterior.rds")
}
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