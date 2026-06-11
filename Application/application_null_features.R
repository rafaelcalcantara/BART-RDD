# Setup--------------------------------------------------------------
if ("Figures" %in% list.files() == FALSE) dir.create("Figures")
library(stochtree)
library(rpart)
library(rpart.plot)
library(xtable)
library(MASS)
library(doParallel)
library(foreach)
run.barddt <- FALSE ## toggle to fit BARDDT
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
## Fitting the models----------------------------------------------------------
w$totcredits_year1 <- factor(w$totcredits_year1,ordered=TRUE)
w$male <- factor(w$male,ordered=FALSE)
w$bpl_north_america <- factor(w$bpl_north_america,ordered=FALSE)
w$loc_campus1 <- factor(w$loc_campus1,ordered=FALSE)
w$loc_campus2 <- factor(w$loc_campus2,ordered=FALSE)
w$loc_campus3 <- factor(w$loc_campus3,ordered=FALSE)
### Add null features
for (i in 1:20)
{
  temp <- rnorm(n)
  w[,ncol(w)+1] <- temp
  colnames(w)[ncol(w)] <- paste0("null.",i)
}
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
  saveRDS(pred,"bart_rdd_posterior_null_features.rds")
} else
{
  pred <- readRDS("bart_rdd_posterior_null_features.rds")
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
pdf("Figures/cate_gpa_null_features.pdf")
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

contour.min <- 1

drawlabels <- FALSE
lwd <- 1.2
lvl <- c(1,15,50)

pdf("Figures/contour_null_features.pdf",width=9,height=4.5)
contour(kdens.barddt,bty='n',xlab="Group A",ylab="Group B",col="black",levels=lvl,drawlabels=drawlabels,lwd=lwd+0.25)
abline(a=0,b=1)
dev.off()
## Correlations
vals <- cor(rowMeans(pred),w[test,9:28])
n   <- length(vals)
row_names <- paste0("$W_{", seq_len(n), "}$")
rounded   <- formatC(round(vals, 3), format = "f", digits = 3)

sink("Tables/null_correlations.tex")
cat("\\begin{table}[ht]\n")
cat("\\centering\n")
cat("\\begin{tabular}{lc}\n")
cat("\\hline\n")
cat("& Value \\\\\n")
cat("\\hline\n")
for (i in seq_len(n)) {
  cat(row_names[i], "&", rounded[i], "\\\\\n")
}
cat("\\hline\n")
cat("\\end{tabular}\n")
cat("\\caption{}\n")
cat("\\label{tab:}\n")
cat("\\end{table}\n")
sink()
