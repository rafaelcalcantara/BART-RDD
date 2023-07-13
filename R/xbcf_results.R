library(XBART)
library(parallel)
library(foreach)
library(doParallel)
setwd("~/Documents/Git/XBCF-RDD/")
s <- 1000
num_sweeps <- 120
burnin     <- 20
sample     <- (burnin+1):num_sweeps
### Function to read results files
readFiles <- function(s,dgp,file,p,ab)
{
    foreach(x=1:s,.multicombine=T) %do%
        readRDS(paste0("Results/",file,"_",dgp,ab,"_",p,"_",x,".rds"))
}
readFiles2 <- function(s,dgp,file)
{
    foreach(x=1:s,.multicombine=T) %do%
        readRDS(paste0("Results/",file,"_",dgp,"_",x,".rds"))
}
#### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## P=4
### DGP1a
results <- readFiles(s,"1","xbcf",4,"a")
#### Obtain ATE posterior
ate.sum.1a.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP1b
results <- readFiles(s,"1","xbcf",4,"b")
#### Obtain ATE posterior
ate.sum.1b.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP2a
results <- readFiles(s,"2","xbcf",4,"a")
#### Obtain ATE posterior
ate.sum.2a.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP2b
results <- readFiles(s,"2","xbcf",4,"b")
#### Obtain ATE posterior
ate.sum.2b.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP3a
results <- readFiles(s,"3","xbcf",4,"a")
#### Obtain ATE posterior
ate.sum.3a.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP3b
results <- readFiles(s,"3","xbcf",4,"b")
#### Obtain ATE posterior
ate.sum.3b.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP4a
results <- readFiles(s,"4","xbcf",4,"a")
#### Obtain ATE posterior
ate.sum.4a.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP4b
results <- readFiles(s,"4","xbcf",4,"b")
#### Obtain ATE posterior
ate.sum.4b.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP5a
results <- readFiles(s,"5","xbcf",4,"a")
#### Obtain ATE posterior
ate.sum.5a.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP5b
results <- readFiles(s,"5","xbcf",4,"b")
#### Obtain ATE posterior
ate.sum.5b.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP6a
results <- readFiles(s,"6","xbcf",4,"a")
#### Obtain ATE posterior
ate.sum.6a.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP6b
results <- readFiles(s,"6","xbcf",4,"b")
#### Obtain ATE posterior
ate.sum.6b.4 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### Read data and save true ATE
data <- readRDS("Data/DGP1_4.rds")
data <- data[1:s]
ate1 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP2_4.rds")
data <- data[1:s]
ate2 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP3_4.rds")
data <- data[1:s]
ate3 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP4_4.rds")
data <- data[1:s]
ate4 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP5_4.rds")
data <- data[1:s]
ate5 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP6_4.rds")
data <- data[1:s]
ate6 <- sapply(data,function(x) x$ate)
### Store results
res.mat.1.4 <- matrix(NA,7,7)
res.mat.2.4 <- matrix(NA,7,7)
res.mat.1.4[1,] <- res.mat.2.4[1,] <- c("","","Controls","MSE","Variance","Coverage","Size")
res.mat.1.4[2,] <- res.mat.2.4[2,] <- c("XBCF","DGP 1","","","","","")
res.mat.1.4[4,] <- res.mat.2.4[4,] <- c("XBCF","DGP 2","","","","","")
res.mat.1.4[6,] <- res.mat.2.4[6,] <- c("XBCF","DGP 3","","","","","")
res.mat.1.4[2:7,3] <- res.mat.2.4[2:7,3] <- c("No","Yes")
#### MSE
##### Table 1
res.mat.1.4[2,4] <- round(mean((ate.sum.1a.4[,1]-ate1)^2),3)
res.mat.1.4[3,4] <- round(mean((ate.sum.1b.4[,1]-ate1)^2),3)
res.mat.1.4[4,4] <- round(mean((ate.sum.2a.4[,1]-ate2)^2),3)
res.mat.1.4[5,4] <- round(mean((ate.sum.2b.4[,1]-ate2)^2),3)
res.mat.1.4[6,4] <- round(mean((ate.sum.3a.4[,1]-ate3)^2),3)
res.mat.1.4[7,4] <- round(mean((ate.sum.3b.4[,1]-ate3)^2),3)
##### Table 2
res.mat.2.4[2,4] <- round(mean((ate.sum.4a.4[,1]-ate4)^2),3)
res.mat.2.4[3,4] <- round(mean((ate.sum.4b.4[,1]-ate4)^2),3)
res.mat.2.4[4,4] <- round(mean((ate.sum.5a.4[,1]-ate5)^2),3)
res.mat.2.4[5,4] <- round(mean((ate.sum.5b.4[,1]-ate5)^2),3)
res.mat.2.4[6,4] <- round(mean((ate.sum.6a.4[,1]-ate6)^2),3)
res.mat.2.4[7,4] <- round(mean((ate.sum.6b.4[,1]-ate6)^2),3)
#### Variance
##### Table 1
res.mat.1.4[2,5] <- round(var(ate.sum.1a.4[,1]),3)
res.mat.1.4[3,5] <- round(var(ate.sum.1b.4[,1]),3)
res.mat.1.4[4,5] <- round(var(ate.sum.2a.4[,1]),3)
res.mat.1.4[5,5] <- round(var(ate.sum.2b.4[,1]),3)
res.mat.1.4[6,5] <- round(var(ate.sum.3a.4[,1]),3)
res.mat.1.4[7,5] <- round(var(ate.sum.3b.4[,1]),3)
##### Table 2
res.mat.2.4[2,5] <- round(var(ate.sum.4a.4[,1]),3)
res.mat.2.4[3,5] <- round(var(ate.sum.4b.4[,1]),3)
res.mat.2.4[4,5] <- round(var(ate.sum.5a.4[,1]),3)
res.mat.2.4[5,5] <- round(var(ate.sum.5b.4[,1]),3)
res.mat.2.4[6,5] <- round(var(ate.sum.6a.4[,1]),3)
res.mat.2.4[7,5] <- round(var(ate.sum.6b.4[,1]),3)
#### Coverage
##### Table 1
res.mat.1.4[2,6] <- 100*mean(ate.sum.1a.4[,2]<=ate1 & ate1<=ate.sum.1a.4[,3])
res.mat.1.4[3,6] <- 100*mean(ate.sum.1b.4[,2]<=ate1 & ate1<=ate.sum.1b.4[,3])
res.mat.1.4[4,6] <- 100*mean(ate.sum.2a.4[,2]<=ate2 & ate2<=ate.sum.2a.4[,3])
res.mat.1.4[5,6] <- 100*mean(ate.sum.2b.4[,2]<=ate2 & ate2<=ate.sum.2b.4[,3])
res.mat.1.4[6,6] <- 100*mean(ate.sum.3a.4[,2]<=ate3 & ate3<=ate.sum.3a.4[,3])
res.mat.1.4[7,6] <- 100*mean(ate.sum.3b.4[,2]<=ate3 & ate3<=ate.sum.3b.4[,3])
##### Table 2
res.mat.2.4[2,6] <- 100*mean(ate.sum.4a.4[,2]<=ate4 & ate4<=ate.sum.4a.4[,3])
res.mat.2.4[3,6] <- 100*mean(ate.sum.4b.4[,2]<=ate4 & ate4<=ate.sum.4b.4[,3])
res.mat.2.4[4,6] <- 100*mean(ate.sum.5a.4[,2]<=ate5 & ate5<=ate.sum.5a.4[,3])
res.mat.2.4[5,6] <- 100*mean(ate.sum.5b.4[,2]<=ate5 & ate5<=ate.sum.5b.4[,3])
res.mat.2.4[6,6] <- 100*mean(ate.sum.6a.4[,2]<=ate6 & ate6<=ate.sum.6a.4[,3])
res.mat.2.4[7,6] <- 100*mean(ate.sum.6b.4[,2]<=ate6 & ate6<=ate.sum.6b.4[,3])
#### Size
##### Table 1
res.mat.1.4[2,7] <- round(mean(-ate.sum.1a.4[,2]+ate.sum.1a.4[,3]),2)
res.mat.1.4[3,7] <- round(mean(-ate.sum.1b.4[,2]+ate.sum.1b.4[,3]),2)
res.mat.1.4[4,7] <- round(mean(-ate.sum.2a.4[,2]+ate.sum.2a.4[,3]),2)
res.mat.1.4[5,7] <- round(mean(-ate.sum.2b.4[,2]+ate.sum.2b.4[,3]),2)
res.mat.1.4[6,7] <- round(mean(-ate.sum.3a.4[,2]+ate.sum.3a.4[,3]),2)
res.mat.1.4[7,7] <- round(mean(-ate.sum.3b.4[,2]+ate.sum.3b.4[,3]),2)
##### Table 2
res.mat.2.4[2,7] <- round(mean(-ate.sum.4a.4[,2]+ate.sum.4a.4[,3]),2)
res.mat.2.4[3,7] <- round(mean(-ate.sum.4b.4[,2]+ate.sum.4b.4[,3]),2)
res.mat.2.4[4,7] <- round(mean(-ate.sum.5a.4[,2]+ate.sum.5a.4[,3]),2)
res.mat.2.4[5,7] <- round(mean(-ate.sum.5b.4[,2]+ate.sum.5b.4[,3]),2)
res.mat.2.4[6,7] <- round(mean(-ate.sum.6a.4[,2]+ate.sum.6a.4[,3]),2)
res.mat.2.4[7,7] <- round(mean(-ate.sum.6b.4[,2]+ate.sum.6b.4[,3]),2)
## P=6
### DGP1a
results <- readFiles(s,"1","xbcf",6,"a")
#### Obtain ATE posterior
ate.sum.1a.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP1b
results <- readFiles(s,"1","xbcf",6,"b")
#### Obtain ATE posterior
ate.sum.1b.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP2a
results <- readFiles(s,"2","xbcf",6,"a")
#### Obtain ATE posterior
ate.sum.2a.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP2b
results <- readFiles(s,"2","xbcf",6,"b")
#### Obtain ATE posterior
ate.sum.2b.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP3a
results <- readFiles(s,"3","xbcf",6,"a")
#### Obtain ATE posterior
ate.sum.3a.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP3b
results <- readFiles(s,"3","xbcf",6,"b")
#### Obtain ATE posterior
ate.sum.3b.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP4a
results <- readFiles(s,"4","xbcf",6,"a")
#### Obtain ATE posterior
ate.sum.4a.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP4b
results <- readFiles(s,"4","xbcf",6,"b")
#### Obtain ATE posterior
ate.sum.4b.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP5a
results <- readFiles(s,"5","xbcf",6,"a")
#### Obtain ATE posterior
ate.sum.5a.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP5b
results <- readFiles(s,"5","xbcf",6,"b")
#### Obtain ATE posterior
ate.sum.5b.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP6a
results <- readFiles(s,"6","xbcf",6,"a")
#### Obtain ATE posterior
ate.sum.6a.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP6b
results <- readFiles(s,"6","xbcf",6,"b")
#### Obtain ATE posterior
ate.sum.6b.6 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### Read data and save true ATE
data <- readRDS("Data/DGP1_6.rds")
data <- data[1:s]
ate1 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP2_6.rds")
data <- data[1:s]
ate2 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP3_6.rds")
data <- data[1:s]
ate3 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP4_6.rds")
data <- data[1:s]
ate4 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP5_6.rds")
data <- data[1:s]
ate5 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP6_6.rds")
data <- data[1:s]
ate6 <- sapply(data,function(x) x$ate)
### Store results
res.mat.1.6 <- matrix(NA,7,7)
res.mat.2.6 <- matrix(NA,7,7)
res.mat.1.6[1,] <- res.mat.2.6[1,] <- c("","","Controls","MSE","Variance","Coverage","Size")
res.mat.1.6[2,] <- res.mat.2.6[2,] <- c("XBCF","DGP 1","","","","","")
res.mat.1.6[4,] <- res.mat.2.6[4,] <- c("XBCF","DGP 2","","","","","")
res.mat.1.6[6,] <- res.mat.2.6[6,] <- c("XBCF","DGP 3","","","","","")
res.mat.1.6[2:7,3] <- res.mat.2.6[2:7,3] <- c("No","Yes")
#### MSE
##### Table 1
res.mat.1.6[2,4] <- round(mean((ate.sum.1a.6[,1]-ate1)^2),3)
res.mat.1.6[3,4] <- round(mean((ate.sum.1b.6[,1]-ate1)^2),3)
res.mat.1.6[4,4] <- round(mean((ate.sum.2a.6[,1]-ate2)^2),3)
res.mat.1.6[5,4] <- round(mean((ate.sum.2b.6[,1]-ate2)^2),3)
res.mat.1.6[6,4] <- round(mean((ate.sum.3a.6[,1]-ate3)^2),3)
res.mat.1.6[7,4] <- round(mean((ate.sum.3b.6[,1]-ate3)^2),3)
##### Table 2
res.mat.2.6[2,4] <- round(mean((ate.sum.4a.6[,1]-ate4)^2),3)
res.mat.2.6[3,4] <- round(mean((ate.sum.4b.6[,1]-ate4)^2),3)
res.mat.2.6[4,4] <- round(mean((ate.sum.5a.6[,1]-ate5)^2),3)
res.mat.2.6[5,4] <- round(mean((ate.sum.5b.6[,1]-ate5)^2),3)
res.mat.2.6[6,4] <- round(mean((ate.sum.6a.6[,1]-ate6)^2),3)
res.mat.2.6[7,4] <- round(mean((ate.sum.6b.6[,1]-ate6)^2),3)
#### Variance
##### Table 1
res.mat.1.6[2,5] <- round(var(ate.sum.1a.6[,1]),3)
res.mat.1.6[3,5] <- round(var(ate.sum.1b.6[,1]),3)
res.mat.1.6[4,5] <- round(var(ate.sum.2a.6[,1]),3)
res.mat.1.6[5,5] <- round(var(ate.sum.2b.6[,1]),3)
res.mat.1.6[6,5] <- round(var(ate.sum.3a.6[,1]),3)
res.mat.1.6[7,5] <- round(var(ate.sum.3b.6[,1]),3)
##### Table 2
res.mat.2.6[2,5] <- round(var(ate.sum.4a.6[,1]),3)
res.mat.2.6[3,5] <- round(var(ate.sum.4b.6[,1]),3)
res.mat.2.6[4,5] <- round(var(ate.sum.5a.6[,1]),3)
res.mat.2.6[5,5] <- round(var(ate.sum.5b.6[,1]),3)
res.mat.2.6[6,5] <- round(var(ate.sum.6a.6[,1]),3)
res.mat.2.6[7,5] <- round(var(ate.sum.6b.6[,1]),3)
#### Coverage
##### Table 1
res.mat.1.6[2,6] <- 100*mean(ate.sum.1a.6[,2]<=ate1 & ate1<=ate.sum.1a.6[,3])
res.mat.1.6[3,6] <- 100*mean(ate.sum.1b.6[,2]<=ate1 & ate1<=ate.sum.1b.6[,3])
res.mat.1.6[4,6] <- 100*mean(ate.sum.2a.6[,2]<=ate2 & ate2<=ate.sum.2a.6[,3])
res.mat.1.6[5,6] <- 100*mean(ate.sum.2b.6[,2]<=ate2 & ate2<=ate.sum.2b.6[,3])
res.mat.1.6[6,6] <- 100*mean(ate.sum.3a.6[,2]<=ate3 & ate3<=ate.sum.3a.6[,3])
res.mat.1.6[7,6] <- 100*mean(ate.sum.3b.6[,2]<=ate3 & ate3<=ate.sum.3b.6[,3])
##### Table 2
res.mat.2.6[2,6] <- 100*mean(ate.sum.4a.6[,2]<=ate4 & ate4<=ate.sum.4a.6[,3])
res.mat.2.6[3,6] <- 100*mean(ate.sum.4b.6[,2]<=ate4 & ate4<=ate.sum.4b.6[,3])
res.mat.2.6[4,6] <- 100*mean(ate.sum.5a.6[,2]<=ate5 & ate5<=ate.sum.5a.6[,3])
res.mat.2.6[5,6] <- 100*mean(ate.sum.5b.6[,2]<=ate5 & ate5<=ate.sum.5b.6[,3])
res.mat.2.6[6,6] <- 100*mean(ate.sum.6a.6[,2]<=ate6 & ate6<=ate.sum.6a.6[,3])
res.mat.2.6[7,6] <- 100*mean(ate.sum.6b.6[,2]<=ate6 & ate6<=ate.sum.6b.6[,3])
#### Size
##### Table 1
res.mat.1.6[2,7] <- round(mean(-ate.sum.1a.6[,2]+ate.sum.1a.6[,3]),2)
res.mat.1.6[3,7] <- round(mean(-ate.sum.1b.6[,2]+ate.sum.1b.6[,3]),2)
res.mat.1.6[4,7] <- round(mean(-ate.sum.2a.6[,2]+ate.sum.2a.6[,3]),2)
res.mat.1.6[5,7] <- round(mean(-ate.sum.2b.6[,2]+ate.sum.2b.6[,3]),2)
res.mat.1.6[6,7] <- round(mean(-ate.sum.3a.6[,2]+ate.sum.3a.6[,3]),2)
res.mat.1.6[7,7] <- round(mean(-ate.sum.3b.6[,2]+ate.sum.3b.6[,3]),2)
##### Table 2
res.mat.2.6[2,7] <- round(mean(-ate.sum.4a.6[,2]+ate.sum.4a.6[,3]),2)
res.mat.2.6[3,7] <- round(mean(-ate.sum.4b.6[,2]+ate.sum.4b.6[,3]),2)
res.mat.2.6[4,7] <- round(mean(-ate.sum.5a.6[,2]+ate.sum.5a.6[,3]),2)
res.mat.2.6[5,7] <- round(mean(-ate.sum.5b.6[,2]+ate.sum.5b.6[,3]),2)
res.mat.2.6[6,7] <- round(mean(-ate.sum.6a.6[,2]+ate.sum.6a.6[,3]),2)
res.mat.2.6[7,7] <- round(mean(-ate.sum.6b.6[,2]+ate.sum.6b.6[,3]),2)
## P=10
### DGP1a
results <- readFiles(s,"1","xbcf",10,"a")
#### Obtain ATE posterior
ate.sum.1a.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP1b
results <- readFiles(s,"1","xbcf",10,"b")
#### Obtain ATE posterior
ate.sum.1b.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP2a
results <- readFiles(s,"2","xbcf",10,"a")
#### Obtain ATE posterior
ate.sum.2a.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP2b
results <- readFiles(s,"2","xbcf",10,"b")
#### Obtain ATE posterior
ate.sum.2b.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP3a
results <- readFiles(s,"3","xbcf",10,"a")
#### Obtain ATE posterior
ate.sum.3a.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP3b
results <- readFiles(s,"3","xbcf",10,"b")
#### Obtain ATE posterior
ate.sum.3b.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP4a
results <- readFiles(s,"4","xbcf",10,"a")
#### Obtain ATE posterior
ate.sum.4a.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP4b
results <- readFiles(s,"4","xbcf",10,"b")
#### Obtain ATE posterior
ate.sum.4b.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP5a
results <- readFiles(s,"5","xbcf",10,"a")
#### Obtain ATE posterior
ate.sum.5a.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP5b
results <- readFiles(s,"5","xbcf",10,"b")
#### Obtain ATE posterior
ate.sum.5b.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP6a
results <- readFiles(s,"6","xbcf",10,"a")
#### Obtain ATE posterior
ate.sum.6a.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### DGP6b
results <- readFiles(s,"6","xbcf",10,"b")
#### Obtain ATE posterior
ate.sum.6b.10 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
### Read data and save true ATE
data <- readRDS("Data/DGP1_10.rds")
data <- data[1:s]
ate1 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP2_10.rds")
data <- data[1:s]
ate2 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP3_10.rds")
data <- data[1:s]
ate3 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP4_10.rds")
data <- data[1:s]
ate4 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP5_10.rds")
data <- data[1:s]
ate5 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP6_10.rds")
data <- data[1:s]
ate6 <- sapply(data,function(x) x$ate)
### Store results
res.mat.1.10 <- matrix(NA,7,7)
res.mat.2.10 <- matrix(NA,7,7)
res.mat.1.10[1,] <- res.mat.2.10[1,] <- c("","","Controls","MSE","Variance","Coverage","Size")
res.mat.1.10[2,] <- c("XBCF","DGP 1","","","","","")
res.mat.2.10[2,] <- c("XBCF","DGP 4","","","","","")
res.mat.1.10[4,] <- c("XBCF","DGP 2","","","","","")
res.mat.2.10[4,] <- c("XBCF","DGP 5","","","","","")
res.mat.1.10[6,] <- c("XBCF","DGP 3","","","","","")
res.mat.2.10[6,] <- c("XBCF","DGP 6","","","","","")
res.mat.1.10[2:7,3] <- res.mat.2.10[2:7,3] <- c("No","Yes")
#### MSE
##### Table 1
res.mat.1.10[2,4] <- round(mean((ate.sum.1a.10[,1]-ate1)^2),3)
res.mat.1.10[3,4] <- round(mean((ate.sum.1b.10[,1]-ate1)^2),3)
res.mat.1.10[4,4] <- round(mean((ate.sum.2a.10[,1]-ate2)^2),3)
res.mat.1.10[5,4] <- round(mean((ate.sum.2b.10[,1]-ate2)^2),3)
res.mat.1.10[6,4] <- round(mean((ate.sum.3a.10[,1]-ate3)^2),3)
res.mat.1.10[7,4] <- round(mean((ate.sum.3b.10[,1]-ate3)^2),3)
##### Table 2
res.mat.2.10[2,4] <- round(mean((ate.sum.4a.10[,1]-ate4)^2),3)
res.mat.2.10[3,4] <- round(mean((ate.sum.4b.10[,1]-ate4)^2),3)
res.mat.2.10[4,4] <- round(mean((ate.sum.5a.10[,1]-ate5)^2),3)
res.mat.2.10[5,4] <- round(mean((ate.sum.5b.10[,1]-ate5)^2),3)
res.mat.2.10[6,4] <- round(mean((ate.sum.6a.10[,1]-ate6)^2),3)
res.mat.2.10[7,4] <- round(mean((ate.sum.6b.10[,1]-ate6)^2),3)
#### Variance
##### Table 1
res.mat.1.10[2,5] <- round(var(ate.sum.1a.10[,1]),3)
res.mat.1.10[3,5] <- round(var(ate.sum.1b.10[,1]),3)
res.mat.1.10[4,5] <- round(var(ate.sum.2a.10[,1]),3)
res.mat.1.10[5,5] <- round(var(ate.sum.2b.10[,1]),3)
res.mat.1.10[6,5] <- round(var(ate.sum.3a.10[,1]),3)
res.mat.1.10[7,5] <- round(var(ate.sum.3b.10[,1]),3)
##### Table 2
res.mat.2.10[2,5] <- round(var(ate.sum.4a.10[,1]),3)
res.mat.2.10[3,5] <- round(var(ate.sum.4b.10[,1]),3)
res.mat.2.10[4,5] <- round(var(ate.sum.5a.10[,1]),3)
res.mat.2.10[5,5] <- round(var(ate.sum.5b.10[,1]),3)
res.mat.2.10[6,5] <- round(var(ate.sum.6a.10[,1]),3)
res.mat.2.10[7,5] <- round(var(ate.sum.6b.10[,1]),3)
#### Coverage
##### Table 1
res.mat.1.10[2,6] <- 100*mean(ate.sum.1a.10[,2]<=ate1 & ate1<=ate.sum.1a.10[,3])
res.mat.1.10[3,6] <- 100*mean(ate.sum.1b.10[,2]<=ate1 & ate1<=ate.sum.1b.10[,3])
res.mat.1.10[4,6] <- 100*mean(ate.sum.2a.10[,2]<=ate2 & ate2<=ate.sum.2a.10[,3])
res.mat.1.10[5,6] <- 100*mean(ate.sum.2b.10[,2]<=ate2 & ate2<=ate.sum.2b.10[,3])
res.mat.1.10[6,6] <- 100*mean(ate.sum.3a.10[,2]<=ate3 & ate3<=ate.sum.3a.10[,3])
res.mat.1.10[7,6] <- 100*mean(ate.sum.3b.10[,2]<=ate3 & ate3<=ate.sum.3b.10[,3])
##### Table 2
res.mat.2.10[2,6] <- 100*mean(ate.sum.4a.10[,2]<=ate4 & ate4<=ate.sum.4a.10[,3])
res.mat.2.10[3,6] <- 100*mean(ate.sum.4b.10[,2]<=ate4 & ate4<=ate.sum.4b.10[,3])
res.mat.2.10[4,6] <- 100*mean(ate.sum.5a.10[,2]<=ate5 & ate5<=ate.sum.5a.10[,3])
res.mat.2.10[5,6] <- 100*mean(ate.sum.5b.10[,2]<=ate5 & ate5<=ate.sum.5b.10[,3])
res.mat.2.10[6,6] <- 100*mean(ate.sum.6a.10[,2]<=ate6 & ate6<=ate.sum.6a.10[,3])
res.mat.2.10[7,6] <- 100*mean(ate.sum.6b.10[,2]<=ate6 & ate6<=ate.sum.6b.10[,3])
#### Size
##### Table 1
res.mat.1.10[2,7] <- round(mean(-ate.sum.1a.10[,2]+ate.sum.1a.10[,3]),2)
res.mat.1.10[3,7] <- round(mean(-ate.sum.1b.10[,2]+ate.sum.1b.10[,3]),2)
res.mat.1.10[4,7] <- round(mean(-ate.sum.2a.10[,2]+ate.sum.2a.10[,3]),2)
res.mat.1.10[5,7] <- round(mean(-ate.sum.2b.10[,2]+ate.sum.2b.10[,3]),2)
res.mat.1.10[6,7] <- round(mean(-ate.sum.3a.10[,2]+ate.sum.3a.10[,3]),2)
res.mat.1.10[7,7] <- round(mean(-ate.sum.3b.10[,2]+ate.sum.3b.10[,3]),2)
##### Table 2
res.mat.2.10[2,7] <- round(mean(-ate.sum.4a.10[,2]+ate.sum.4a.10[,3]),2)
res.mat.2.10[3,7] <- round(mean(-ate.sum.4b.10[,2]+ate.sum.4b.10[,3]),2)
res.mat.2.10[4,7] <- round(mean(-ate.sum.5a.10[,2]+ate.sum.5a.10[,3]),2)
res.mat.2.10[5,7] <- round(mean(-ate.sum.5b.10[,2]+ate.sum.5b.10[,3]),2)
res.mat.2.10[6,7] <- round(mean(-ate.sum.6a.10[,2]+ate.sum.6a.10[,3]),2)
res.mat.2.10[7,7] <- round(mean(-ate.sum.6b.10[,2]+ate.sum.6b.10[,3]),2)
## Heterogeneous effects
### Read data and save true ATE and CATE
data <- readRDS("Data/DGP7.rds")
ate7 <- sapply(data,function(x) x$ate)
ate7 <- ate7[1:s]
### Obtain ATE and CATE posterior
results <- readFiles2(s,"7","xbcf")
ate.sum.7 <- t(sapply(results, function(x) c(mean(colMeans(x$ate.post)),quantile(colMeans(x$ate.post),c(0.025,0.975)))))
cate.sum.7 <- lapply(results, function(i) t(apply(i$ate.post,1,function(x) c(mean(x),quantile(x,c(0.025,0.975))))))
cate <- mapply(function(i,j) i$cate[-j$Owidth<=i$x & i$x<=j$Owidth],data,results,SIMPLIFY=F)[1:s]
### Table results
res.mat.het <- matrix("",3,9)
res.mat.het[1,] <- c("","","ATE","","","","CATE","","")
res.mat.het[2,] <- c("","","MSE","Coverage","Size","","MSE","Coverage","Size")
res.mat.het[3,1] <- "XBCF"
res.mat.het[3,3] <- round(mean((ate.sum.7[,1]-ate7)^2),2)
res.mat.het[3,4] <- mean(ate.sum.7[,2]<=ate7 & ate7<=ate.sum.7[,3])
res.mat.het[3,5] <- round(mean(-ate.sum.7[,2]+ate.sum.7[,3]),2)
res.mat.het[3,7] <- round(mean(mapply(function(x,y) mean((x[,1]-y)^2),cate.sum.7,cate)),2)
res.mat.het[3,8] <- round(mean(mapply(function(x,y) mean(x[,2]<=y & y<=x[,3]),cate.sum.7,cate)),2)
res.mat.het[3,9] <- round(mean(sapply(cate.sum.7,function(i) mean(i[,3]-i[,2]))),2)
###
rm(data,results,no_cores,readFiles,readFiles2,ate1,ate2,ate3,ate4,ate5,ate6,ate7,cate.sum.7,cate)
save.image("Tables/xbcf_results.RData")
####
stopImplicitCluster()
####
