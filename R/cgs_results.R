## library(XBART)
library(parallel)
library(foreach)
library(doParallel)
setwd("~/Git/XBCF-RDD/")
s <- 1000
### Function to read results files
readFiles <- function(s,dgp,file,p,ab)
{
    foreach(x=1:s,.multicombine=T) %do%
        readRDS(paste0("Results/",file,dgp,"_",ab,"_",p,"_",x,".rds"))
}
#### Parallelization
no_cores <- detectCores() - 1
registerDoParallel(no_cores)
## P=4
### DGP1a
results <- readFiles(s,"1","cgs",4,"a")
#### Obtain ATE posterior
ate.sum.1a.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP1b
results <- readFiles(s,"1","cgs",4,"b")
#### Obtain ATE posterior
ate.sum.1b.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP2a
results <- readFiles(s,"2","cgs",4,"a")
#### Obtain ATE posterior
ate.sum.2a.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP2b
results <- readFiles(s,"2","cgs",4,"b")
#### Obtain ATE posterior
ate.sum.2b.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP3a
results <- readFiles(s,"3","cgs",4,"a")
#### Obtain ATE posterior
ate.sum.3a.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP3b
results <- readFiles(s,"3","cgs",4,"b")
#### Obtain ATE posterior
ate.sum.3b.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP4a
results <- readFiles(s,"4","cgs",4,"a")
#### Obtain ATE posterior
ate.sum.4a.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP4b
results <- readFiles(s,"4","cgs",4,"b")
#### Obtain ATE posterior
ate.sum.4b.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP5a
results <- readFiles(s,"5","cgs",4,"a")
#### Obtain ATE posterior
ate.sum.5a.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP5b
results <- readFiles(s,"5","cgs",4,"b")
#### Obtain ATE posterior
ate.sum.5b.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP6a
results <- readFiles(s,"6","cgs",4,"a")
#### Obtain ATE posterior
ate.sum.6a.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP6b
results <- readFiles(s,"6","cgs",4,"b")
#### Obtain ATE posterior
ate.sum.6b.4 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### Read data and save true ATE
data <- readRDS("Data/DGP1_4.rds")
ate1 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP2_4.rds")
ate2 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP3_4.rds")
ate3 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP4_4.rds")
ate4 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP5_4.rds")
ate5 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP6_4.rds")
ate6 <- sapply(data,function(x) x$ate)
### Store results
res.mat.1.4 <- matrix(NA,7,7)
res.mat.2.4 <- matrix(NA,7,7)
res.mat.1.4[1,] <- res.mat.2.4[1,] <- c("","","Controls","MSE","Variance","Coverage","Size")
res.mat.1.4[2,] <- res.mat.2.4[2,] <- c("CGS","DGP 1","","","","","")
res.mat.1.4[4,] <- res.mat.2.4[4,] <- c("CGS","DGP 2","","","","","")
res.mat.1.4[6,] <- res.mat.2.4[6,] <- c("CGS","DGP 3","","","","","")
res.mat.1.4[2:7,3] <- res.mat.2.4[2:7,3] <- c("No","Yes")
#### MSE
##### Table 1
res.mat.1.4[2,4] <- round(mean((ate.sum.1a.4[,1]-ate1)^2),2)
res.mat.1.4[3,4] <- round(mean((ate.sum.1b.4[,1]-ate1)^2),2)
res.mat.1.4[4,4] <- round(mean((ate.sum.2a.4[,1]-ate2)^2),2)
res.mat.1.4[5,4] <- round(mean((ate.sum.2b.4[,1]-ate2)^2),2)
res.mat.1.4[6,4] <- round(mean((ate.sum.3a.4[,1]-ate3)^2),2)
res.mat.1.4[7,4] <- round(mean((ate.sum.3b.4[,1]-ate3)^2),2)
##### Table 2
res.mat.2.4[2,4] <- round(mean((ate.sum.4a.4[,1]-ate4)^2),2)
res.mat.2.4[3,4] <- round(mean((ate.sum.4b.4[,1]-ate4)^2),2)
res.mat.2.4[4,4] <- round(mean((ate.sum.5a.4[,1]-ate5)^2),2)
res.mat.2.4[5,4] <- round(mean((ate.sum.5b.4[,1]-ate5)^2),2)
res.mat.2.4[6,4] <- round(mean((ate.sum.6a.4[,1]-ate6)^2),2)
res.mat.2.4[7,4] <- round(mean((ate.sum.6b.4[,1]-ate6)^2),2)
#### Variance
##### Table 1
res.mat.1.4[2,5] <- round(var(ate.sum.1a.4[,1]),2)
res.mat.1.4[3,5] <- round(var(ate.sum.1b.4[,1]),2)
res.mat.1.4[4,5] <- round(var(ate.sum.2a.4[,1]),2)
res.mat.1.4[5,5] <- round(var(ate.sum.2b.4[,1]),2)
res.mat.1.4[6,5] <- round(var(ate.sum.3a.4[,1]),2)
res.mat.1.4[7,5] <- round(var(ate.sum.3b.4[,1]),2)
##### Table 2
res.mat.2.4[2,5] <- round(var(ate.sum.4a.4[,1]),2)
res.mat.2.4[3,5] <- round(var(ate.sum.4b.4[,1]),2)
res.mat.2.4[4,5] <- round(var(ate.sum.5a.4[,1]),2)
res.mat.2.4[5,5] <- round(var(ate.sum.5b.4[,1]),2)
res.mat.2.4[6,5] <- round(var(ate.sum.6a.4[,1]),2)
res.mat.2.4[7,5] <- round(var(ate.sum.6b.4[,1]),2)
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
results <- readFiles(s,"1","cgs",6,"a")
#### Obtain ATE posterior
ate.sum.1a.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP1b
results <- readFiles(s,"1","cgs",6,"b")
#### Obtain ATE posterior
ate.sum.1b.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP2a
results <- readFiles(s,"2","cgs",6,"a")
#### Obtain ATE posterior
ate.sum.2a.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP2b
results <- readFiles(s,"2","cgs",6,"b")
#### Obtain ATE posterior
ate.sum.2b.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP3a
results <- readFiles(s,"3","cgs",6,"a")
#### Obtain ATE posterior
ate.sum.3a.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP3b
results <- readFiles(s,"3","cgs",6,"b")
#### Obtain ATE posterior
ate.sum.3b.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP4a
results <- readFiles(s,"4","cgs",6,"a")
#### Obtain ATE posterior
ate.sum.4a.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP4b
results <- readFiles(s,"4","cgs",6,"b")
#### Obtain ATE posterior
ate.sum.4b.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP5a
results <- readFiles(s,"5","cgs",6,"a")
#### Obtain ATE posterior
ate.sum.5a.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP5b
results <- readFiles(s,"5","cgs",6,"b")
#### Obtain ATE posterior
ate.sum.5b.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP6a
results <- readFiles(s,"6","cgs",6,"a")
#### Obtain ATE posterior
ate.sum.6a.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP6b
results <- readFiles(s,"6","cgs",6,"b")
#### Obtain ATE posterior
ate.sum.6b.6 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### Read data and save true ATE
data <- readRDS("Data/DGP1_6.rds")
ate1 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP2_6.rds")
ate2 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP3_6.rds")
ate3 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP4_6.rds")
ate4 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP5_6.rds")
ate5 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP6_6.rds")
ate6 <- sapply(data,function(x) x$ate)
### Store results
res.mat.1.6 <- matrix(NA,7,7)
res.mat.2.6 <- matrix(NA,7,7)
res.mat.1.6[1,] <- res.mat.2.6[1,] <- c("","","Controls","MSE","Variance","Coverage","Size")
res.mat.1.6[2,] <- res.mat.2.6[2,] <- c("CGS","DGP 1","","","","","")
res.mat.1.6[4,] <- res.mat.2.6[4,] <- c("CGS","DGP 2","","","","","")
res.mat.1.6[6,] <- res.mat.2.6[6,] <- c("CGS","DGP 3","","","","","")
res.mat.1.6[2:7,3] <- res.mat.2.6[2:7,3] <- c("No","Yes")
#### MSE
##### Table 1
res.mat.1.6[2,4] <- round(mean((ate.sum.1a.6[,1]-ate1)^2),2)
res.mat.1.6[3,4] <- round(mean((ate.sum.1b.6[,1]-ate1)^2),2)
res.mat.1.6[4,4] <- round(mean((ate.sum.2a.6[,1]-ate2)^2),2)
res.mat.1.6[5,4] <- round(mean((ate.sum.2b.6[,1]-ate2)^2),2)
res.mat.1.6[6,4] <- round(mean((ate.sum.3a.6[,1]-ate3)^2),2)
res.mat.1.6[7,4] <- round(mean((ate.sum.3b.6[,1]-ate3)^2),2)
##### Table 2
res.mat.2.6[2,4] <- round(mean((ate.sum.4a.6[,1]-ate4)^2),2)
res.mat.2.6[3,4] <- round(mean((ate.sum.4b.6[,1]-ate4)^2),2)
res.mat.2.6[4,4] <- round(mean((ate.sum.5a.6[,1]-ate5)^2),2)
res.mat.2.6[5,4] <- round(mean((ate.sum.5b.6[,1]-ate5)^2),2)
res.mat.2.6[6,4] <- round(mean((ate.sum.6a.6[,1]-ate6)^2),2)
res.mat.2.6[7,4] <- round(mean((ate.sum.6b.6[,1]-ate6)^2),2)
#### Variance
##### Table 1
res.mat.1.6[2,5] <- round(var(ate.sum.1a.6[,1]),2)
res.mat.1.6[3,5] <- round(var(ate.sum.1b.6[,1]),2)
res.mat.1.6[4,5] <- round(var(ate.sum.2a.6[,1]),2)
res.mat.1.6[5,5] <- round(var(ate.sum.2b.6[,1]),2)
res.mat.1.6[6,5] <- round(var(ate.sum.3a.6[,1]),2)
res.mat.1.6[7,5] <- round(var(ate.sum.3b.6[,1]),2)
##### Table 2
res.mat.2.6[2,5] <- round(var(ate.sum.4a.6[,1]),2)
res.mat.2.6[3,5] <- round(var(ate.sum.4b.6[,1]),2)
res.mat.2.6[4,5] <- round(var(ate.sum.5a.6[,1]),2)
res.mat.2.6[5,5] <- round(var(ate.sum.5b.6[,1]),2)
res.mat.2.6[6,5] <- round(var(ate.sum.6a.6[,1]),2)
res.mat.2.6[7,5] <- round(var(ate.sum.6b.6[,1]),2)
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
results <- readFiles(s,"1","cgs",10,"a")
#### Obtain ATE posterior
ate.sum.1a.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP1b
results <- readFiles(s,"1","cgs",10,"b")
#### Obtain ATE posterior
ate.sum.1b.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP2a
results <- readFiles(s,"2","cgs",10,"a")
#### Obtain ATE posterior
ate.sum.2a.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP2b
results <- readFiles(s,"2","cgs",10,"b")
#### Obtain ATE posterior
ate.sum.2b.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP3a
results <- readFiles(s,"3","cgs",10,"a")
#### Obtain ATE posterior
ate.sum.3a.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP3b
results <- readFiles(s,"3","cgs",10,"b")
#### Obtain ATE posterior
ate.sum.3b.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP4a
results <- readFiles(s,"4","cgs",10,"a")
#### Obtain ATE posterior
ate.sum.4a.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP4b
results <- readFiles(s,"4","cgs",10,"b")
#### Obtain ATE posterior
ate.sum.4b.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP5a
results <- readFiles(s,"5","cgs",10,"a")
#### Obtain ATE posterior
ate.sum.5a.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP5b
results <- readFiles(s,"5","cgs",10,"b")
#### Obtain ATE posterior
ate.sum.5b.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP6a
results <- readFiles(s,"6","cgs",10,"a")
#### Obtain ATE posterior
ate.sum.6a.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### DGP6b
results <- readFiles(s,"6","cgs",10,"b")
#### Obtain ATE posterior
ate.sum.6b.10 <- t(sapply(results, function(x) c(mean(x$ate.post),quantile(x$ate.post,c(0.025,0.975)))))
### Read data and save true ATE
data <- readRDS("Data/DGP1_10.rds")
ate1 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP2_10.rds")
ate2 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP3_10.rds")
ate3 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP4_10.rds")
ate4 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP5_10.rds")
ate5 <- sapply(data,function(x) x$ate)
data <- readRDS("Data/DGP6_10.rds")
ate6 <- sapply(data,function(x) x$ate)
### Store results
res.mat.1.10 <- matrix(NA,7,7)
res.mat.2.10 <- matrix(NA,7,7)
res.mat.1.10[1,] <- res.mat.2.10[1,] <- c("","","Controls","MSE","Variance","Coverage","Size")
res.mat.1.10[2,] <- c("CGS","DGP 1","","","","","")
res.mat.2.10[2,] <- c("CGS","DGP 4","","","","","")
res.mat.1.10[4,] <- c("CGS","DGP 2","","","","","")
res.mat.2.10[4,] <- c("CGS","DGP 5","","","","","")
res.mat.1.10[6,] <- c("CGS","DGP 3","","","","","")
res.mat.2.10[6,] <- c("CGS","DGP 6","","","","","")
res.mat.1.10[2:7,3] <- res.mat.2.10[2:7,3] <- c("No","Yes")
#### MSE
##### Table 1
res.mat.1.10[2,4] <- round(mean((ate.sum.1a.10[,1]-ate1)^2),2)
res.mat.1.10[3,4] <- round(mean((ate.sum.1b.10[,1]-ate1)^2),2)
res.mat.1.10[4,4] <- round(mean((ate.sum.2a.10[,1]-ate2)^2),2)
res.mat.1.10[5,4] <- round(mean((ate.sum.2b.10[,1]-ate2)^2),2)
res.mat.1.10[6,4] <- round(mean((ate.sum.3a.10[,1]-ate3)^2),2)
res.mat.1.10[7,4] <- round(mean((ate.sum.3b.10[,1]-ate3)^2),2)
##### Table 2
res.mat.2.10[2,4] <- round(mean((ate.sum.4a.10[,1]-ate4)^2),2)
res.mat.2.10[3,4] <- round(mean((ate.sum.4b.10[,1]-ate4)^2),2)
res.mat.2.10[4,4] <- round(mean((ate.sum.5a.10[,1]-ate5)^2),2)
res.mat.2.10[5,4] <- round(mean((ate.sum.5b.10[,1]-ate5)^2),2)
res.mat.2.10[6,4] <- round(mean((ate.sum.6a.10[,1]-ate6)^2),2)
res.mat.2.10[7,4] <- round(mean((ate.sum.6b.10[,1]-ate6)^2),2)
#### Variance
##### Table 1
res.mat.1.10[2,5] <- round(var(ate.sum.1a.10[,1]),2)
res.mat.1.10[3,5] <- round(var(ate.sum.1b.10[,1]),2)
res.mat.1.10[4,5] <- round(var(ate.sum.2a.10[,1]),2)
res.mat.1.10[5,5] <- round(var(ate.sum.2b.10[,1]),2)
res.mat.1.10[6,5] <- round(var(ate.sum.3a.10[,1]),2)
res.mat.1.10[7,5] <- round(var(ate.sum.3b.10[,1]),2)
##### Table 2
res.mat.2.10[2,5] <- round(var(ate.sum.4a.10[,1]),2)
res.mat.2.10[3,5] <- round(var(ate.sum.4b.10[,1]),2)
res.mat.2.10[4,5] <- round(var(ate.sum.5a.10[,1]),2)
res.mat.2.10[5,5] <- round(var(ate.sum.5b.10[,1]),2)
res.mat.2.10[6,5] <- round(var(ate.sum.6a.10[,1]),2)
res.mat.2.10[7,5] <- round(var(ate.sum.6b.10[,1]),2)
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
###
rm(data,results,no_cores,sreadFiles,ate1,ate2,ate3,ate4,ate5,ate6)
save.image("Tables/cgs_results.RData")
####
stopImplicitCluster()
####
