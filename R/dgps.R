devtools::install_github("JingyuHe/XBART@XBCF-RDD")
library(XBART)
## Read the data
data1 <- readRDS("Data/DGP1_4.rds")
data1 <- data1[c(1,10,100,1000)]
data2 <- readRDS("Data/DGP2_4.rds")
data2 <- data2[c(1,10,100,1000)]
data3 <- readRDS("Data/DGP3_4.rds")
data3 <- data3[c(1,10,100,1000)]
data4 <- readRDS("Data/DGP4_4.rds")
data4 <- data4[c(1,10,100,1000)]
data5 <- readRDS("Data/DGP5_4.rds")
data5 <- data5[c(1,10,100,1000)]
data6 <- readRDS("Data/DGP6_4.rds")
data6 <- data6[c(1,10,100,1000)]
###
n <- 500
c             <- 0
Omin          <- 10
Opct          <- 0.9
m             <- 10
Nmin          <- 10
num_sweeps    <- 150
burnin        <- 50
p_categorical <- 0
num_cutpoints <- n
attach(data4[[3]])
h <- 0.2
fit <- XBCF.rd(y, w, x, c, Owidth = h, Omin = Omin, Opct = Opct,
               num_trees_mod = m, num_trees_con = m,
               num_cutpoints = n, num_sweeps = num_sweeps,
               burnin = burnin, Nmin = Nmin,
               p_categorical_con = p_categorical, p_categorical_mod = p_categorical,
               tau_con = 2*var(y)/m,
               tau_mod = 0.5*var(y)/m, parallel=F)
test <- -h<=x & x<=h
pred <- predict.XBCFrd(fit,w[test,],rep(0,sum(test)))
post <- pred$tau.adj[,(burnin+1):num_sweeps]
summary(colMeans(post))
fit$importance_treatment
plot(cumsum(colMeans(fit$sigma0))/1:150,type="l")
