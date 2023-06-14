setwd("~/Documents/Git/XBCF-RDD")
## Main tables
### Previous simulations
xbcf <- read.csv("Tables/xbcf_simulations_1.csv",col.names=c("","DGP1a","DGP1b","DGP2","DGP3"))
kr <- read.csv("Tables/kr_simulations_1.csv",col.names=c("","DGP1a","DGP1b","DGP2","DGP3"))
cgs <- read.csv("Tables/cgs_simulations_1.csv",col.names=c("","DGP1a","DGP1b","DGP2","DGP3"))
###
sim1 <- cbind(xbcf[,2],kr[,2],cgs[,2],
              xbcf[,3],kr[,3],cgs[,3],
              xbcf[,4],kr[,4],cgs[,4],
              xbcf[,5],kr[,5],cgs[,5])
sim1[1,] <- round(sim1[1,],3)
sim1[2,] <- sim1[2,]*100
sim1[2,] <- round(sim1[2,],1)
sim1[3,] <- round(sim1[3,],3)
sim1 <- rbind(c("XBCF","KR","CGS"),sim1)
sim1 <- rbind(c("","1(a)","","","1(b)","",
                "","2","","","3",""),sim1)
rownames(sim1) <- c("","","MSE","Coverage","Length")
###
print(xtable::xtable(sim1,caption="Simulations",label="tab:sim1"), include.colnames = F)
### New simulations
xbcf <- read.csv("Tables/xbcf_simulations_2.csv",col.names=c("","DGP4","DGP5","DGP6"))
kr <- read.csv("Tables/kr_simulations_2.csv",col.names=c("","DGP4","DGP5","DGP6"))
cgs <- read.csv("Tables/cgs_simulations_2.csv",col.names=c("","DGP4","DGP5","DGP6"))
###
sim2 <- cbind(xbcf[,2],kr[,2],cgs[,2],
              xbcf[,3],kr[,3],cgs[,3],
              xbcf[,4],kr[,4],cgs[,4])
sim2[1,] <- round(sim2[1,],3)
sim2[2,] <- sim2[2,]*100
sim2[2,] <- round(sim2[2,],1)
sim2[3,] <- round(sim2[3,],3)
sim2 <- rbind(c("XBCF","KR","CGS"),sim2)
sim2 <- rbind(c("","4","","","5","","","6",""),sim2)
rownames(sim2) <- c("","","MSE","Coverage","Length")
###
print(xtable::xtable(sim2,caption="Simulations",label="tab:sim2"), include.colnames = F)
## Plot of point estimates
xbcf <- readRDS("Tables/xbcf_simulations.rds")
kr <- readRDS("Tables/kr_simulations.rds")
cgs <- readRDS("Tables/cgs_simulations.rds")
box1 <- data.frame(tau=unlist(lapply(xbcf, function(i) i[,1])),DGP=c(rep("1(a)",1000),rep("1(b)",1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000),rep(6,1000)),Est="XBCF")
box2 <- data.frame(tau=unlist(lapply(kr, function(i) i[,1])),DGP=c(rep("1(a)",1000),rep("1(b)",1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000),rep(6,1000)),Est="KR")
box3 <- data.frame(tau=unlist(lapply(cgs, function(i) i[,1])),DGP=c(rep("1(a)",1000),rep("1(b)",1000),rep(2,1000),rep(3,1000),rep(4,1000),rep(5,1000),rep(6,1000)),Est="CGS")
box <- rbind(box1,box2)
box <- rbind(box,box3)
box.pre <- subset(box,DGP %in% c("1(a)","1(b)","2","3"))
box.new <- subset(box,DGP %in% c("4","5","6"))
###
cols <- rainbow(3,s=0.5)
png("Figures/point_estimates.png")
par(mfrow=c(1,2))
boxplot(tau~Est+DGP,data=box.pre, col = cols, axes=F,
        at = c(1:3,5:7,9:11,13:15),
        ylab = expression(tau), xlab = "DGP")
axis(2)
axis(1, at=c(2,6,10,14), labels = c("1(a)","1(b)","2","3"))
legend("topleft",c("CGS","KR","XBCF"),fill=cols, cex=0.6, ncol=2)
boxplot(tau~Est+DGP,data=box.new, col = cols, axes=F,
        at = c(1:3,5:7,9:11),
        ylab = expression(tau), xlab = "DGP")
axis(2)
axis(1, at=c(2,6,10), labels = c("4","5","6"))
legend("topleft",c("CGS","KR","XBCF"),fill=cols, cex=0.6, ncol=2)
par(mfrow=c(1,1))
dev.off()
## Plot of bias and variance
dgp1 <- readRDS("Data/DGP1.rds")
ate1 <- sapply(dgp1, function(i) i$ate)
dgp2 <- readRDS("Data/DGP2.rds")
ate2 <- sapply(dgp2, function(i) i$ate)
dgp3 <- readRDS("Data/DGP3.rds")
ate3 <- sapply(dgp3, function(i) i$ate)
dgp4 <- readRDS("Data/DGP4.rds")
ate4 <- sapply(dgp4, function(i) i$ate)
dgp5 <- readRDS("Data/DGP5.rds")
ate5 <- sapply(dgp5, function(i) i$ate)
dgp6 <- readRDS("Data/DGP6.rds")
ate6 <- sapply(dgp6, function(i) i$ate)
###
bv <- data.frame(Est=c(rep("XBCF",7),rep("KR",7),rep("CGS",7)),
                 DGP=rep(c("1(a)","1(b)","2","3","4","5","6"),3),
                 Bias=0,Var=0)
#### XBCF
bv[1,3:4] <- c(mean(xbcf[[1]][,1]-ate1)^2,var(xbcf[[1]][,1]))
bv[2,3:4] <- c(mean(xbcf[[2]][,1]-ate1)^2,var(xbcf[[2]][,1]))
bv[3,3:4] <- c(mean(xbcf[[3]][,1]-ate2)^2,var(xbcf[[3]][,1]))
bv[4,3:4] <- c(mean(xbcf[[4]][,1]-ate3)^2,var(xbcf[[4]][,1]))
bv[5,3:4] <- c(mean(xbcf[[5]][,1]-ate4)^2,var(xbcf[[5]][,1]))
bv[6,3:4] <- c(mean(xbcf[[6]][,1]-ate5)^2,var(xbcf[[6]][,1]))
bv[7,3:4] <- c(mean(xbcf[[7]][,1]-ate6)^2,var(xbcf[[7]][,1]))
#### KR
bv[8,3:4] <- c(mean(kr[[1]][,1]-ate1)^2,var(kr[[1]][,1]))
bv[9,3:4] <- c(mean(kr[[2]][,1]-ate1)^2,var(kr[[2]][,1]))
bv[10,3:4] <- c(mean(kr[[3]][,1]-ate2)^2,var(kr[[3]][,1]))
bv[11,3:4] <- c(mean(kr[[4]][,1]-ate3)^2,var(kr[[4]][,1]))
bv[12,3:4] <- c(mean(kr[[5]][,1]-ate4)^2,var(kr[[5]][,1]))
bv[13,3:4] <- c(mean(kr[[6]][,1]-ate5)^2,var(kr[[6]][,1]))
bv[14,3:4] <- c(mean(kr[[7]][,1]-ate6)^2,var(kr[[7]][,1]))
#### CGS
bv[15,3:4] <- c(mean(cgs[[1]][,1]-ate1)^2,var(cgs[[1]][,1]))
bv[16,3:4] <- c(mean(cgs[[2]][,1]-ate1)^2,var(cgs[[2]][,1]))
bv[17,3:4] <- c(mean(cgs[[3]][,1]-ate2)^2,var(cgs[[3]][,1]))
bv[18,3:4] <- c(mean(cgs[[4]][,1]-ate3)^2,var(cgs[[4]][,1]))
bv[19,3:4] <- c(mean(cgs[[5]][,1]-ate4)^2,var(cgs[[5]][,1]))
bv[20,3:4] <- c(mean(cgs[[6]][,1]-ate5)^2,var(cgs[[6]][,1]))
bv[21,3:4] <- c(mean(cgs[[7]][,1]-ate6)^2,var(cgs[[7]][,1]))
###
png("Figures/bias_variance.png")
par(mfrow=c(1,2))
barplot(bv$Bias~bv$Est+bv$DGP,beside=T,col=cols,
        xlab="DGP",ylab="Squared Bias")
legend(x="topleft",legend=c("CGS","KR","XBCF"),fill=cols,
       cex=0.6, ncol=2)
barplot(bv$Var~bv$Est+bv$DGP,beside=T,col=cols,
        xlab="DGP",ylab="Variance")
legend(x="topleft",legend=c("CGS","KR","XBCF"),fill=cols,
       cex=0.6, ncol=2)
par(mfrow=c(1,1))
dev.off()
