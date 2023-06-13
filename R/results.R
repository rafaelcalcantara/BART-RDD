setwd("~/Documents/Git/XBCF-RDD")
## Previous simulations
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
## New simulations
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
