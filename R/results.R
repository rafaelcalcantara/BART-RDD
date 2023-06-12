setwd("~/Documents/Git/XBCF-RDD")
xbcf <- read.csv("Tables/xbcf_simulations_1.csv",col.names=c("DGP1a","DGP1b","DGP2","DGP3"))
kr <- read.csv("Tables/kr_simulations_1.csv",col.names=c("DGP1a","DGP1b","DGP2","DGP3"))
cgs <- read.csv("Tables/cgs_simulations_1.csv",col.names=c("DGP1a","DGP1b","DGP2","DGP3"))
###
sim1 <- cbind(xbcf[,1],kr[,1],cgs[,1],
              xbcf[,2],kr[,2],cgs[,2],
              xbcf[,3],kr[,3],cgs[,3],
              xbcf[,4],kr[,4],cgs[,4])
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
