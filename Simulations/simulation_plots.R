if ("Figures" %in% list.files("Results") == FALSE) dir.create("Results/Figures/")
# Table 1----------------------------------------------------------------------
Owidth <- 0.1
res.list <- c("k1_1_k2_1_k3_0_k4_0.1_k5_0_p_2_rho_0.5",
              "k1_1_k2_1_k3_0_k4_0.1_k5_0_p_4_rho_0",
              "k1_1_k2_1_k3_0_k4_0.1_k5_1_p_2_rho_0",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_0_p_4_rho_0.5",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_1_p_2_rho_0.5",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_1_p_4_rho_0")
rmse.mean <- matrix(0,6,14)
colnames(rmse.mean) <- c("k1","k2","k3","k4","k5","p","rho",
                         "BARDDT","T-BART","S-BART","OLS","HS","Ridge","RD-Tree")
for (i in 1:length(res.list))
{
  dgp <- res.list[i]
  out <- matrix(0,100,8)
  colnames(out) <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate","reguly")
  for (sample in 1:100)
  {
    file <- paste0("sample_",sample,".csv")
    methods <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate")
    file.name <- paste0(paste("Results/RMSE/",dgp,"/",methods,sep=""),"_",file)
    data <- readRDS(paste0("Data/dgp_",dgp,"_",gsub("\\.csv","\\.rds",file)))
    cate <- data$cate[data$x>-Owidth & data$x<Owidth]
    for (j in 1:length(methods))
    {
      out[sample,methods[j]] <- as.numeric(read.table(file.name[j]))
    }
    file.reguly <- paste0("RD-Tree/Results/",dgp,"/",file)
    out.reguly <- as.matrix(read.table(file.reguly))
    out.reguly <- sqrt(mean((out.reguly[,1]-cate)^2))
    out[sample,"reguly"] <- out.reguly
    out[sample,] <- out[sample,]/out[sample,"ate"]
  }
  out <- out[,colnames(out) != "ate"]
  rmse.mean[i,] <- as.numeric(c(strsplit(dgp,"_")[[1]][(1:7)*2],round(colMeans(out),3)))
}
rmse.mean <- as.data.frame(rmse.mean)
print(xtable::xtable(subset(rmse.mean,(k1==1 & k2==1 & k3==0 & k4==0.1) | (k1==5 & k2==0.25 & k3==1 & k4==0.5)),
                     caption="Average RMSE per DGP, also divided here by the RMSE of the naive ATE estimator",
                     label="tab:results",
                     digits=2,
                     align=rep("c",1+ncol(rmse.mean))),
      include.rownames=FALSE)
# Figure 5---------------------------------------------------------------------
axis.font.size <- 0.7
estimators <- c("BARDDT","T-BART","OLS","HS","Ridge","RD-Tree")
angle <- 90
ind <- 0
pdf(paste0("Results/Figures/boxplots.pdf"))
par(mfrow=c(2,3),mai=c(0.7,0.4,0.4,0.1))
for (dgp in res.list[1:3])
{
  ind <- ind+1
  out <- matrix(0,100,8)
  colnames(out) <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate","reguly")
  for (sample in 1:100)
  {
    file <- paste0("sample_",sample,".csv")
    methods <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate")
    file.name <- paste0(paste("Results/RMSE/",dgp,"/",methods,sep=""),"_",file)
    for (j in 1:length(methods))
    {
      data <- readRDS(paste0("Data/dgp_",dgp,"_",gsub("\\.csv","\\.rds",file)))
      cate <- data$cate[data$x>-0.1 & data$x<0.1]
      ##
      out[sample,methods[j]] <- as.numeric(read.table(file.name[j]))
    }
    file.reguly <- paste0("RD-Tree/Results/",dgp,"/",file)
    out.reguly <- as.matrix(read.table(file.reguly))
    out.reguly <- sqrt(mean((out.reguly[,1]-cate)^2))
    out[sample,"reguly"] <- out.reguly
    out[sample,] <- out[sample,]/out[sample,"ate"]
  }
  out <- out[,colnames(out) != "ate"]
  colnames(out) <- c("BARDDT","T-BART","S-BART","OLS","HS","Ridge","RD-Tree")
  if (ind==3)
  {
    boxplot(out[,-3], main = "",frame=FALSE,ylim=c(0,1.3),xaxt="n")
    axis(1,1:length(estimators),labels=FALSE)
    text(x = 1:length(estimators), y = par("usr")[3] - 0.1, labels = estimators, srt = angle, adj = 1, xpd = TRUE, cex=axis.font.size)
  } else if (ind==2)
  {
    boxplot(out[,-3], main = "",frame=FALSE,ylim=c(0,1.3),xaxt="n")
    mtext(bquote("High signal, separable "~mu),line=1,font=2)
    axis(1,1:length(estimators),labels=FALSE)
    text(x = 1:length(estimators), y = par("usr")[3] - 0.1, labels = estimators, srt = angle, adj = 1, xpd = TRUE, cex=axis.font.size)
  } else
  {
    boxplot(out[,-3], main = "",frame=FALSE,ylim=c(0,1.3),xaxt="n")
    axis(1,1:length(estimators),labels=FALSE)
    text(x = 1:length(estimators), y = par("usr")[3] - 0.1, labels = estimators, srt = angle, adj = 1, xpd = TRUE, cex=axis.font.size)
  }
  abline(h=1,col="red")
}
for (dgp in res.list[4:6])
{
  ind <- ind+1
  out <- matrix(0,100,8)
  colnames(out) <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate","reguly")
  for (sample in 1:100)
  {
    file <- paste0("sample_",sample,".csv")
    methods <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate")
    file.name <- paste0(paste("Results/RMSE/",dgp,"/",methods,sep=""),"_",file)
    for (j in 1:length(methods))
    {
      data <- readRDS(paste0("Data/dgp_",dgp,"_",gsub("\\.csv","\\.rds",file)))
      cate <- data$cate[data$x>-0.1 & data$x<0.1]
      ##
      out[sample,methods[j]] <- as.numeric(read.table(file.name[j]))
    }
    file.reguly <- paste0("RD-Tree/Results/",dgp,"/",file)
    out.reguly <- as.matrix(read.table(file.reguly))
    out.reguly <- sqrt(mean((out.reguly[,1]-cate)^2))
    out[sample,"reguly"] <- out.reguly
    out[sample,] <- out[sample,]/out[sample,"ate"]
  }
  out <- out[,colnames(out) != "ate"]
  colnames(out) <- c("BARDDT","T-BART","S-BART","OLS","HS","Ridge","RD-Tree")
  if (ind==5)
  {
    boxplot(out[,-3], main = "",frame=FALSE,ylim=c(0,3),xaxt="n")
    mtext(bquote("Low signal, non-separable "~mu),line=1,font=2)
    axis(1,1:length(estimators),labels=FALSE)
    text(x = 1:length(estimators), y = par("usr")[3] - 0.2, labels = estimators, srt = angle, adj = 1, xpd = TRUE, cex=axis.font.size)
  } else
  {
    boxplot(out[,-3], main = "",frame=FALSE,ylim=c(0,3),xaxt="n")
    axis(1,1:length(estimators),labels=FALSE)
    text(x = 1:length(estimators), y = par("usr")[3] - 0.2, labels = estimators, srt = angle, adj = 1, xpd = TRUE, cex=axis.font.size)
  }
  abline(h=1,col="red")
}
dev.off()
# Figure 6---------------------------------------------------------------------
sample <- 4
## Easy fits
dgp <- res.list[1]
file <- paste0("sample_",sample,".csv")
methods <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate")
file.name <- paste0(paste("Results/Fits/",dgp,"/cate/",methods,sep=""),"_",file)
file.reguly <- paste0("RD-Tree/Results/",dgp,"/",file)
reguly <- read.table(file.reguly)
###
mat <- merge(read.table(file.name[1]),read.table(file.name[2]),by=1:2)
colnames(mat) <- c("W","CATE","BARDDT","T-BART")
mat <- merge(mat,read.table(file.name[3]),by=1:2)
colnames(mat)[5] <- "S-BART"
mat <- merge(mat,read.table(file.name[4]),by=1:2)
colnames(mat)[6] <- "Polynomial"
mat <- merge(mat,read.table(file.name[5]),by=1:2)
colnames(mat)[7] <- "Horseshoe"
mat <- merge(mat,read.table(file.name[6]),by=1:2)
colnames(mat)[8] <- "Ridge"
mat <- cbind(mat,reguly)
colnames(mat)[9] <- "RD-Tree"
saveRDS(mat,"easy_fits.rds")
## Hard fits
dgp <- res.list[4]
file <- paste0("sample_",sample,".csv")
methods <- c("barddt","tbart","sbart","polynomial","horseshoe","ridge","ate")
file.name <- paste0(paste("Results/Fits/",dgp,"/cate/",methods,sep=""),"_",file)
file.reguly <- paste0("RD-Tree/Results/",dgp,"/",file)
reguly <- read.table(file.reguly)
####
mat <- merge(read.table(file.name[1]),read.table(file.name[2]),by=1:2)
colnames(mat) <- c("W","CATE","BARDDT","T-BART")
mat <- merge(mat,read.table(file.name[3]),by=1:2)
colnames(mat)[5] <- "S-BART"
mat <- merge(mat,read.table(file.name[4]),by=1:2)
colnames(mat)[6] <- "Polynomial"
mat <- merge(mat,read.table(file.name[5]),by=1:2)
colnames(mat)[7] <- "Horseshoe"
mat <- merge(mat,read.table(file.name[6]),by=1:2)
colnames(mat)[8] <- "Ridge"
mat <- cbind(mat,reguly)
colnames(mat)[9] <- "RD-Tree"
saveRDS(mat,"hard_fits.rds")
## Plot setup
colors <- c("red","gold1","magenta3","dodgerblue","maroon","green3","darkgreen","grey")
bgcolors <- colors
lincol <- adjustcolor('red',alpha=0.75)
colvec <- c(lincol,rep("black",7))
cex = 1.25
pch <- 21
lwd <- 1.25
methods <- c("BARDDT","T-BART","S-BART","Horseshoe","RD-Tree")
pdf("Results/Figures/fits.pdf",width=7,height=10)
layout(matrix(c(1:4,9,5:8,9),ncol=2),heights = c(2,2,2,2,1))
par(mar=c(4,4,0.1,0.1))
## Easy setup
mat <- readRDS("easy_fits.rds")
ylim <- c(min(mat[,methods])-0.1,max(mat[,methods])+0.1)
mat <- mat[order(mat[,1]),]
###
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex = cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","p","n","n","n","n","n"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col= colvec, bg=colors, cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","p","n","n","n","n"),lty=1,lwd=lwd)
# matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
#         ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","p","n","n","n"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","n","p","n","n"),lty=1,lwd=lwd)
# matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
#         ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","n","n","p","n"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","n","n","n","p"),lty=1,lwd=lwd)
## Hard setup
mat <- readRDS("hard_fits.rds")
ylim <- c(min(mat[,methods])-0.1,max(mat[,methods])+0.1)
mat <- mat[order(mat[,1]),]
###
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex = cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","p","n","n","n","n","n"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col= colvec, bg=colors, cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","p","n","n","n","n"),lty=1,lwd=lwd)
# matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
#         ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","p","n","n","n"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","n","p","n","n"),lty=1,lwd=lwd)
# matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
#         ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","n","n","p","n"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("l","p","n","n","n","n","n","p"),lty=1,lwd=lwd)
par(mar=c(0.1,4,0.05,0.1))
plot.new()
legend("center",legend=methods,col="black",
       pt.bg=colors[colnames(mat)[-1] %in% methods],pch=21,cex=1.25,title="Estimator",title.font = 2,ncol=5)
dev.off()
