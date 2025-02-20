## Check best estimator for all 128 DGPs
rmse <- matrix(0,128,11)
colnames(rmse) <- c("k1","k2","k3","k4","k5","p","rho",
                    "BARDDT","T-BART","S-BART","Polynomial")
for (i in 1:length(list.files("Results/RMSE")))
{
  dgp <- list.files("Results/RMSE")[i]
  out <- matrix(0,length(list.files(paste0("Results/RMSE/",dgp))),5)
  for (sample in 1:length(list.files(paste0("Results/RMSE/",dgp))))
  {
    file.name <- list.files(paste0("Results/RMSE/",dgp))[sample]
    file.name <- paste0("Results/RMSE/",dgp,"/",file.name)
    out[sample,] <- as.matrix(read.table(file.name))
    out[sample,] <- out[sample,]/out[sample,5]
  }
  rmse[i,] <- as.numeric(c(strsplit(dgp,"_")[[1]][(1:7)*2],round(colMeans(out)[1:4],3)))
}
table(apply(rmse[,8:11],1,function(i) which(i==min(i))))
rmse[apply(rmse[,8:11],1,function(i) which(i==min(i)))!=1,]
## Average for each DGP
res.list <- c("k1_1_k2_1_k3_0_k4_0.1_k5_0_p_2_rho_0.5",
              "k1_1_k2_1_k3_0_k4_0.1_k5_0_p_4_rho_0",
              "k1_1_k2_1_k3_0_k4_0.1_k5_1_p_2_rho_0",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_0_p_4_rho_0.5",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_1_p_2_rho_0.5",
              "k1_5_k2_0.25_k3_1_k4_0.5_k5_1_p_4_rho_0")
rmse.mean <- matrix(0,6,11)
colnames(rmse.mean) <- c("k1","k2","k3","k4","k5","p","rho",
                         "BARDDT","T-BART","S-BART","Polynomial")
for (i in 1:length(res.list))
{
  dgp <- res.list[i]
  out <- matrix(0,length(list.files(paste0("Results/RMSE/",dgp))),5)
  for (sample in 1:length(list.files(paste0("Results/RMSE/",dgp))))
  {
    file.name <- list.files(paste0("Results/RMSE/",dgp))[sample]
    file.name <- paste0("Results/RMSE/",dgp,"/",file.name)
    out[sample,] <- as.matrix(read.table(file.name))
    out[sample,] <- out[sample,]/out[sample,5]
  }
  rmse.mean[i,] <- as.numeric(c(strsplit(dgp,"_")[[1]][(1:7)*2],round(colMeans(out)[1:4],3)))
}
rmse.mean <- as.data.frame(rmse.mean)
subset(rmse.mean,k3==1)
subset(rmse.mean,k3==0)
rmse.mean
apply(rmse.mean[,8:11],1,function(i) which(i==min(i)))
table(apply(rmse.mean[,8:11],1,function(i) which(i==min(i))))
print(xtable::xtable(subset(rmse.mean,(k1==1 & k2==1 & k3==0 & k4==0.1) | (k1==5 & k2==0.25 & k3==1 & k4==0.5)),
                     caption="Average RMSE per DGP, also divided here by the RMSE of the naive ATE estimator",
                     label="tab:rmse.average",
                     digits=2,
                     align=rep("c",1+ncol(rmse.mean))),
      include.rownames=FALSE)
## Boxplot grid
ind <- 0
pdf(paste0("Results/Figures/boxplots.pdf"))
par(mfrow=c(2,3),mai=c(0.7,0.4,0.4,0.1))
for (dgp in res.list[1:3])
{
  ind <- ind+1
  out <- matrix(0,length(list.files(paste0("Results/RMSE/",dgp))),5)
  for (sample in 1:length(list.files(paste0("Results/RMSE/",dgp))))
  {
    file.name <- list.files(paste0("Results/RMSE/",dgp))[sample]
    file.name <- paste0("Results/RMSE/",dgp,"/",file.name)
    out[sample,] <- as.matrix(read.table(file.name))
    out[sample,] <- out[sample,]/out[sample,5]
  }
  out <- out[,1:4]
  colnames(out) <- c("BARDDT","T-BART","S-BART","Polynomial")
  if (ind==3)
  {
    boxplot(out[,-3], main = "",frame=FALSE,ylim=c(0,1.01),xaxt="n")
    axis(1,1:3,colnames(out)[-3],cex.axis=0.8)
  } else if (ind==2)
  {
    boxplot(out[,-3], main = "",frame=FALSE,xaxt="n")
    mtext(bquote("High signal, separable "~mu),line=1,font=2)
    axis(1,1:3,colnames(out)[-3],cex.axis=0.8)
  } else
  {
    boxplot(out[,-3], main = "",frame=FALSE,xaxt="n")
    axis(1,1:3,colnames(out)[-3],cex.axis=0.8)
  }
  abline(h=1,col="red")
}
for (dgp in res.list[4:6])
{
  ind <- ind+1
  out <- matrix(0,length(list.files(paste0("Results/RMSE/",dgp))),5)
  for (sample in 1:length(list.files(paste0("Results/RMSE/",dgp))))
  {
    file.name <- list.files(paste0("Results/RMSE/",dgp))[sample]
    file.name <- paste0("Results/RMSE/",dgp,"/",file.name)
    out[sample,] <- as.matrix(read.table(file.name))
    out[sample,] <- out[sample,]/out[sample,5]
  }
  out <- out[,1:4]
  colnames(out) <- c("BARDDT","T-BART","S-BART","Polynomial")
  if (ind==5)
  {
    boxplot(out[,-3], main = "",frame=FALSE,xaxt="n")
    mtext(bquote("Low signal, non-separable "~mu),line=1,font=2)
    axis(1,1:3,colnames(out)[-3],cex.axis=0.8)
  } else
  {
    boxplot(out[,-3], main = "",frame=FALSE,xaxt="n")
    axis(1,1:3,colnames(out)[-3],cex.axis=0.8)
  }
  abline(h=1,col="red")
}
dev.off()
##
for (dgp in res.list[1:3])
{
  if (strsplit(dgp,"_")[[1]][2]!="1" |strsplit(dgp,"_")[[1]][4]!="1" |
      strsplit(dgp,"_")[[1]][6]!="0" | strsplit(dgp,"_")[[1]][8]!="0.1" |
      strsplit(dgp,"_")[[1]][10]!="0" | strsplit(dgp,"_")[[1]][12]!="2" |
      strsplit(dgp,"_")[[1]][14]!="0.5") next
  file.name <- paste0("Results/Fits/",dgp,"/cate/sample_1.csv")
  mat <- read.table(file.name)
  colnames(mat) <- c("W","CATE","BARDDT","T-BART","S-BART","Polynomial")
  saveRDS(mat,"easy_fits.rds")
}
for (dgp in res.list[4:6])
{
  if (strsplit(dgp,"_")[[1]][2]!="5" |strsplit(dgp,"_")[[1]][4]!="0.25" |
      strsplit(dgp,"_")[[1]][6]!="1" | strsplit(dgp,"_")[[1]][8]!="0.5" |
      strsplit(dgp,"_")[[1]][10]!="0" | strsplit(dgp,"_")[[1]][12]!="4" |
      strsplit(dgp,"_")[[1]][14]!="0.5") next
  file.name <- paste0("Results/Fits/",dgp,"/cate/sample_8.csv")
  mat <- read.table(file.name)
  colnames(mat) <- c("W","CATE","BARDDT","T-BART","S-BART","Polynomial")
  saveRDS(mat,"hard_fits.rds")
}
####
colors <- c("magenta3","dodgerblue", "green3","gold1")
colors[1] <- grDevices::adjustcolor(colors[1],alpha=1)

bgcolors <- colors

lincol <- adjustcolor('red',alpha=0.75)

colvec <- c(rep("black",3),"black",lincol)
cex = 1.75
pch <- c(rep(21,3),rep(21,2))
lwd <- c(rep(1,4),2)
pdf("Results/Figures/fits.pdf",width=12,height=10)
layout(matrix(c(1,4,7,2,5,7,3,6,7),ncol=3),heights = c(2,2,1))
## Easy setup
mat <- readRDS("easy_fits.rds")
ylim <- c(min(mat[,-1])-0.1,max(mat[,-1])+0.1)
mat <- mat[order(mat[,1]),]
mat <- mat[,c(1,4,5,6,3,2)]
###
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex = cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("p","n","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col= colvec, bg=colors, cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","p","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","n","p","p","l"),lty=1,lwd=lwd)
## Hard setup
mat <- readRDS("hard_fits.rds")
mat <- mat[order(mat[,1]),]
mat <- mat[,c(1,4,5,6,3,2)]
###
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex = cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("p","n","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col= colvec, bg=colors, cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","p","n","p","l"),lty=1,lwd=lwd)
matplot(mat[,1],mat[,-1],pch=pch,col=colvec, bg=colors,cex=cex,
        ylab="CATE",xlab=bquote(W[1]),bty="n",ylim=ylim,type=c("n","n","p","p","l"),lty=1,lwd=lwd)
plot.new()
legend("center",legend=c("BARDDT","T-BART","S-BART","Polynomial"),
       col="black",pt.bg=colors[c(4,1,2,3)],pch=21,cex=1.25,title="Estimator",title.font = 2,ncol=4)
dev.off()