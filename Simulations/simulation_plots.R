# Owidth <- 0.1
# c <- 0
## Boxplots per DGP
if ("Figures" %in% list.files("Results") == FALSE) dir.create("Results/Figures")
pdf(paste0("Results/Figures/boxplots_.pdf"))
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
  out <- out[,1:4]
  colnames(out) <- c("BARDDT","T-BART","S-BART","Polynomial")
  ####
  # pdf(paste0("Results/Figures/boxplots_",dgp,".pdf"))
  boxplot(out, main = dgp,frame=FALSE)
  abline(h=1,col="red")
  # dev.off()
}
dev.off()
## Average fits per W
pdf(paste0("Results/Figures/fits_.pdf"))
for (i in 1:length(list.files("Results/Fits")))
{
  dgp <- list.files("Results/Fits")[i]
  file.name <- paste0("Results/Fits/",dgp,"/sample_1.csv")
  mat <- read.table(file.name)
  colnames(mat) <- c("W","BARDDT","T-BART","S-BART","Polynomial")
  ###
  # pdf(paste0("Results/Figures/fits_",dgp,".pdf"))
  par(mfrow=c(2,2))
  matplot(mat[,1],mat[,c(2,3,4)],pch=19,cex=0.6,col=c("black","orange","maroon"),
          ylab="CATE",xlab=bquote(W[1]),bty="n")
  matplot(mat[,1],mat[,c(2,3,5)],pch=19,cex=0.6,col=c("black","orange","red"),
          ylab="CATE",xlab=bquote(W[1]),bty="n")
  matplot(mat[,1],mat[,c(2,3,6)],pch=19,cex=0.6,col=c("black","orange","pink"),
          ylab="CATE",xlab=bquote(W[1]),bty="n")
  plot.new()
  legend("center",legend=c("BARDDT","T-BART","S-BART","Polynomial"),
         col=c("orange","maroon","red","pink"),pch=19,cex=0.75,title="Estimator",title.font = 2)
  mtext(paste(strsplit(dgp,"_")[[1]],collapse=" "), side = 3, line = - 2, outer = TRUE, font = 2)
  # dev.off()
}
dev.off()
## Average for each DGP
rmse.mean <- matrix(0,14,11)
colnames(rmse.mean) <- c("k1","k2","k3","k4","k5","p","rho",
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
  rmse.mean[i,] <- as.numeric(c(strsplit(dgp,"_")[[1]][(1:7)*2],round(colMeans(out)[1:4],3)))
}
rmse.mean <- as.data.frame(rmse.mean)
subset(rmse.mean,k3==1)
subset(rmse.mean,k3==0)