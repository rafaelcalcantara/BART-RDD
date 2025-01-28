Owidth <- 0.1
c <- 0
## Boxplots per DGP
if ("Figures" %in% list.files("Results") == FALSE) dir.create("Results/Figures")
pdf("Results/Figures/boxplots.pdf")
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
  boxplot(out, main = paste(strsplit(dgp,"_")[[1]],collapse=" "),frame=FALSE)
  abline(h=1,col="red")
}
dev.off()
## Average fits per W
pdf("Results/Figures/fits.pdf",width = 9, height = 5)
for (i in 1:length(list.files("Results/Fits")))
{
  dgp <- list.files("Results/Fits")[i]
  file.name <- paste0("Results/Fits/",dgp,"/sample_1.csv")
  mat <- read.table(file.name)
  colnames(mat) <- c("W","BARDDT","T-BART","S-BART","Polynomial")
  par(mfrow=c(1,3))
  matplot(mat[,1],mat[,2:3],pch=19,col=c("orange","maroon"),ylab="CATE",xlab=bquote(W[1]))
}
dev.off()