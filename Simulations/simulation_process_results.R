# Function to read data, results, calculate and store RMSE
calc.rmse <- function(sample,ate,test,pred.barddt,pred.tbart,pred.sbart,pred.polynomial)
{
  cate <- cate[test]
  barddt <- pred.barddt$post
  tbart <- pred.tbart$post
  sbart <- pred.sbart$post
  polynomial <- pred.polynomial$cate
  rmse <- function(fit,true) sqrt(mean((true-fit)^2))
  out <- c(rmse(rowMeans(barddt),cate),rmse(rowMeans(tbart),cate),
           rmse(rowMeans(sbart),cate),rmse(polynomial,cate),rmse(ate,cate))
  write.table(out,paste0("Results/RMSE/",dgp,"/sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
}
# Function to generate screenshots (important to keep track while running on cluster)
screenshot <- function(s0,s1)
{
  out <- matrix(0,length(s0:s1),5)
  for (i in 1:length(s0:s1))
  {
    out[i,] <- as.matrix(read.table(paste0("Results/RMSE/",dgp,"/sample_",s0+(i-1),".csv")))
  }
  out <- sweep(out,1,out[,5],"/")
  txt <- paste(c(rep("Mean",5),rep("SD",5)),rep(c("BARDDT","T-BART","S-BART","Polynomial","ATE")),c(round(colMeans(out),2),round(apply(out,2,sd),2)))
  writeLines(txt,paste0("Results/Screenshots/",dgp,"/sample_",s0,"_",s1,".txt"))
}
# Function to collect point estimates (posterior means and polynomial fit)
point.est <- function(sample,test,h,c,pred.barddt,pred.tbart,pred.sbart,pred.polynomial)
{
  ## CATE
  cate <- cate[test]
  w <- w[test,1]
  barddt <- pred.barddt$post
  tbart <- pred.tbart$post
  sbart <- pred.sbart$post
  polynomial <- pred.polynomial$cate
  out <- cbind(W=w,CATE=cate,BARDDT=rowMeans(barddt),`T-BART`=rowMeans(tbart),`S-BART`=rowMeans(sbart),Polynomial=polynomial)
  write.table(out,paste0("Results/Fits/",dgp,"/cate/sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  ## Yhat
  barddt <- pred.barddt$yhat
  tbart1 <- pred.tbart$yhat1
  tbart0 <- pred.tbart$yhat0
  sbart <- pred.sbart$yhat
  polynomial <- pred.polynomial$yhat
  ### Make polynomial match dimensions
  poly.sample <- c-h<=x & x<=c+h
  ind <- 1
  poly <- rep(NA,length(poly.sample))
  for (i in 1:length(poly.sample))
  {
    if (isTRUE(poly.sample[i]))
    {
      poly[i] <- polynomial[ind]
      ind <- ind+1
    }
  }
  ## Make T-BART match dimensions
  ind0 <- ind1 <- 1
  tbart <- rep(NA,length(z))
  tbart1 <- rowMeans(tbart1)
  tbart0 <- rowMeans(tbart0)
  for (i in 1:length(z))
  {
    if (z[i]==1)
    {
      tbart[i] <- tbart1[ind1]
      ind1 <- ind1 + 1
    } else
    {
      tbart[i] <- tbart0[ind0]
      ind0 <- ind0 + 1
    }
  }
  ##
  out <- cbind(X=x,Y=y,BARDDT=rowMeans(barddt),`T-BART`=tbart,`S-BART`=rowMeans(sbart),Polynomial=poly)
  write.table(out,paste0("Results/Fits/",dgp,"/yhat/sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
}