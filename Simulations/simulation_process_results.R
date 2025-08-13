# Function to read data, results, calculate and store RMSE
calc.rmse <- function(sample,ate,test)
{
  rmse <- function(fit,true) sqrt(mean((true-fit)^2))
  cate <- cate[test]
  pred.ate <- rep(ate,length(cate))
  out <- rmse(pred.ate,cate)
  write.table(out,paste0("Results/RMSE/",dgp,"/ate_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  if ("barddt" %in% models)
  {
    pred.barddt <- get("pred.barddt",envir=parent.frame())
    barddt <- pred.barddt$post
    out <- rmse(rowMeans(barddt),cate)
    write.table(out,paste0("Results/RMSE/",dgp,"/barddt_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("tbart" %in% models)
  {
    pred.tbart <- get("pred.tbart",envir=parent.frame())
    tbart <- pred.tbart$post
    out <- rmse(rowMeans(tbart),cate)
    write.table(out,paste0("Results/RMSE/",dgp,"/tbart_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("sbart" %in% models)
  {
    pred.sbart <- get("pred.sbart",envir=parent.frame())
    sbart <- pred.sbart$post
    out <- rmse(rowMeans(sbart),cate)
    write.table(out,paste0("Results/RMSE/",dgp,"/sbart_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("polynomial" %in% models)
  {
    pred.polynomial <- get("pred.polynomial",envir=parent.frame())
    polynomial <- pred.polynomial$cate
    out <- rmse(polynomial,cate)
    write.table(out,paste0("Results/RMSE/",dgp,"/polynomial_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("horseshoe" %in% models)
  {
    pred.horseshoe <- get("pred.horseshoe",envir=parent.frame())
    horseshoe <- pred.horseshoe$post
    out <- rmse(rowMeans(horseshoe),cate)
    write.table(out,paste0("Results/RMSE/",dgp,"/horseshoe_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("ridge" %in% models)
  {
    pred.ridge <- get("pred.ridge",envir=parent.frame())
    ridge <- pred.ridge$post
    out <- rmse(rowMeans(ridge),cate)
    write.table(out,paste0("Results/RMSE/",dgp,"/ridge_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
}
# Function to generate screenshots (important to keep track while running on cluster)
screenshot <- function(s0,s1)
{
  out <- matrix(0,length(s0:s1),5)
  for (i in 1:length(s0:s1))
  {
    out[i,] <- as.matrix(read.table(paste0("Results/RMSE/",dgp,"/sample_",s0+(i-1),".csv")))
  }
  out <- sweep(out,1,out[,ncol(out)],"/")
  txt <- paste(c(rep("Mean",5),rep("SD",5)),rep(c("BARDDT","T-BART","S-BART","Polynomial","Horseshoe","Ridge","ATE")),c(round(colMeans(out),2),round(apply(out,2,sd),2)))
  writeLines(txt,paste0("Results/Screenshots/",dgp,"/sample_",s0,"_",s1,".txt"))
}
# Function to collect point estimates (posterior means and polynomial fit)
point.est <- function(sample,ate,test)
{
  ## CATE
  cate <- cate[test]
  w <- w[test,1]
  temp <- cbind(W=w,CATE=cate)
  pred.ate <- rep(ate,length(cate))
  out <- cbind(temp,ATE=pred.ate)
  write.table(out,paste0("Results/Fits/",dgp,"/cate/ate_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  if ("barddt" %in% models)
  {
    pred.barddt <- get("pred.barddt",envir=parent.frame())
    barddt <- pred.barddt$post
    out <- cbind(temp,BARDDT=rowMeans(barddt))
    write.table(out,paste0("Results/Fits/",dgp,"/cate/barddt_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("tbart" %in% models)
  {
    pred.tbart <- get("pred.tbart",envir=parent.frame())
    tbart <- pred.tbart$post
    out <- cbind(temp,`T-BART`=rowMeans(tbart))
    write.table(out,paste0("Results/Fits/",dgp,"/cate/tbart_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("sbart" %in% models)
  {
    pred.sbart <- get("pred.sbart",envir=parent.frame())
    sbart <- pred.sbart$post
    out <- cbind(temp,`S-BART`=rowMeans(sbart))
    write.table(out,paste0("Results/Fits/",dgp,"/cate/sbart_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("polynomial" %in% models)
  {
    pred.polynomial <- get("pred.polynomial",envir=parent.frame())
    polynomial <- pred.polynomial$cate
    out <- cbind(temp,polynomial=polynomial)
    write.table(out,paste0("Results/Fits/",dgp,"/cate/polynomial_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("horseshoe" %in% models)
  {
    pred.horseshoe <- get("pred.horseshoe",envir=parent.frame())
    horseshoe <- pred.horseshoe$post
    out <- cbind(temp,Horseshoe=rowMeans(horseshoe))
    write.table(out,paste0("Results/Fits/",dgp,"/cate/horseshoe_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  if ("ridge" %in% models)
  {
    pred.ridge <- get("pred.ridge",envir=parent.frame())
    ridge <- pred.ridge$post
    out <- cbind(temp,Ridge=rowMeans(ridge))
    write.table(out,paste0("Results/Fits/",dgp,"/cate/ridge_sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
  }
  # ## Yhat
  # barddt <- pred.barddt$yhat
  # tbart1 <- pred.tbart$yhat1
  # tbart0 <- pred.tbart$yhat0
  # sbart <- pred.sbart$yhat
  # polynomial <- pred.polynomial$yhat
  # horseshoe <- pred.horseshoe$yhat
  # ridge <- pred.ridge$yhat
  # ### Make polynomial match dimensions
  # poly.sample <- c-bw<=x & x<=c+h
  # ind <- 1
  # poly <- rep(NA,length(poly.sample))
  # for (i in 1:length(poly.sample))
  # {
  #   if (isTRUE(poly.sample[i]))
  #   {
  #     poly[i] <- polynomial[ind]
  #     ind <- ind+1
  #   }
  # }
  # ## Make T-BART match dimensions
  # ind0 <- ind1 <- 1
  # tbart <- rep(NA,length(z))
  # tbart1 <- rowMeans(tbart1)
  # tbart0 <- rowMeans(tbart0)
  # for (i in 1:length(z))
  # {
  #   if (z[i]==1)
  #   {
  #     tbart[i] <- tbart1[ind1]
  #     ind1 <- ind1 + 1
  #   } else
  #   {
  #     tbart[i] <- tbart0[ind0]
  #     ind0 <- ind0 + 1
  #   }
  # }
  # ##
  # out <- cbind(X=x,Y=y,BARDDT=rowMeans(barddt),`T-BART`=tbart,`S-BART`=rowMeans(sbart),Polynomial=poly)
  # write.table(out,paste0("Results/Fits/",dgp,"/yhat/sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
}
