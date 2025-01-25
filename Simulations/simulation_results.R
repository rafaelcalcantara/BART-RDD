options(warn=2)
## List DGP per parameter values
dgp <- gsub("dgp_|.rds","",list.files("Data"))
## Results for each DGP
rmse.dgp <- matrix(0,length(dgp),4,dimnames = list(dgp,c("BARDDT","T-BART","S-BART","Polynomial")))
for (i in dgp)
{
  ## Read data
  data <- readRDS(paste0("Data/dgp_",i,".rds"))
  ate <- data$ate
  ## Keep CATE only for test sample
  test <- apply(data$x,2, function(xs) c-Owidth <= xs & xs <= c+Owidth)
  s <- ncol(test)
  cate <- lapply(1:s, function(rep) data$cate[test[,rep],rep])
  ## Load results
  files <- list.files("Results") ## get result files names
  methods <- gsub("_.*","",files) ## separate method name in file name
  ind <- grep(paste0(i,".rds"),files) ## pick out files relative to current DGP
  names(ind) <- methods[ind] ## so results are named in result list
  files <- paste0("Results/",files)
  results <- sapply(ind, function(i) readRDS(files[i]), simplify = FALSE, USE.NAMES = TRUE)
  ## Calculate RMSE per sample for each method
  rmse.fun <- function(fit,true) sqrt(mean((fit-true)^2))
  rmse <- data.frame(BARDDT = sapply(1:s, function(sample) rmse.fun(rowMeans(results$barddt[[sample]]),cate[[sample]])),
                     `T-BART` = sapply(1:s, function(sample) rmse.fun(rowMeans(results$tbart[[sample]]),cate[[sample]])),
                     `S-BART` = sapply(1:s, function(sample) rmse.fun(rowMeans(results$sbart[[sample]]),cate[[sample]])),
                     Polynomial = sapply(1:s, function(sample) rmse.fun(results$polynomial[[sample]],cate[[sample]])),
                     ATE = sapply(1:s, function(sample) rmse.fun(ate,cate[[sample]])))
  ## Visualize results
  # plot(rmse$BARDDT/rmse$ATE,rmse$T.BART/rmse$ATE, col="navy", pch=19, bty="n",
  #      xlab="BARDDT", ylab="T-BART", ylim=c(0,max(rmse$BARDDT/rmse$ATE,rmse$T.BART/rmse$ATE)))
  # abline(a=0,b=1,col="red")
  # plot(rmse$BARDDT/rmse$ATE,rmse$S.BART/rmse$ATE, col="navy", pch=19, bty="n",
  #      xlab="BARDDT", ylab="S-BART", ylim=c(0,max(rmse$BARDDT/rmse$ATE,rmse$S.BART/rmse$ATE)))
  # abline(a=0,b=1,col="red")
  # plot(rmse$BARDDT/rmse$ATE,rmse$Polynomial/rmse$ATE, col="navy", pch=19, bty="n",
  #      xlab="BARDDT", ylab="Polynomial", ylim=c(0,max(rmse$BARDDT/rmse$ATE,rmse$Polynomial/rmse$ATE)))
  # abline(a=0,b=1,col="red")
  ### Plotting tau fits averaged over w[,1]
  aggregate.per.w <- function(w,object,obj.name)
  {
    out <- lapply(1:s, function(sample) cbind(w[test[,sample],1],object[[sample]]))
    out <- do.call("rbind",out)
    out <- aggregate(out[,2],by=list(w=out[,1]),mean)
    colnames(out)[2] <- obj.name
    return(out)
  }
  aggregates <- cbind(rowSums(data$w[,1]*test/rowSums(test)),rowSums(data$cate*test/rowSums(test)))
  aggregates <- aggregates[complete.cases(aggregates[,1]),]
  colnames(aggregates) <- c("w","CATE")
  aggregates <- merge(aggregates,aggregate.per.w(data$w,sapply(results$barddt,rowMeans),"BARDDT"))
  aggregates <- merge(aggregates,aggregate.per.w(data$w,sapply(results$tbart,rowMeans),"T-BART"))
  aggregates <- merge(aggregates,aggregate.per.w(data$w,sapply(results$sbart,rowMeans),"S-BART"))
  aggregates <- merge(aggregates,aggregate.per.w(data$w,results$polynomial,"Polynomial"))
  ####
  pdf(paste0("plot",which(dgp==i),".pdf"),width = 10, height = 20)
  par(mfrow=c(2,2))
  matplot(aggregates[,1],aggregates[,c(2,3,4)],col=c("black","orange","maroon"),
          pch=19,bty="n",xlab=bquote(W[1]),ylab="CATE",main=gsub("_"," ",i))
  matplot(aggregates[,1],aggregates[,c(2,3,5)],col=c("black","orange","red"),
          pch=19,bty="n",xlab=bquote(W[1]),ylab="CATE",main=gsub("_"," ",i))
  matplot(aggregates[,1],aggregates[,c(2,3,6)],col=c("black","orange","pink"),
          pch=19,bty="n",xlab=bquote(W[1]),ylab="CATE",main=gsub("_"," ",i))
  dev.off()
  ## Summarizing results
  rmse.dgp[i,] <- colMeans(rmse/rmse$ATE)[1:(ncol(rmse)-1)]
}
## Best method per DGP
best <- apply(rmse.dgp,1, function(rmse) which(rmse==min(rmse)))
print(which(best!=1))
print(table(best))
## Add parameters to RMSE matrix
parameters <- do.call("rbind",strsplit(dgp,"_"))
param.names <- parameters[1,seq(1,14,2)]
parameters <- apply(parameters[,seq(2,14,2)],2,as.numeric)
colnames(parameters) <- param.names
rmse.dgp <- cbind(parameters,rmse.dgp)
rownames(rmse.dgp) <- NULL
print(round(rmse.dgp,2))