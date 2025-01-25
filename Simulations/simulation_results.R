## List DGP per parameter values
dgp <- gsub("dgp_|.rds","",list.files("Data"))
## Results for each DGP
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
  ind <- grep(i,files) ## pick out files relative to current DGP
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
  plot(rmse$BARDDT/rmse$ATE,rmse$T.BART/rmse$ATE, col="navy", pch=19, bty="n",
       xlab="BARDDT", ylab="T-BART", ylim=c(min()))
  abline(a=0,b=1,col="red")
}