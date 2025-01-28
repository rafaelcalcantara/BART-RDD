# Function to read data, results, calculate and store RMSE
calc.rmse <- function(sample)
{
  data <- readRDS(paste0("Data/",dgp,"/sample_",sample,".rds"))
  test <- c-Owidth <= data$x & data$x <= c+Owidth
  cate <- data$cate[test]
  barddt <- read.table(paste0("Results/",dgp,"/barddt_sample_",sample,".csv"),sep=",")
  tbart <- read.table(paste0("Results/",dgp,"/tbart_sample_",sample,".csv"),sep=",")
  sbart <- read.table(paste0("Results/",dgp,"/sbart_sample_",sample,".csv"),sep=",")
  polynomial <- read.table(paste0("Results/",dgp,"/polynomial_sample_",sample,".csv"),sep=",")
  rmse <- function(fit,true) sqrt(mean((true-fit)^2))
  out <- cbind(rmse(rowMeans(barddt),cate),rmse(rowMeans(tbart),cate),
               rmse(rowMeans(sbart),cate),rmse(polynomial[,1],cate),rmse(data$ate,cate))
  write.table(out,paste0("Results/RMSE/",dgp,"/sample_",sample,".csv"), row.names = FALSE, col.names = FALSE)
}
# Function to generate screenshots (important to keep track while running on cluster)
screenshot <- function(s0,s1)
{
  out <- matrix(0,length(s0:s1),5)
  for (i in s0:s1)
  {
    out[i,] <- as.matrix(read.table(paste0("Results/RMSE/",dgp,"/sample_",i,".csv")))
  }
  out <- sweep(out,1,out[,5],"/")
  txt <- paste(c(rep("Mean",5),rep("SD",5)),rep(c("BARDDT","T-BART","S-BART","Polynomial","ATE")),c(round(colMeans(out),2),round(apply(out,2,sd),2)))
  writeLines(txt,paste0("Results/Screenshots/",dgp,"/sample_",s0,"_",s1,".txt"))
}
# Function to collect point estimates (posterior means and polynomial fit)
point.est <- function(sample)
{
  data <- readRDS(paste0("Data/",dgp,"/sample_",sample,".rds"))
  test <- c-Owidth <= data$x & data$x <= c+Owidth
  cate <- data$cate[test]
  w <- data$w[test,1]
  barddt <- read.table(paste0("Results/",dgp,"/barddt_sample_",sample,".csv"),sep=",")
  tbart <- read.table(paste0("Results/",dgp,"/tbart_sample_",sample,".csv"),sep=",")
  sbart <- read.table(paste0("Results/",dgp,"/sbart_sample_",sample,".csv"),sep=",")
  polynomial <- read.table(paste0("Results/",dgp,"/polynomial_sample_",sample,".csv"),sep=",")
  out <- cbind(W=w,BARDDT=rowMeans(barddt),`T-BART`=rowMeans(tbart),`S-BART`=rowMeans(sbart),Polynomial=polynomial)
  write.table(out,paste0("Results/Fits/",dgp,"/sample_",sample,".csv"), row.names = FALSE, col.names = TRUE)
}
# Create RMSE subfolder for current DGP if it doesn't exist
if ("RMSE" %in% list.files("Results") == FALSE) dir.create("Results/RMSE")
## Running locally
# cl <- makeCluster(no_cores,type="SOCK")
# registerDoParallel(cl)
# clusterExport(cl,varlist=ls())
# time <- system.time({
#   out <- parLapply(cl,1:s,calc.rmse)
# })
# stopCluster(cl)
## Running from cluster
### Change values to whatever batch of samples to run
s0 <- 1
s1 <- 200
c <- 0
Owidth <- 0.1
# for (dgp in list.files("Results")[list.files("Results") %in% c("RMSE","Screenshots","Fits")==F])
# {
#   if (dgp %in% list.files("Results/RMSE") == FALSE) dir.create(paste0("Results/RMSE/",dgp))
#   for (i in s0:s1)
#   {
#     calc.rmse(i)
#   }
# }
# # Create screenshots subfolder for current DGP if it doesn't exist
# if ("Screenshots" %in% list.files("Results") == FALSE) dir.create("Results/Screenshots")
# for (dgp in list.files("Results")[list.files("Results") %in% c("RMSE","Screenshots","Fits")==F])
# {
#   if (dgp %in% list.files("Results/Screenshots") == FALSE) dir.create(paste0("Results/Screenshots/",dgp))
#   screenshot(s0,s1)
# }
# Create Fits subfolder if it doesn't exist
if ("Fits" %in% list.files("Results") == FALSE) dir.create("Results/Fits")
for (dgp in list.files("Results")[list.files("Results") %in% c("RMSE","Screenshots","Fits")==F])
{
  if (dgp %in% list.files("Results/Fits") == FALSE) dir.create(paste0("Results/Fits/",dgp))
  point.est(1)
}