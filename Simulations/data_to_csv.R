files <- list.files("Data")
files <- files[files %in% c("w_2.rds","w_4.rds") == FALSE]
folders <- gsub(".*dgp_(.+)_sample.*","\\1",files)
samples <- gsub(".*sample","sample",files)
samples <- gsub(".rds",".mat",samples)
## Create data folder for Reguly estimation
if (dir.exists("RD-Tree/Data")==FALSE)
{
  dir.create("RD-Tree/Data")
}
for (i in 1:length(files))
{
  dgp <- files[i]
  subfolder <- paste0("RD-Tree/Data/",folders[i])
  sample <- samples[i]
  if (dir.exists(subfolder)==FALSE)
  {
    dir.create(subfolder)
  }
  data <- readRDS(paste0("Data/",dgp))
  R.matlab::writeMat(paste0(subfolder,"/",sample),a=data)
}
