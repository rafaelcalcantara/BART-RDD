setwd("~/../Git/BART-RDD/")
## Get simulation script names
scripts <- list.files("R")[grep("win",list.files("R"))]
scripts <- scripts[-grep("master|data",scripts)]
scripts <- paste0("R/",scripts)
for (i in scripts)
{
  print(paste0("Script: ",i))
  time <- system.time({
    source(i)
  })
  print(paste0("Elapsed time: ",time))
}