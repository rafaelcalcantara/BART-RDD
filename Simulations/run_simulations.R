cluster <- TRUE
n <- 500
s <- 10
# Create necessary folders if non-existent
## Data
if (!dir.exists("Data")) dir.create("Data") ## Create data folder, if non-existent
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
## Logs
if ("Logs" %in% list.files() == FALSE) dir.create("Logs")
## Results and subfolders
if ("Results" %in% list.files() == FALSE) dir.create("Results")
if ("RMSE" %in% list.files("Results") == FALSE) dir.create("Results/RMSE")
if ("Fits" %in% list.files("Results") == FALSE) dir.create("Results/Fits")
if ("Screenshots" %in% list.files("Results") == FALSE) dir.create("Results/Screenshots")
## Times
if ("Time" %in% list.files() == FALSE) dir.create("Time")
# Generate scripts for each DGP
source("generate_master_scripts.R")
# Run all regressions
for (i in 1:1) source(paste0("simulation_master",i,".R"))