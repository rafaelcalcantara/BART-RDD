## Create/clean data folder
if (!dir.exists("Data")) dir.create("Data") ## Create data folder, if non-existent
if (length(list.files("Data"))!=0) ## Clean up folder
{
  files <- paste0("Data/",list.files("Data"))
  for (i in files) file.remove(i)
}
## Create/clean results folder
if (!dir.exists("Results")) dir.create("Results") ## Create results folder, if non-existent
if (length(list.files("Results"))!=0) ## Clean up folder
{
  files <- paste0("Results/",list.files("Results"))
  for (i in files) file.remove(i)
}