## Master script for simulations with DGPs from the other papers
rm(list=ls())
n <- 500
p <- 10
c <- 0
s <- 1 # samples
if (Sys.info()["sysname"] == "Linux")
{
    ## XBCF-RDD simulations
    source("R/simulations_xbcf.R")
    rm(list=ls())
    gc()
    ## KR simulations
    source("R/simulations_kr.R")
    rm(list=ls())
    gc()
    ## FH simulations
    source("R/simulations_fh.R")
    rm(list=ls())
    gc()
} else
{
    ## CGS simulations
    source("R/simulations_cgs.R")
    rm(list=ls())
    gc()
}
